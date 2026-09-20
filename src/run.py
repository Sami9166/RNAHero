"""One-command GEO collection and safe handoff to the edgeR workflow."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path
from typing import Any

from dotenv import load_dotenv

from agent import run_critic, run_development_critic, run_summarizer_agent
from ncbi_mcp import fetch_geo_sample_metadata, fetch_raw_count_matrix, find_bulk_rnaseq, match_count_columns_to_samples
from summarizer import run_summarizer
from workflow import run_pipeline


AI_FIELD_PREFIXES = (
    "Sample_title",
    "Sample_source_name",
    "Sample_characteristics",
    "Sample_description",
    "Sample_treatment_protocol",
    "Sample_growth_protocol",
)
MAX_AI_FIELD_CHARS = 500
COHORT_BATCH_SIZE = 4
DEFAULT_MAX_LOOP_ROUNDS = 2
BASE_AGENT_CALLS = 6  # Orchestrator + three development critics + external critic + summarizer.
RETRYABLE_REASON_CODES = {
    "cohort_selection",
    "insufficient_internal_validation",
    "no_locked_candidates",
    "reproducibility_gap",
}


class ApiCallBudgetExceeded(RuntimeError):
    """Raised when an agent workflow would exceed its declared API budget."""


@dataclass
class ApiCallLedger:
    """Track Gemini/ADK calls and the proposal's loop budget formula."""

    candidate_count: int
    max_loop_rounds: int
    batch_size: int = COHORT_BATCH_SIZE
    events: list[dict[str, Any]] = field(default_factory=list)
    loop_rounds: int = 0
    budget: int = field(init=False)
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    def __post_init__(self) -> None:
        if self.max_loop_rounds < 0:
            raise ValueError("max_loop_rounds must be zero or greater")
        self.budget = math.ceil(max(0, self.candidate_count) / self.batch_size) + BASE_AGENT_CALLS

    @property
    def api_call_count(self) -> int:
        return sum(1 for event in self.events if event.get("kind") != "local")

    def reserve(self, kind: str, phase: str, loop_round: int = 0) -> None:
        with self._lock:
            if self.api_call_count >= self.budget:
                raise ApiCallBudgetExceeded(
                    f"API call budget exceeded ({self.api_call_count}/{self.budget}) before {phase}"
                )
            self.events.append({"kind": kind, "phase": phase, "loop_round": loop_round})

    def record_local(self, phase: str, loop_round: int = 0) -> None:
        with self._lock:
            self.events.append({"kind": "local", "phase": phase, "loop_round": loop_round})

    def open_loop(self, new_study_count: int) -> None:
        """Extend the budget by K_l+4 for one retry round."""
        if self.loop_rounds >= self.max_loop_rounds:
            raise ApiCallBudgetExceeded("max_loop_rounds reached")
        batches = max(1, math.ceil(max(1, new_study_count) / self.batch_size))
        with self._lock:
            self.budget += batches + 4
            self.loop_rounds += 1

    def snapshot(self) -> dict[str, Any]:
        return {
            "candidate_count": self.candidate_count,
            "batch_size": self.batch_size,
            "base_budget": math.ceil(max(0, self.candidate_count) / self.batch_size) + BASE_AGENT_CALLS,
            "max_loop_rounds": self.max_loop_rounds,
            "loop_rounds": self.loop_rounds,
            "budget": self.budget,
            "api_call_count": self.api_call_count,
            "events": list(self.events),
            "formula": "ceil(N/4) + 6 + sum(K_l + 4)",
        }


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


RANKING_FIELDS = [
    "rank", "gene_id", "discovery", "direction", "min_auc", "mean_auc", "min_sensitivity", "min_specificity", "passed",
    "validation_1_cohort", "validation_1_auc", "validation_1_sensitivity", "validation_1_specificity",
    "validation_2_cohort", "validation_2_auc", "validation_2_sensitivity", "validation_2_specificity",
]
EXTERNAL_VALIDATION_FIELDS = ["gene_id", "validation_gene_id", "discovery", "cohort", "direction", "auc", "sensitivity", "specificity", "cutoff"]


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RANKING_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_validation_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXTERNAL_VALIDATION_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _candidate_preview(study: dict[str, Any]) -> dict[str, str]:
    return {"gse_id": study["gse_id"], "title": study["title"]}


def _cohort_summary(study: dict[str, Any], prepared: dict[str, Any]) -> dict[str, Any]:
    assessment = study["cohort_assessment"]
    return {
        "gse_id": study["gse_id"],
        "role": prepared["role"],
        "title": study["title"],
        "sample_count": study["sample_count"],
        "raw_count_file": study["raw_count_files"][0]["filename"],
        "case_count": assessment["case_count"],
        "control_count": assessment["control_count"],
        "selection_reason": prepared.get("selection_reason", ""),
    }


def clear_output(output_dir: Path) -> dict[str, str]:
    """Delete only a run folder below the project directory."""
    root = Path.cwd().resolve()
    target = output_dir.resolve()
    if target == root or root not in target.parents:
        raise ValueError("only a subfolder of the RNAHero project can be cleared")
    if target.exists():
        shutil.rmtree(target)
    return {"stage": "cleared", "output": str(target)}


load_dotenv(Path(__file__).resolve().parent.parent / ".env")


def _retry_seconds(message: str) -> int | None:
    match = re.search(r"(?:retryDelay|retry in).{0,40}?(\d+(?:\.\d+)?)s", message, flags=re.IGNORECASE)
    return min(60, max(1, int(float(match.group(1))))) if match else None


def _ask_gemini(
    prompt: str,
    *,
    ledger: ApiCallLedger | None = None,
    phase: str = "curator",
    loop_round: int = 0,
) -> tuple[dict[str, Any] | None, str | None]:
    if not os.getenv("GEMINI_API_KEY"):
        return None, "GEMINI_API_KEY is not configured"
    if ledger is not None:
        try:
            ledger.reserve("gemini", phase, loop_round)
        except ApiCallBudgetExceeded as error:
            return None, str(error)
    from google import genai

    client = genai.Client(api_key=os.environ["GEMINI_API_KEY"])
    for attempt in range(3):
        try:
            response = client.models.generate_content(
                model=os.getenv("GEMINI_MODEL", "gemini-3.5-flash-lite"),
                contents=prompt,
                config={"response_mime_type": "application/json"},
            )
            return json.loads(response.text), None
        except Exception as error:  # The run manifest preserves the exact safe stop reason.
            delay = _retry_seconds(str(error))
            if delay is None or attempt == 2:
                return None, str(error)
            time.sleep(delay)
    return None, "Gemini retries exhausted"


def _select_cohorts(
    disease: str,
    studies: list[dict[str, Any]],
    *,
    ledger: ApiCallLedger | None = None,
    loop_round: int = 0,
) -> tuple[list[dict[str, str]] | None, str | None]:
    compact = [
        {
            "gse_id": study["gse_id"],
            "title": study["title"],
            "sample_count": study["sample_count"],
            "cohort_assessment": study["cohort_assessment"],
        }
        for study in studies
    ]
    answer, error = _ask_gemini(
        "You select research cohorts for a bulk RNA-seq biomarker workflow. "
        f"Disease: {disease}. Use only the supplied GEO metadata and its AI assessment. Select exactly three development cohorts and one external cohort only if every chosen study has bulk-like independent samples and comparable tissue/disease definition. "
        "Each study must have cohort_assessment.status=eligible. Do not assign sample groups or override the assessment. Return JSON only: {\"selected\":[{\"gse_id\":\"GSE...\",\"role\":\"development|external\",\"reason\":\"short evidence\"}],\"reason\":\"...\"}. "
        f"STUDIES={json.dumps(compact, ensure_ascii=False)}",
        ledger=ledger,
        phase="orchestrator_cohort_selection",
        loop_round=loop_round,
    )
    if error:
        return None, error
    selected = answer.get("selected", []) if answer else []
    identifiers = {study["gse_id"] for study in studies}
    if (
        len(selected) != 4
        or len({item.get("gse_id") for item in selected}) != 4
        or any(item.get("gse_id") not in identifiers for item in selected)
        or sum(item.get("role") == "development" for item in selected) != 3
        or sum(item.get("role") == "external" for item in selected) != 1
    ):
        return None, (answer or {}).get("reason", "Gemini could not establish four suitable cohorts")
    return selected, None


def _columns(path: Path) -> list[str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        second = handle.readline()
    delimiter = "," if first.count(",") > first.count("\t") else "\t"
    header = next(csv.reader([first], delimiter=delimiter))
    first_row = next(csv.reader([second], delimiter=delimiter), [])
    if len(first_row) == len(header) + 1:
        return header  # GEO matrices sometimes omit the gene-ID header cell.
    if len(first_row) == len(header):
        return header[1:]
    raise ValueError("count matrix header does not match its first data row")


def _require_integer_counts(path: Path) -> None:
    """Reject normalized or fractional matrices before handing them to edgeR."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        delimiter = "," if handle.readline().count(",") > 0 else "\t"
        for row in csv.reader(handle, delimiter=delimiter):
            if len(row) > 1 and any(not value.isdigit() for value in row[1:]):
                raise ValueError("matrix contains non-integer values and is not a raw count matrix")


def _grounded_evidence(sample: dict[str, Any], evidence: list[dict[str, Any]]) -> list[dict[str, str]]:
    available = {
        str(item.get("field", "")): str(item.get("value", ""))
        for item in sample.get("fields", [])
        if isinstance(item, dict)
    }
    grounded: list[dict[str, str]] = []
    for item in evidence:
        if not isinstance(item, dict):
            continue
        field, value = str(item.get("field", "")), str(item.get("value", ""))
        if field in available and value and value.casefold() in available[field].casefold():
            grounded.append({"field": field, "value": value})
    return grounded


def _ai_sample_metadata(samples: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep only per-sample fields that can establish a biological comparison."""
    compact: list[dict[str, Any]] = []
    for sample in samples:
        fields: list[dict[str, str]] = []
        seen: set[tuple[str, str]] = set()
        for item in sample.get("fields", []):
            if not isinstance(item, dict):
                continue
            field = str(item.get("field", ""))
            if not field.startswith(AI_FIELD_PREFIXES):
                continue
            value = " ".join(str(item.get("value", "")).split())[:MAX_AI_FIELD_CHARS]
            if value and (field, value) not in seen:
                fields.append({"field": field, "value": value})
                seen.add((field, value))
        compact.append({"accession": sample["accession"], "fields": fields})
    return compact


def _assessment_from_answer(study: dict[str, Any], answer: dict[str, Any] | None) -> dict[str, Any]:
    metadata = {sample["accession"]: sample for sample in study["sample_metadata"]}
    reviewed: list[dict[str, Any]] = []
    seen: set[str] = set()
    for item in (answer or {}).get("samples", []):
        if not isinstance(item, dict):
            continue
        accession = str(item.get("accession", ""))
        if accession not in metadata or accession in seen:
            continue
        seen.add(accession)
        grounded = _grounded_evidence(metadata[accession], item.get("evidence", []))
        group = str(item.get("group", "unknown"))
        flags = {
            "cell_line": item.get("cell_line") is True,
            "treatment_experiment": item.get("treatment_experiment") is True,
            "not_target_disease": group == "case" and item.get("disease_match") is False,
            "not_clinical_tissue": item.get("clinical_tissue") is False,
            "contradictory_metadata": bool(item.get("contradictions")),
            "ungrounded_evidence": not grounded,
        }
        included = group in {"case", "control"} and not any(flags.values())
        reviewed.append({
            "accession": accession,
            "group": group if included else "unknown",
            "included": included,
            "excluded_by": [name for name, flagged in flags.items() if flagged],
            "evidence": grounded,
            "contradictions": item.get("contradictions", []),
        })
    case_count = sum(item["group"] == "case" for item in reviewed)
    control_count = sum(item["group"] == "control" for item in reviewed)
    status = "eligible" if case_count >= 2 and control_count >= 2 else "insufficient_usable_groups"
    return {
        "status": status,
        "case_count": case_count,
        "control_count": control_count,
        "samples": reviewed,
        "reason": (answer or {}).get("reason", "") if status != "eligible" else "",
    }


def _assess_cohort(
    disease: str,
    study: dict[str, Any],
    *,
    ledger: ApiCallLedger | None = None,
    loop_round: int = 0,
) -> tuple[dict[str, Any], str | None]:
    answer, error = _ask_gemini(
        "Review every supplied GEO sample for a bulk RNA-seq disease-biomarker comparison. "
        f"Target disease: {disease}. STUDY_TITLE={study.get('title', '')}. STUDY_SUMMARY={study.get('summary', '')}. Use only this study context and the supplied GEO fields; do not use papers or outside knowledge. "
        "For each accession return group (case, control, or unknown), disease_match (true, false, or unknown), clinical_tissue (true, false, or unknown), cell_line (true, false, or unknown), treatment_experiment (true, false, or unknown), contradictions (array), and evidence (array of exact {field,value} pairs copied from GEO_FIELDS). "
        "A control may be normal or adjacent non-tumor tissue and need not name the target disease; mark disease_match false only for an explicitly different disease or model. Do not call case-versus-normal labels contradictory. A treatment contrast, cell line, non-target case, non-clinical material, or genuinely conflicting sample label must be marked explicitly. Unknown is allowed; never guess. "
        "Return JSON only: {\"samples\":[{\"accession\":\"GSM...\",\"group\":\"case|control|unknown\",\"disease_match\":true,\"clinical_tissue\":true,\"cell_line\":false,\"treatment_experiment\":false,\"contradictions\":[],\"evidence\":[{\"field\":\"Sample_title\",\"value\":\"...\"}]}],\"reason\":\"...\"}. "
        f"GEO_FIELDS={json.dumps(_ai_sample_metadata(study['sample_metadata']), ensure_ascii=False)}",
        ledger=ledger,
        phase="curator_single_study",
        loop_round=loop_round,
    )
    if error:
        return {"status": "assessment_failed", "reason": error, "samples": []}, error
    metadata = {sample["accession"]: sample for sample in study["sample_metadata"]}
    reviewed: list[dict[str, Any]] = []
    seen: set[str] = set()
    for item in (answer or {}).get("samples", []):
        if not isinstance(item, dict):
            continue
        accession = str(item.get("accession", ""))
        if accession not in metadata or accession in seen:
            continue
        seen.add(accession)
        grounded = _grounded_evidence(metadata[accession], item.get("evidence", []))
        group = str(item.get("group", "unknown"))
        flags = {
            "cell_line": item.get("cell_line") is True,
            "treatment_experiment": item.get("treatment_experiment") is True,
            "not_target_disease": group == "case" and item.get("disease_match") is False,
            "not_clinical_tissue": item.get("clinical_tissue") is False,
            "contradictory_metadata": bool(item.get("contradictions")),
            "ungrounded_evidence": not grounded,
        }
        included = group in {"case", "control"} and not any(flags.values())
        reviewed.append({
            "accession": accession,
            "group": group if included else "unknown",
            "included": included,
            "excluded_by": [name for name, flagged in flags.items() if flagged],
            "evidence": grounded,
            "contradictions": item.get("contradictions", []),
        })
    case_count = sum(item["group"] == "case" for item in reviewed)
    control_count = sum(item["group"] == "control" for item in reviewed)
    status = "eligible" if case_count >= 2 and control_count >= 2 else "insufficient_usable_groups"
    return {
        "status": status,
        "case_count": case_count,
        "control_count": control_count,
        "samples": reviewed,
        "reason": (answer or {}).get("reason", "") if status != "eligible" else "",
    }, None


def _assess_batch(
    disease: str,
    batch: list[dict[str, Any]],
    *,
    ledger: ApiCallLedger | None = None,
    loop_round: int = 0,
) -> dict[str, dict[str, Any]]:
    payload = [
            {
                "gse_id": study["gse_id"],
                "title": study.get("title", ""),
                "summary": study.get("summary", ""),
                "samples": _ai_sample_metadata(study["sample_metadata"]),
            }
            for study in batch
    ]
    answer, error = _ask_gemini(
            "Review supplied GEO studies for a bulk RNA-seq disease-biomarker comparison. "
            f"Target disease: {disease}. Use only the supplied study context and GEO fields; do not use papers or outside knowledge. "
            "Return only confidently usable case/control samples (never unknown); include at most 20 cases and 20 controls per study. For each returned sample include group, disease_match, clinical_tissue, cell_line, treatment_experiment, contradictions, and exact GEO field/value evidence. "
            "Controls may be normal or adjacent tissue and need not name the target disease. Mark disease_match false only for an explicitly different disease or model. "
            "Return JSON only: {\"studies\":[{\"gse_id\":\"GSE...\",\"samples\":[{\"accession\":\"GSM...\",\"group\":\"case|control|unknown\",\"disease_match\":true,\"clinical_tissue\":true,\"cell_line\":false,\"treatment_experiment\":false,\"contradictions\":[],\"evidence\":[{\"field\":\"Sample_title\",\"value\":\"...\"}]}],\"reason\":\"...\"}]}. "
            f"STUDIES={json.dumps(payload, ensure_ascii=False)}",
            ledger=ledger,
            phase="curator_batch",
            loop_round=loop_round,
    )
    reply_items = answer if isinstance(answer, list) else (answer or {}).get("studies", [])
    replies = {str(item.get("gse_id", "")): item for item in reply_items if isinstance(item, dict)}
    assessments: dict[str, dict[str, Any]] = {}
    for study in batch:
        reply = replies.get(study["gse_id"])
        if error or reply is None:
            assessments[study["gse_id"]] = {"status": "assessment_failed", "reason": error or "missing Gemini batch response", "samples": []}
        else:
            assessments[study["gse_id"]] = _assessment_from_answer(study, reply)
    return assessments


def _assess_cohorts(
    disease: str,
    studies: list[dict[str, Any]],
    *,
    ledger: ApiCallLedger | None = None,
    loop_round: int = 0,
) -> dict[str, dict[str, Any]]:
    """Review four studies per Gemini call, with four calls in flight."""
    batches = [studies[start : start + COHORT_BATCH_SIZE] for start in range(0, len(studies), COHORT_BATCH_SIZE)]
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = executor.map(
            lambda batch: _assess_batch(disease, batch, ledger=ledger, loop_round=loop_round),
            batches,
        )
        assessments = {gse_id: assessment for result in results for gse_id, assessment in result.items()}
    if ledger is not None:
        # The proposal's initial budget counts one Curator call per batch. A failed
        # batch is preserved for inspection instead of silently issuing unbounded retries.
        return assessments
    for study in studies:
        assessment = assessments[study["gse_id"]]
        if assessment["status"] == "assessment_failed" or (assessment["status"] == "insufficient_usable_groups" and max(assessment["case_count"], assessment["control_count"]) >= 2):
            assessments[study["gse_id"]], _ = _assess_cohort(disease, study)
    return assessments


def _groups_from_assessment(assessment: dict[str, Any], columns: list[str], crosswalk: list[dict[str, str]]) -> tuple[list[dict[str, Any]] | None, str | None]:
    labels = {item["accession"]: item for item in assessment.get("samples", []) if item.get("included")}
    matched = [item for item in crosswalk if item.get("accession") in labels]
    accession_counts: dict[str, int] = {}
    for item in matched:
        accession_counts[item["accession"]] = accession_counts.get(item["accession"], 0) + 1
    groups = [
        {"sample_id": item["raw_column"], "group": labels[item["accession"]]["group"], "evidence": labels[item["accession"]]["evidence"]}
        for item in matched
        if item["raw_column"] in columns and accession_counts[item["accession"]] == 1
    ]
    if sum(item["group"] == "case" for item in groups) < 2 or sum(item["group"] == "control" for item in groups) < 2:
        return None, "fewer than two uniquely mapped usable case and control samples"
    return groups, None


def _prepare_matrix(raw_path: Path, groups: list[dict[str, str]], destination: Path) -> None:
    wanted = [item["sample_id"] for item in groups]
    opener = gzip.open if raw_path.suffix.lower() == ".gz" else open
    with opener(raw_path, "rt", encoding="utf-8-sig", newline="") as source:
        first = source.readline()
        delimiter = "," if first.count(",") > first.count("\t") else "\t"
        reader = csv.reader(source, delimiter=delimiter)
        header = next(csv.reader([first], delimiter=delimiter))
        first_row = next(reader, [])
        if len(first_row) == len(header) + 1:
            indexes = [header.index(sample) + 1 for sample in wanted]
        elif len(first_row) == len(header):
            indexes = [header.index(sample) for sample in wanted]
        else:
            raise ValueError("count matrix header does not match its first data row")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target)
            writer.writerow(["gene_id", *wanted])
            for row in chain([first_row], reader):
                if len(row) > max(indexes):
                    writer.writerow([row[0], *[row[index] for index in indexes]])


def _clinical_candidates(disease: str, organism: str, limit: int) -> tuple[list[str], list[dict[str, Any]]]:
    """Prefer explicit normal/paired studies, then fill the review pool with clinical studies."""
    queries = [
        f"({disease}) AND (patient OR patients OR tissue) AND (normal OR adjacent OR paired)",
        f"({disease}) AND (patient OR patients OR tissue)",
    ]
    studies: list[dict[str, Any]] = []
    seen: set[str] = set()
    for query in queries:
        for study in find_bulk_rnaseq(query, organism, limit):
            if study["gse_id"] not in seen:
                studies.append(study)
                seen.add(study["gse_id"])
            if len(studies) == limit:
                return queries, studies
    return queries, studies


def _prepare_cohorts(
    output_dir: Path,
    selected: list[dict[str, str]],
    assessed_studies: list[dict[str, Any]],
    eligible: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, str]], dict[str, dict[str, Any]]]:
    """Download and crosswalk selected studies, filling failed slots safely."""
    by_id = {study["gse_id"]: study for study in assessed_studies}
    prepared: list[dict[str, Any]] = []
    rejected: list[dict[str, str]] = []
    pending = list(selected)
    reserved = {choice["gse_id"] for choice in selected}
    while pending:
        choice = pending.pop(0)
        study = by_id[choice["gse_id"]]
        cohort_dir = output_dir / "cohorts" / study["gse_id"]
        raw_dir = cohort_dir / "raw"
        try:
            filename = study["raw_count_files"][0]["filename"]
            download = fetch_raw_count_matrix(study["gse_id"], filename, raw_dir)
            raw_path = Path(download["path"])
            columns = _columns(raw_path)
            _require_integer_counts(raw_path)
            crosswalk = match_count_columns_to_samples(columns, study["sample_metadata"])
            _write_json(cohort_dir / "input" / "crosswalk.json", crosswalk)
            groups, group_error = _groups_from_assessment(study["cohort_assessment"], columns, crosswalk)
            if not groups:
                raise ValueError(group_error)
            samples_path = cohort_dir / "input" / "samples.csv"
            samples_path.parent.mkdir(parents=True, exist_ok=True)
            with samples_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_id", "group"])
                writer.writeheader()
                writer.writerows({key: row[key] for key in ("sample_id", "group")} for row in groups)
            counts_path = cohort_dir / "input" / "counts.csv"
            _prepare_matrix(raw_path, groups, counts_path)
            prepared.append({
                "name": study["gse_id"],
                "role": choice["role"],
                "counts": str(counts_path),
                "samples": str(samples_path),
                "selection_reason": choice.get("reason", ""),
            })
        except Exception as error:
            rejected.append({"gse_id": study["gse_id"], "reason": str(error)})
            if cohort_dir.exists():
                shutil.rmtree(cohort_dir)
            backup = next((item for item in eligible if item["gse_id"] not in reserved), None)
            if backup:
                reserved.add(backup["gse_id"])
                pending.append({
                    "gse_id": backup["gse_id"],
                    "role": choice["role"],
                    "reason": f"automatic fallback after {study['gse_id']} failed",
                })
    return prepared, rejected, by_id


def _analysis_config(
    analysis_dir: Path,
    prepared: list[dict[str, Any]],
    *,
    max_loop_rounds: int | None = None,
) -> dict[str, Any]:
    external = next(item for item in prepared if item["role"] == "external")
    config = {
        "fdr": 0.05,
        "min_validation_auc": 0.8,
        "min_validation_sensitivity": 0.7,
        "min_validation_specificity": 0.7,
        "output_dir": str(analysis_dir / "cohorts"),
        "validation_output_dir": str(analysis_dir / "validation"),
        "development": [{key: item[key] for key in ("name", "counts", "samples")} for item in prepared if item["role"] == "development"],
        "external": {key: external[key] for key in ("name", "counts", "samples")},
    }
    if max_loop_rounds is not None:
        config["max_loop_rounds"] = max_loop_rounds
    return config


def _persist_cohort_selection(
    output_dir: Path,
    provenance: Path,
    prepared: list[dict[str, Any]],
    by_id: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    cohorts = [_cohort_summary(by_id[item["name"]], item) for item in prepared]
    for item in prepared:
        _write_json(output_dir / "cohorts" / item["name"] / "metadata.json", by_id[item["name"]])
    _write_json(provenance / "cohorts.json", cohorts)
    return cohorts


def _deterministic_development_critic(analysis: dict[str, Any], cohort_id: str) -> dict[str, Any]:
    """Produce a local safety check when an ADK critic is unavailable."""
    rows = [row for row in analysis.get("development_validation", []) if str(row.get("discovery")) == cohort_id]
    passed = [row for row in analysis.get("candidate_ranking", []) if str(row.get("discovery")) == cohort_id and row.get("passed")]
    reason_codes: list[str] = []
    if not rows:
        reason_codes.append("insufficient_internal_validation")
    if not passed:
        reason_codes.append("no_locked_candidates")
    return {
        "engine": "deterministic-fallback",
        "scope": "development",
        "cohort_id": cohort_id,
        "verdict": "interpretable" if not reason_codes else "caution",
        "summary": f"{cohort_id}: 내부 검증 행 {len(rows)}개, 기준 통과 후보 {len(passed)}개",
        "reason_codes": reason_codes,
    }


def _run_development_critics(
    output_dir: Path,
    analysis: dict[str, Any],
    config: dict[str, Any],
    validation_dir: Path,
    ledger: ApiCallLedger,
    loop_round: int,
) -> list[dict[str, Any]]:
    reports: list[dict[str, Any]] = []
    for item in config.get("development", []):
        cohort_id = str(item["name"])
        try:
            ledger.reserve("adk", "development_critic", loop_round)
            report = run_development_critic(output_dir, cohort_id)
        except Exception:
            report = _deterministic_development_critic(analysis, cohort_id)
        reports.append(report)
        _write_json(validation_dir / f"critic_{cohort_id}.json", report)
    return reports


def _retry_reason_codes(analysis: dict[str, Any], critic_reports: list[dict[str, Any]]) -> list[str]:
    reason_codes = {
        str(code)
        for report in critic_reports
        for code in report.get("reason_codes", [])
    }
    if not any(row.get("passed") for row in analysis.get("candidate_ranking", [])):
        reason_codes.add("no_locked_candidates")
    return sorted(reason_codes)


def _record_analysis(
    output_dir: Path,
    manifest: dict[str, Any],
    analysis: dict[str, Any],
    *,
    finalize: bool = True,
    ledger: ApiCallLedger | None = None,
    loop_round: int = 0,
) -> dict[str, Any]:
    analysis_dir = output_dir / "analysis"
    validation_dir = analysis_dir / "validation"
    ranking = analysis["candidate_ranking"]
    biomarkers = [row for row in ranking if row["passed"]]
    _write_json(analysis_dir / "biomarkers.json", biomarkers)
    _write_csv(analysis_dir / "biomarkers.csv", biomarkers)
    _write_csv(analysis_dir / "candidate_scores.csv", ranking)
    internal_validation = analysis.get("development_validation", [])
    _write_json(validation_dir / "internal_validation.json", internal_validation)
    _write_validation_csv(validation_dir / "internal_validation.csv", internal_validation)
    external_validation = analysis.get("external_validation", [])
    _write_json(validation_dir / "external_validation.json", external_validation)
    _write_validation_csv(validation_dir / "external_validation.csv", external_validation)
    manifest["analysis"] = {"candidate_score_count": len(ranking), "biomarker_count": len(biomarkers), "external_validation_count": len(external_validation), "top_biomarkers": biomarkers[:5], "fdr": analysis["fdr"], "min_validation_auc": analysis["min_validation_auc"], "min_validation_sensitivity": analysis["min_validation_sensitivity"], "min_validation_specificity": analysis["min_validation_specificity"], "validation_engine": analysis.get("validation_engine", "unknown"), "validation_metrics": analysis.get("validation_metrics", [])}
    result: dict[str, Any] = {"critic": None, "summary": None}
    if not finalize:
        return result
    try:
        if ledger is not None:
            ledger.reserve("adk", "external_critic", loop_round)
        critic = run_critic(output_dir)
        _write_json(validation_dir / "critic_report.json", critic)
        result["critic"] = critic
        manifest["critic"] = {key: critic.get(key) for key in ("engine", "verdict", "summary", "reason_codes")}
    except Exception as error:
        manifest["critic"] = {"status": "unavailable", "reason": str(error)}
    try:
        if ledger is not None:
            ledger.reserve("adk", "summarizer", loop_round)
        summary = run_summarizer(output_dir)
        report = run_summarizer_agent(summary)
        report_dir = output_dir / "report"
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "summary_report.md").write_text(report, encoding="utf-8")
        result["summary"] = summary
        manifest["summarizer"] = {"engine": "google-adk", "primary_development_cohort": summary["primary_development_cohort"], "top_gene_count": len(summary["top_five"])}
    except Exception as error:
        manifest["summarizer"] = {"status": "unavailable", "reason": str(error)}
    return result


def resume_analysis(output_dir: Path = Path("output")) -> dict[str, Any]:
    """Finish a run from retained edgeR artifacts without repeating GEO collection."""
    output_dir = output_dir.resolve()
    provenance = output_dir / "provenance"
    manifest = json.loads((provenance / "run_manifest.json").read_text(encoding="utf-8"))
    try:
        _record_analysis(output_dir, manifest, run_pipeline(provenance / "analysis_config.json"))
        manifest["stage"] = "complete"
        manifest.pop("reason", None)
    except Exception as error:
        manifest.update({"stage": "analysis_failed", "reason": str(error)})
    _write_json(provenance / "run_manifest.json", manifest)
    return manifest


def run_disease(
    disease: str,
    output_dir: Path = Path("output"),
    organism: str = "Homo sapiens",
    limit: int = 10,
    max_loop_rounds: int = 0,
) -> dict[str, Any]:
    """Collect, validate, analyse, and optionally retry a disease workflow.

    The library default keeps the original one-pass behaviour for callers that
    already manage retries. The CLI enables the proposal's two bounded retry
    rounds explicitly.
    """
    if max_loop_rounds < 0:
        raise ValueError("max_loop_rounds must be zero or greater")
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = output_dir / "provenance"
    manifest: dict[str, Any] = {
        "disease": disease,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "searching",
        "studies": [],
        "max_loop_rounds": max_loop_rounds,
        "loops": [],
        "agent": {
            "provider": "google-gemini",
            "model": os.getenv("GEMINI_MODEL", "gemini-3.5-flash-lite"),
            "curator_batch_size": COHORT_BATCH_SIZE,
        },
    }
    _write_json(provenance / "run_manifest.json", manifest)
    candidate_limit = min(25, max(20, limit))
    try:
        clinical_queries, studies = _clinical_candidates(disease, organism, candidate_limit)
    except Exception as error:
        manifest.update({"stage": "search_failed", "reason": str(error), "candidate_limit": candidate_limit})
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest
    manifest["clinical_queries"] = clinical_queries
    manifest["candidate_limit"] = candidate_limit
    manifest["studies"] = [_candidate_preview(study) for study in studies]
    _write_json(provenance / "search_results.json", manifest["studies"])
    ledger = ApiCallLedger(len(studies), max_loop_rounds)
    manifest["api_budget"] = ledger.snapshot()
    try:
        studies_with_metadata = [{**study, **fetch_geo_sample_metadata(study["gse_id"])} for study in studies]
        assessments = _assess_cohorts(disease, studies_with_metadata, ledger=ledger)
    except TypeError as error:
        # Keep compatibility with small test doubles and older integrations that
        # still expose the two-argument Curator function.
        if "ledger" not in str(error):
            raise
        studies_with_metadata = [{**study, **fetch_geo_sample_metadata(study["gse_id"])} for study in studies]
        assessments = _assess_cohorts(disease, studies_with_metadata)
    except Exception as error:
        manifest.update({"stage": "cohort_assessment_failed", "reason": str(error), "api_budget": ledger.snapshot()})
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest
    assessed_studies = [{**study, "cohort_assessment": assessments[study["gse_id"]]} for study in studies_with_metadata]
    _write_json(provenance / "curator_assessments.json", assessed_studies)
    manifest["stage"] = "selecting_cohorts"
    manifest["api_budget"] = ledger.snapshot()
    _write_json(provenance / "run_manifest.json", manifest)

    eligible = [study for study in assessed_studies if study["cohort_assessment"]["status"] == "eligible"]
    if len(eligible) < 4:
        failed = [study for study in assessed_studies if study["cohort_assessment"]["status"] == "assessment_failed"]
        if failed and len(failed) == len(assessed_studies):
            manifest.update({"stage": "cohort_assessment_failed", "reason": failed[0]["cohort_assessment"].get("reason", "Curator assessment failed")})
        else:
            manifest.update({"stage": "needs_cohort_review", "reason": f"only {len(eligible)} studies have at least two usable case and control samples after AI metadata exclusions"})
        manifest["api_budget"] = ledger.snapshot()
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest

    try:
        selected, error = _select_cohorts(disease, eligible, ledger=ledger)
    except TypeError as type_error:
        if "ledger" not in str(type_error):
            raise
        selected, error = _select_cohorts(disease, eligible)
    if not selected:
        manifest.update({"stage": "needs_cohort_review", "reason": error})
        manifest["api_budget"] = ledger.snapshot()
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest

    analysis_dir = output_dir / "analysis"
    current_selected = selected
    current_assessed = assessed_studies
    current_eligible = eligible
    finalised = False

    def finish_analysis(analysis: dict[str, Any], loop_round: int) -> bool:
        """Evaluate the held-out cohort once, then run the final critic and summary."""
        try:
            if not analysis.get("external_evaluated", False):
                analysis = run_pipeline(provenance / "analysis_config.json", include_external=True)
        except Exception as error:
            manifest.update({"stage": "analysis_failed", "reason": str(error)})
            return False
        _record_analysis(output_dir, manifest, analysis, ledger=ledger, loop_round=loop_round)
        return True

    for round_index in range(max_loop_rounds + 1):
        prepared, rejected, by_id = _prepare_cohorts(output_dir, current_selected, current_assessed, current_eligible)
        if len(prepared) != 4:
            manifest.update({"stage": "needs_sample_review", "reason": "fewer than four cohorts could be prepared", "rejected": rejected, "prepared": prepared})
            break
        cohorts = _persist_cohort_selection(output_dir, provenance, prepared, by_id)
        config = _analysis_config(analysis_dir, prepared, max_loop_rounds=max_loop_rounds)
        _write_json(provenance / "analysis_config.json", config)
        manifest.update({"stage": "analysing", "cohorts": cohorts, "rejected_count": len(rejected), "active_loop_round": round_index})
        manifest["api_budget"] = ledger.snapshot()
        _write_json(provenance / "run_manifest.json", manifest)
        try:
            analysis = run_pipeline(
                provenance / "analysis_config.json",
                include_external=max_loop_rounds == 0,
            )
        except Exception as error:
            manifest.update({"stage": "analysis_failed", "reason": str(error)})
            break

        if max_loop_rounds == 0:
            _record_analysis(output_dir, manifest, analysis, ledger=ledger, loop_round=round_index)
            manifest.update({"stage": "complete", "stop_reason": "single_pass_complete"})
            finalised = True
            break

        _record_analysis(output_dir, manifest, analysis, finalize=False, ledger=ledger, loop_round=round_index)
        critic_reports = _run_development_critics(
            output_dir,
            analysis,
            config,
            analysis_dir / "validation",
            ledger,
            round_index,
        )
        reason_codes = _retry_reason_codes(analysis, critic_reports)
        loop_record = {
            "round": round_index,
            "reason_codes": reason_codes,
            "development_critics": critic_reports,
            "action": "retry" if set(reason_codes) & RETRYABLE_REASON_CODES and round_index < max_loop_rounds else "stop",
        }
        manifest.setdefault("loops", []).append(loop_record)
        if not set(reason_codes) & RETRYABLE_REASON_CODES or round_index >= max_loop_rounds:
            if finish_analysis(analysis, round_index):
                manifest.update({
                    "stage": "complete",
                    "stop_reason": "max_loop_rounds_reached" if round_index >= max_loop_rounds and set(reason_codes) & RETRYABLE_REASON_CODES else "critic_passed",
                })
                finalised = True
            break

        current_ids = {choice["gse_id"] for choice in current_selected}
        unused = [study for study in current_assessed if study["gse_id"] not in current_ids]
        if not unused:
            loop_record["action"] = "stop_no_replacement_cohort"
            if finish_analysis(analysis, round_index):
                manifest.update({"stage": "complete", "stop_reason": "no_replacement_cohort"})
                finalised = True
            break
        replacement_batch = unused[:COHORT_BATCH_SIZE]
        try:
            ledger.open_loop(len(replacement_batch))
            try:
                replacement_assessments = _assess_cohorts(
                    disease,
                    replacement_batch,
                    ledger=ledger,
                    loop_round=round_index + 1,
                )
            except TypeError as type_error:
                if "ledger" not in str(type_error):
                    raise
                replacement_assessments = _assess_cohorts(disease, replacement_batch)
        except Exception as error:
            loop_record.update({"action": "stop_budget_or_curator_error", "error": str(error)})
            if finish_analysis(analysis, round_index):
                manifest.update({"stage": "complete", "stop_reason": "curator_or_budget_error"})
                finalised = True
            break
        updated_by_id = {study["gse_id"]: study for study in current_assessed}
        for study in replacement_batch:
            updated_by_id[study["gse_id"]] = {**study, "cohort_assessment": replacement_assessments[study["gse_id"]]}
        current_assessed = list(updated_by_id.values())
        current_eligible = [study for study in current_assessed if study["cohort_assessment"]["status"] == "eligible"]
        assessed_replacements = [updated_by_id[study["gse_id"]] for study in replacement_batch]
        _write_json(provenance / f"loop_{round_index + 1}_curator_assessments.json", assessed_replacements)
        replacement_ids = {study["gse_id"] for study in replacement_batch}
        candidate_pool = [
            study for study in current_assessed
            if study["gse_id"] in current_ids or study["gse_id"] in replacement_ids
        ]
        try:
            next_selected, selection_error = _select_cohorts(
                disease,
                [study for study in candidate_pool if study in current_eligible],
                ledger=ledger,
                loop_round=round_index + 1,
            )
        except TypeError as type_error:
            if "ledger" not in str(type_error):
                raise
            next_selected, selection_error = _select_cohorts(disease, [study for study in candidate_pool if study in current_eligible])
        if not next_selected or {choice["gse_id"] for choice in next_selected} == current_ids:
            backup = next((study for study in replacement_batch if study["gse_id"] in {item["gse_id"] for item in current_eligible}), None)
            if backup:
                next_selected = [dict(choice) for choice in current_selected]
                development_index = next(index for index, choice in enumerate(next_selected) if choice["role"] == "development")
                next_selected[development_index] = {"gse_id": backup["gse_id"], "role": "development", "reason": "replacement after critic review"}
            else:
                loop_record.update({"action": "stop_no_eligible_replacement", "error": selection_error or "no replacement cohort"})
                if finish_analysis(analysis, round_index):
                    manifest.update({"stage": "complete", "stop_reason": "no_eligible_replacement"})
                    finalised = True
                break
        current_selected = next_selected
        manifest["api_budget"] = ledger.snapshot()
        _write_json(provenance / "run_manifest.json", manifest)

    if not finalised and manifest.get("stage") not in {"complete", "analysis_failed", "needs_sample_review"}:
        manifest["stage"] = "complete"
    manifest["api_budget"] = ledger.snapshot()
    _write_json(provenance / "run_manifest.json", manifest)
    return manifest
