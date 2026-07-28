"""One-command GEO collection and safe handoff to the edgeR workflow."""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path
from typing import Any

from dotenv import load_dotenv

from agent import run_critic, run_summarizer_agent
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


def _ask_gemini(prompt: str) -> tuple[dict[str, Any] | None, str | None]:
    if not os.getenv("GEMINI_API_KEY"):
        return None, "GEMINI_API_KEY is not configured"
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


def _select_cohorts(disease: str, studies: list[dict[str, Any]]) -> tuple[list[dict[str, str]] | None, str | None]:
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
        f"STUDIES={json.dumps(compact, ensure_ascii=False)}"
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


def _assess_cohort(disease: str, study: dict[str, Any]) -> tuple[dict[str, Any], str | None]:
    answer, error = _ask_gemini(
        "Review every supplied GEO sample for a bulk RNA-seq disease-biomarker comparison. "
        f"Target disease: {disease}. STUDY_TITLE={study.get('title', '')}. STUDY_SUMMARY={study.get('summary', '')}. Use only this study context and the supplied GEO fields; do not use papers or outside knowledge. "
        "For each accession return group (case, control, or unknown), disease_match (true, false, or unknown), clinical_tissue (true, false, or unknown), cell_line (true, false, or unknown), treatment_experiment (true, false, or unknown), contradictions (array), and evidence (array of exact {field,value} pairs copied from GEO_FIELDS). "
        "A control may be normal or adjacent non-tumor tissue and need not name the target disease; mark disease_match false only for an explicitly different disease or model. Do not call case-versus-normal labels contradictory. A treatment contrast, cell line, non-target case, non-clinical material, or genuinely conflicting sample label must be marked explicitly. Unknown is allowed; never guess. "
        "Return JSON only: {\"samples\":[{\"accession\":\"GSM...\",\"group\":\"case|control|unknown\",\"disease_match\":true,\"clinical_tissue\":true,\"cell_line\":false,\"treatment_experiment\":false,\"contradictions\":[],\"evidence\":[{\"field\":\"Sample_title\",\"value\":\"...\"}]}],\"reason\":\"...\"}. "
        f"GEO_FIELDS={json.dumps(_ai_sample_metadata(study['sample_metadata']), ensure_ascii=False)}"
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


def _assess_batch(disease: str, batch: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
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
            f"STUDIES={json.dumps(payload, ensure_ascii=False)}"
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


def _assess_cohorts(disease: str, studies: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Review four studies per Gemini call, with four calls in flight."""
    batches = [studies[start : start + COHORT_BATCH_SIZE] for start in range(0, len(studies), COHORT_BATCH_SIZE)]
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = executor.map(lambda batch: _assess_batch(disease, batch), batches)
        assessments = {gse_id: assessment for result in results for gse_id, assessment in result.items()}
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


def _record_analysis(output_dir: Path, manifest: dict[str, Any], analysis: dict[str, Any]) -> None:
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
    manifest["analysis"] = {"candidate_score_count": len(ranking), "biomarker_count": len(biomarkers), "external_validation_count": len(external_validation), "top_biomarkers": biomarkers[:5], "fdr": analysis["fdr"], "min_validation_auc": analysis["min_validation_auc"], "min_validation_sensitivity": analysis["min_validation_sensitivity"], "min_validation_specificity": analysis["min_validation_specificity"]}
    try:
        critic = run_critic(output_dir)
        _write_json(validation_dir / "critic_report.json", critic)
        manifest["critic"] = {key: critic.get(key) for key in ("engine", "verdict", "summary")}
    except Exception as error:
        manifest["critic"] = {"status": "unavailable", "reason": str(error)}
    try:
        summary = run_summarizer(output_dir)
        report = run_summarizer_agent(summary)
        report_dir = output_dir / "report"
        report_dir.mkdir(parents=True, exist_ok=True)
        (report_dir / "summary_report.md").write_text(report, encoding="utf-8")
        manifest["summarizer"] = {"engine": "google-adk", "primary_development_cohort": summary["primary_development_cohort"], "top_gene_count": len(summary["top_five"])}
    except Exception as error:
        manifest["summarizer"] = {"status": "unavailable", "reason": str(error)}


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


def run_disease(disease: str, output_dir: Path = Path("output"), organism: str = "Homo sapiens", limit: int = 10) -> dict[str, Any]:
    """Collect, preserve, select, and analyse a disease in one command when evidence permits."""
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = output_dir / "provenance"
    manifest: dict[str, Any] = {
        "disease": disease,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stage": "searching",
        "studies": [],
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
    studies_with_metadata = [{**study, **fetch_geo_sample_metadata(study["gse_id"])} for study in studies]
    assessments = _assess_cohorts(disease, studies_with_metadata)
    assessed_studies = [{**study, "cohort_assessment": assessments[study["gse_id"]]} for study in studies_with_metadata]
    manifest["stage"] = "selecting_cohorts"
    _write_json(provenance / "run_manifest.json", manifest)

    eligible = [study for study in assessed_studies if study["cohort_assessment"]["status"] == "eligible"]
    if len(eligible) < 4:
        failed = [study for study in assessed_studies if study["cohort_assessment"]["status"] == "assessment_failed"]
        if len(failed) == len(assessed_studies):
            manifest.update({"stage": "cohort_assessment_failed", "reason": failed[0]["cohort_assessment"]["reason"]})
            _write_json(provenance / "run_manifest.json", manifest)
            return manifest
        manifest.update({"stage": "needs_cohort_review", "reason": f"only {len(eligible)} studies have at least two usable case and control samples after AI metadata exclusions"})
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest
    selected, error = _select_cohorts(disease, eligible)
    if not selected:
        manifest.update({"stage": "needs_cohort_review", "reason": error})
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest

    by_id = {study["gse_id"]: study for study in assessed_studies}
    analysis_dir = output_dir / "analysis"
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
            prepared.append({"name": study["gse_id"], "role": choice["role"], "counts": str(counts_path), "samples": str(samples_path), "selection_reason": choice.get("reason", "")})
        except Exception as error:
            rejected.append({"gse_id": study["gse_id"], "reason": str(error)})
            if cohort_dir.exists():
                shutil.rmtree(cohort_dir)
            backup = next((item for item in eligible if item["gse_id"] not in reserved), None)
            if backup:
                reserved.add(backup["gse_id"])
                pending.append({"gse_id": backup["gse_id"], "role": choice["role"], "reason": f"automatic fallback after {study['gse_id']} failed"})
            continue

    if len(prepared) != 4:
        manifest.update({"stage": "needs_sample_review", "reason": "fewer than four cohorts could be prepared", "rejected": rejected, "prepared": prepared})
        _write_json(provenance / "run_manifest.json", manifest)
        return manifest

    cohorts = [_cohort_summary(by_id[item["name"]], item) for item in prepared]
    for item in prepared:
        _write_json(output_dir / "cohorts" / item["name"] / "metadata.json", by_id[item["name"]])
    _write_json(provenance / "cohorts.json", cohorts)
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
    _write_json(provenance / "analysis_config.json", config)
    manifest.update({"stage": "analysing", "cohorts": cohorts, "rejected_count": len(rejected)})
    _write_json(provenance / "run_manifest.json", manifest)
    try:
        analysis = run_pipeline(provenance / "analysis_config.json")
        _record_analysis(output_dir, manifest, analysis)
        manifest["stage"] = "complete"
    except Exception as error:
        manifest.update({"stage": "analysis_failed", "reason": str(error)})
    _write_json(provenance / "run_manifest.json", manifest)
    return manifest
