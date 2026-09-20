"""Deterministic edgeR execution and locked-candidate validation."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from statistics import mean
from typing import Any
from urllib.request import Request, urlopen

try:
    from sklearn.metrics import confusion_matrix, roc_auc_score
except ImportError:  # pragma: no cover - pyproject installs scikit-learn in normal runs.
    confusion_matrix = None
    roc_auc_score = None


ENSEMBL_BATCH_SIZE = 200
ENSEMBL_CACHE_NAME = "ensembl_gene_id_cache.json"
ENSEMBL_TIMEOUT_SECONDS = 5


@dataclass(frozen=True)
class Cohort:
    name: str
    counts: Path
    samples: Path


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _rscript() -> str:
    return shutil.which("Rscript") or r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"


def run_edger(cohort: Cohort, output_dir: Path) -> Path:
    script = Path(__file__).with_name("edger.R")
    subprocess.run([_rscript(), str(script), str(cohort.counts), str(cohort.samples), str(output_dir)], check=True)
    return output_dir


def _edger_complete(output_dir: Path) -> bool:
    return (output_dir / "edger_results.csv").is_file() and (output_dir / "logcpm.csv").is_file()


def _cohort_output(output_dir: Path, name: str) -> Path:
    return output_dir / name / "edgeR"


@lru_cache(maxsize=16)
def _logcpm(path: Path) -> dict[str, dict[str, float]]:
    rows = _rows(path)
    return {row["gene_id"]: {sample: float(value) for sample, value in row.items() if sample != "gene_id"} for row in rows}


def _groups(path: Path) -> dict[str, str]:
    return {row["sample_id"]: row["group"] for row in _rows(path)}


def _auc(scores: list[float], labels: list[bool]) -> float:
    positives, negatives = sum(labels), len(labels) - sum(labels)
    if not positives or not negatives:
        raise ValueError("validation cohort needs case and control samples")
    if roc_auc_score is not None:
        return float(roc_auc_score(labels, scores))
    ordered = sorted(enumerate(scores), key=lambda pair: pair[1])
    ranks = [0.0] * len(scores)
    index = 0
    while index < len(ordered):
        end = index
        while end + 1 < len(ordered) and ordered[end + 1][1] == ordered[index][1]:
            end += 1
        rank = (index + end + 2) / 2
        for position in range(index, end + 1):
            ranks[ordered[position][0]] = rank
        index = end + 1
    return (sum(rank for rank, label in zip(ranks, labels) if label) - positives * (positives + 1) / 2) / (positives * negatives)


def _threshold(values: dict[str, float], groups: dict[str, str], direction: float) -> float:
    case = [value for sample, value in values.items() if groups.get(sample) == "case"]
    control = [value for sample, value in values.items() if groups.get(sample) == "control"]
    if not case or not control:
        raise ValueError("discovery cohort needs case and control samples")
    return direction * (mean(case) + mean(control)) / 2


def _ensembl_cache_path() -> Path:
    return Path.cwd() / ".rnahero" / ENSEMBL_CACHE_NAME


def _ensembl_cache() -> dict[str, dict[str, str]]:
    path = _ensembl_cache_path()
    try:
        cache = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        cache = {}
    return {
        "symbol_to_ensembl": cache.get("symbol_to_ensembl", {}),
        "ensembl_to_symbol": cache.get("ensembl_to_symbol", {}),
    }


def _save_ensembl_cache(cache: dict[str, dict[str, str]]) -> None:
    path = _ensembl_cache_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def _ensembl_for_symbols(symbols: list[str]) -> dict[str, str]:
    """Map human gene symbols to Ensembl gene IDs in bounded API batches."""
    cache = _ensembl_cache()
    mappings = cache["symbol_to_ensembl"]
    missing = [symbol for symbol in symbols if symbol not in mappings]
    for start in range(0, len(missing), ENSEMBL_BATCH_SIZE):
        payload = json.dumps({"symbols": missing[start : start + ENSEMBL_BATCH_SIZE]}).encode("utf-8")
        request = Request(
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "RNAHero/0.1"},
            method="POST",
        )
        try:
            with urlopen(request, timeout=ENSEMBL_TIMEOUT_SECONDS) as response:
                result = json.loads(response.read())
        except Exception:
            # ponytail: mapping is optional; unmapped candidates are recorded as skipped downstream.
            break
        for symbol, item in result.items():
            if isinstance(item, dict) and item.get("object_type") == "Gene" and str(item.get("id", "")).startswith("ENSG"):
                mappings[symbol] = item["id"]
    if missing:
        _save_ensembl_cache(cache)
    return {symbol: mappings[symbol] for symbol in symbols if symbol in mappings}


def _symbols_for_ensembl(ids: list[str]) -> dict[str, str]:
    """Map Ensembl gene IDs to their human display symbols in bounded API batches."""
    cache = _ensembl_cache()
    mappings = cache["ensembl_to_symbol"]
    missing = [identifier for identifier in ids if identifier not in mappings]
    for start in range(0, len(missing), ENSEMBL_BATCH_SIZE):
        payload = json.dumps({"ids": missing[start : start + ENSEMBL_BATCH_SIZE]}).encode("utf-8")
        request = Request(
            "https://rest.ensembl.org/lookup/id",
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "RNAHero/0.1"},
            method="POST",
        )
        try:
            with urlopen(request, timeout=ENSEMBL_TIMEOUT_SECONDS) as response:
                result = json.loads(response.read())
        except Exception:
            # ponytail: mapping is optional; unmapped candidates are recorded as skipped downstream.
            break
        for identifier, item in result.items():
            if isinstance(item, dict) and item.get("object_type") == "Gene" and item.get("display_name"):
                mappings[identifier] = str(item["display_name"])
    if missing:
        _save_ensembl_cache(cache)
    return {identifier: mappings[identifier] for identifier in ids if identifier in mappings}


def _validation_gene(gene: str, values: dict[str, dict[str, float]], symbol_to_ensembl: dict[str, str], ensembl_to_symbol: dict[str, str], ensembl_values: dict[str, str] | None = None) -> str | None:
    if gene in values:
        return gene
    ensembl = gene.split(".", 1)[0] if gene.startswith("ENSG") else symbol_to_ensembl.get(gene)
    ensembl_values = ensembl_values or {identifier.split(".", 1)[0]: identifier for identifier in values}
    if ensembl and ensembl in ensembl_values:
        return ensembl_values[ensembl]
    symbol = ensembl_to_symbol.get(ensembl, "") if ensembl else ""
    return symbol if symbol in values else None


def score_gene(gene: str, direction: float, discovery: Cohort, validation: Cohort, outputs: Path, validation_gene: str | None = None, discovery_matrix: dict[str, dict[str, float]] | None = None, validation_matrix: dict[str, dict[str, float]] | None = None, discovery_groups: dict[str, str] | None = None, validation_groups: dict[str, str] | None = None) -> dict[str, float | str]:
    discovery_values = (discovery_matrix or _logcpm(_cohort_output(outputs, discovery.name) / "logcpm.csv"))[gene]
    validation_gene = validation_gene or gene
    validation_values = (validation_matrix or _logcpm(_cohort_output(outputs, validation.name) / "logcpm.csv")).get(validation_gene)
    if validation_values is None:
        raise ValueError(f"{validation_gene} is missing from {validation.name}")
    discovery_groups = discovery_groups or _groups(discovery.samples)
    groups = validation_groups or _groups(validation.samples)
    cutoff = _threshold(discovery_values, discovery_groups, direction)
    samples = [sample for sample in validation_values if sample in groups]
    scores = [direction * validation_values[sample] for sample in samples]
    labels = [groups[sample] == "case" for sample in samples]
    predicted = [score >= cutoff for score in scores]
    if confusion_matrix is not None:
        tn, fp, fn, tp = confusion_matrix(labels, predicted, labels=[False, True]).ravel()
    else:
        tp = sum(prediction and label for prediction, label in zip(predicted, labels))
        tn = sum(not prediction and not label for prediction, label in zip(predicted, labels))
        fp = sum(prediction and not label for prediction, label in zip(predicted, labels))
        fn = sum(not prediction and label for prediction, label in zip(predicted, labels))
    return {"gene_id": gene, "validation_gene_id": validation_gene, "discovery": discovery.name, "cohort": validation.name, "direction": direction, "auc": _auc(scores, labels), "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp), "cutoff": cutoff}


def _candidates(results: Path, fdr: float) -> list[tuple[str, float]]:
    if not 0 < fdr <= 0.1:
        raise ValueError("fdr must be in (0, 0.1]; use 0.05 for primary and label 0.1 exploratory")
    return [(row["gene_id"], 1.0 if float(row["logFC"]) >= 0 else -1.0) for row in _rows(results) if float(row["FDR"]) <= fdr]


def _candidate_ranking(grouped: dict[tuple[str, str], list[dict[str, Any]]], min_auc: float, min_sensitivity: float, min_specificity: float) -> list[dict[str, Any]]:
    ranking: list[dict[str, Any]] = []
    for (discovery, gene_id), scores in grouped.items():
        if len(scores) != 2:
            continue
        scores = sorted(scores, key=lambda score: str(score["cohort"]))
        row = {
            "gene_id": gene_id,
            "discovery": discovery,
            "direction": scores[0]["direction"],
            "min_auc": min(float(score["auc"]) for score in scores),
            "mean_auc": mean(float(score["auc"]) for score in scores),
            "min_sensitivity": min(float(score["sensitivity"]) for score in scores),
            "min_specificity": min(float(score["specificity"]) for score in scores),
            "passed": all(float(score["auc"]) >= min_auc and float(score["sensitivity"]) >= min_sensitivity and float(score["specificity"]) >= min_specificity for score in scores),
        }
        for index, score in enumerate(scores, 1):
            for field in ("cohort", "auc", "sensitivity", "specificity"):
                row[f"validation_{index}_{field}"] = score[field]
        ranking.append(row)
    ranking.sort(key=lambda row: (float(row["min_auc"]), float(row["mean_auc"]), float(row["min_sensitivity"]), float(row["min_specificity"])), reverse=True)
    for index, row in enumerate(ranking, 1):
        row["rank"] = index
    return ranking


def run_pipeline(config_path: Path, *, include_external: bool = True) -> dict[str, Any]:
    """Run development validation and, when requested, the held-out evaluation."""
    config = json.loads(config_path.read_text(encoding="utf-8"))
    development = [Cohort(item["name"], Path(item["counts"]), Path(item["samples"])) for item in config["development"]]
    if len(development) != 3:
        raise ValueError("config requires exactly three development cohorts")
    external = Cohort(**{key: Path(value) if key in {"counts", "samples"} else value for key, value in config["external"].items()})
    development_names = [cohort.name for cohort in development]
    if len(set(development_names)) != len(development_names):
        raise ValueError("development cohort names must be unique")
    if external.name in set(development_names):
        raise ValueError("external cohort must be held out from development")
    output_dir = Path(config.get("output_dir", "outputs"))
    validation_output_dir = Path(config.get("validation_output_dir", output_dir))
    fdr = float(config.get("fdr", 0.05))
    with ThreadPoolExecutor(max_workers=3) as executor:
        list(executor.map(
            lambda cohort: _cohort_output(output_dir, cohort.name) if _edger_complete(_cohort_output(output_dir, cohort.name)) else run_edger(cohort, _cohort_output(output_dir, cohort.name)),
            development,
        ))
    reports: list[dict[str, Any]] = []
    skipped_candidates: list[dict[str, str]] = []
    candidates_by_discovery = {cohort.name: _candidates(_cohort_output(output_dir, cohort.name) / "edger_results.csv", fdr) for cohort in development}
    development_values = {cohort.name: _logcpm(_cohort_output(output_dir, cohort.name) / "logcpm.csv") for cohort in development}
    development_groups = {cohort.name: _groups(cohort.samples) for cohort in development}
    development_ensembl = {name: {identifier.split(".", 1)[0]: identifier for identifier in values} for name, values in development_values.items()}
    candidate_symbols = sorted({gene for candidates in candidates_by_discovery.values() for gene, _ in candidates if not gene.startswith("ENSG")})
    candidate_ensembl = sorted({gene.split(".", 1)[0] for candidates in candidates_by_discovery.values() for gene, _ in candidates if gene.startswith("ENSG")})
    symbol_to_ensembl = _ensembl_for_symbols(candidate_symbols) if candidate_symbols else {}
    ensembl_to_symbol = _symbols_for_ensembl(candidate_ensembl) if candidate_ensembl else {}
    for discovery in development:
        validations = [cohort for cohort in development if cohort != discovery]
        candidates = candidates_by_discovery[discovery.name]
        for gene, direction in candidates:
            validation_genes = {validation.name: _validation_gene(gene, development_values[validation.name], symbol_to_ensembl, ensembl_to_symbol, development_ensembl[validation.name]) for validation in validations}
            missing = [validation.name for validation in validations if validation_genes[validation.name] is None]
            if missing:
                skipped_candidates.append({"gene_id": gene, "discovery": discovery.name, "missing_from": ",".join(missing)})
                continue
            reports.extend(score_gene(gene, direction, discovery, validation, output_dir, str(validation_genes[validation.name]), development_values[discovery.name], development_values[validation.name], development_groups[discovery.name], development_groups[validation.name]) for validation in validations)
    min_auc = float(config.get("min_validation_auc", 0.7))
    min_sensitivity = float(config.get("min_validation_sensitivity", 0.0))
    min_specificity = float(config.get("min_validation_specificity", 0.0))
    if not 0 < min_auc <= 1 or not 0 <= min_sensitivity <= 1 or not 0 <= min_specificity <= 1:
        raise ValueError("validation thresholds must be between 0 and 1")
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in reports:
        grouped.setdefault((str(row["discovery"]), str(row["gene_id"])), []).append(row)
    candidate_ranking = _candidate_ranking(grouped, min_auc, min_sensitivity, min_specificity)
    locked_candidates = [
        scores[0]
        for scores in grouped.values()
        if len(scores) == 2 and all(
            float(score["auc"]) >= min_auc
            and float(score["sensitivity"]) >= min_sensitivity
            and float(score["specificity"]) >= min_specificity
            for score in scores
        )
    ]
    external_reports = []
    external_symbol_to_ensembl: dict[str, str] = {}
    if include_external:
        if not _edger_complete(_cohort_output(output_dir, external.name)):
            run_edger(external, _cohort_output(output_dir, external.name))
        external_values = _logcpm(_cohort_output(output_dir, external.name) / "logcpm.csv")
        external_groups = _groups(external.samples)
        missing_symbols = [str(row["gene_id"]) for row in locked_candidates if row["gene_id"] not in external_values]
        external_symbol_to_ensembl = _ensembl_for_symbols(missing_symbols)
        ensembl_to_external = {gene_id.split(".", 1)[0]: gene_id for gene_id in external_values}
        for row in locked_candidates:
            validation_gene = str(row["gene_id"])
            if validation_gene not in external_values:
                validation_gene = ensembl_to_external.get(external_symbol_to_ensembl.get(validation_gene, ""), "")
            if not validation_gene:
                skipped_candidates.append({"gene_id": row["gene_id"], "discovery": str(row["discovery"]), "missing_from": external.name})
                continue
            discovery = next(cohort for cohort in development if cohort.name == row["discovery"])
            external_reports.append(score_gene(row["gene_id"], float(row["direction"]), discovery, external, output_dir, validation_gene, development_values[discovery.name], external_values, development_groups[discovery.name], external_groups))
    report = {"development_validation": reports, "candidate_ranking": candidate_ranking, "locked_candidates": locked_candidates, "external_validation": external_reports, "external_evaluated": include_external, "skipped_candidates": skipped_candidates, "gene_id_standardization": {"source": "Ensembl REST lookup/symbol and lookup/id", "internal_symbol_mappings": len(symbol_to_ensembl), "internal_id_mappings": len(ensembl_to_symbol), "external_symbol_mappings": len(external_symbol_to_ensembl)}, "fdr": fdr, "min_validation_auc": min_auc, "min_validation_sensitivity": min_sensitivity, "min_validation_specificity": min_specificity, "validation_engine": "scikit-learn" if roc_auc_score is not None else "builtin-fallback", "validation_metrics": ["roc_auc_score", "confusion_matrix"]}
    validation_output_dir.mkdir(parents=True, exist_ok=True)
    (validation_output_dir / "validation_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report
