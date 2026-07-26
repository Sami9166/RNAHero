"""Deterministic evidence packet and PCA/heatmap artifacts for locked biomarkers."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path
from typing import Any

from workflow import _rscript


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _top_five(output_dir: Path) -> tuple[str, list[dict[str, Any]], dict[str, Any]]:
    analysis_dir = output_dir / "analysis"
    report = json.loads((analysis_dir / "validation_report.json").read_text(encoding="utf-8"))
    fdr = float(report["fdr"])
    candidates: list[dict[str, Any]] = []
    for candidate in report.get("candidate_ranking", []):
        if not candidate.get("passed"):
            continue
        result = next(
            (row for row in _rows(analysis_dir / str(candidate["discovery"]) / "edger_results.csv")
             if row["gene_id"] == candidate["gene_id"] and float(row["FDR"]) <= fdr),
            None,
        )
        if result:
            candidates.append({**candidate, "internal_logFC": float(result["logFC"])})
    if not candidates:
        raise ValueError("no internally locked FDR-significant candidates are available for summarization")
    top_five = sorted(candidates, key=lambda row: abs(row["internal_logFC"]), reverse=True)[:5]
    return str(top_five[0]["discovery"]), top_five, report


def create_pca_heatmap(output_dir: Path, cohort_id: str, gene_ids: list[str]) -> dict[str, Any]:
    """Create base-R PCA and heatmap images for selected genes in one cohort."""
    output_dir = output_dir.resolve()
    root = Path.cwd().resolve()
    if root not in output_dir.parents:
        raise ValueError("output_dir must be inside the RNAHero project")
    analysis_dir = output_dir / "analysis"
    available = {row["gene_id"] for row in _rows(analysis_dir / cohort_id / "logcpm.csv")}
    genes = [gene for gene in dict.fromkeys(gene_ids) if gene in available]
    if len(genes) < 2:
        return {"cohort": cohort_id, "status": "skipped", "reason": "at least two selected genes are required"}
    figures = analysis_dir / "summary" / "figures" / cohort_id
    gene_file = figures / "selected_genes.txt"
    figures.mkdir(parents=True, exist_ok=True)
    gene_file.write_text("\n".join(genes) + "\n", encoding="utf-8")
    result = subprocess.run(
        [
            _rscript(), str(Path(__file__).with_name("visualize.R")),
            str(analysis_dir / cohort_id / "logcpm.csv"), str(output_dir / "samples" / f"{cohort_id}.csv"),
            str(gene_file), str(figures),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "R visualization failed")
    return {
        "cohort": cohort_id,
        "status": "complete",
        "gene_ids": genes,
        "pca": str(figures / "top5_pca.png"),
        "heatmap": str(figures / "top5_heatmap.png"),
    }


def run_summarizer(output_dir: Path) -> dict[str, Any]:
    """Prepare locked evidence and visual artifacts for the ADK Markdown writer."""
    output_dir = output_dir.resolve()
    primary, top_five, report = _top_five(output_dir)
    external = json.loads((output_dir / "analysis" / "analysis_config.json").read_text(encoding="utf-8"))["external"]["name"]
    external_scores = {
        (str(row["gene_id"]), str(row["discovery"])): row
        for row in report.get("external_validation", [])
    }
    for row in top_five:
        row["external"] = external_scores.get((str(row["gene_id"]), str(row["discovery"])))
    development_visual = create_pca_heatmap(output_dir, primary, [str(row["gene_id"]) for row in top_five])
    external_visual = create_pca_heatmap(
        output_dir,
        external,
        [str(row["external"]["validation_gene_id"]) for row in top_five if row.get("external")],
    )
    critic_path = output_dir / "analysis" / "critic_report.json"
    critic = json.loads(critic_path.read_text(encoding="utf-8")) if critic_path.is_file() else {}
    return {
        "selection_basis": "Internally locked, FDR-significant candidates ranked by absolute logFC across development-cohort discoveries; external data are not used for ranking.",
        "primary_development_cohort": primary,
        "top_five": top_five,
        "visualizations": {"development": development_visual, "external": external_visual},
        "critic": critic,
    }
