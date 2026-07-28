"""Deterministic evidence packet and PCA/heatmap artifacts for locked biomarkers."""

from __future__ import annotations

import csv
import json
import ssl
import subprocess
from urllib.error import URLError
from urllib.request import Request, urlopen
from pathlib import Path
from typing import Any

import certifi

from workflow import _rscript


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["direction", "source", "term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _go_gene_sets(rows: list[dict[str, str]], fdr: float) -> dict[str, list[str]]:
    """Split FDR-significant differential genes by direction for GO enrichment."""
    return {
        "up": [row["gene_id"] for row in rows if float(row["FDR"]) <= fdr and float(row["logFC"]) > 0],
        "down": [row["gene_id"] for row in rows if float(row["FDR"]) <= fdr and float(row["logFC"]) < 0],
    }


def _gost(genes: list[str]) -> list[dict[str, Any]]:
    payload = json.dumps({
        "organism": "hsapiens",
        "query": genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": True,
    }).encode()
    request = Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", data=payload, headers={"Content-Type": "application/json"})
    with urlopen(request, timeout=30, context=ssl.create_default_context(cafile=certifi.where())) as response:
        return json.loads(response.read().decode("utf-8")).get("result", [])


def create_go_enrichment(output_dir: Path, cohort_id: str, fdr: float) -> dict[str, Any]:
    """Write GO enrichment table and dot plot for one development cohort's DEGs."""
    analysis_dir = output_dir / "analysis"
    rows = _rows(analysis_dir / "cohorts" / cohort_id / "edgeR" / "edger_results.csv")
    gene_sets = _go_gene_sets(rows, fdr)
    output = analysis_dir / "go" / cohort_id
    result_csv = output / "go_enrichment.csv"
    go_rows: list[dict[str, Any]] = []
    error = ""
    try:
        for direction, genes in gene_sets.items():
            if len(genes) < 5:
                continue
            for term in _gost(genes):
                if term.get("significant"):
                    go_rows.append({
                        "direction": direction,
                        "source": term.get("source", ""),
                        "term_id": term.get("native", ""),
                        "term_name": term.get("name", ""),
                        "p_value": term.get("p_value", ""),
                        "term_size": term.get("term_size", ""),
                        "query_size": term.get("query_size", ""),
                        "intersection_size": term.get("intersection_size", ""),
                    })
    except (URLError, TimeoutError, ValueError) as exc:
        error = str(exc)
    go_rows.sort(key=lambda row: float(row["p_value"]))
    _write_rows(result_csv, go_rows)
    output.mkdir(parents=True, exist_ok=True)
    figure = output / "go_dotplot.png"
    result = subprocess.run(
        [_rscript(), str(Path(__file__).with_name("go_visualize.R")), str(result_csv), str(figure), cohort_id],
        capture_output=True,
        text=True,
    )
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "GO visualization failed")
    return {
        "cohort": cohort_id,
        "status": "complete" if not error else "unavailable",
        "fdr": fdr,
        "input_gene_counts": {key: len(value) for key, value in gene_sets.items()},
        "term_count": len(go_rows),
        "results": str(result_csv),
        "dotplot": str(figure),
        "reason": error or None,
    }


def _top_five(output_dir: Path) -> tuple[str, list[dict[str, Any]], dict[str, Any]]:
    analysis_dir = output_dir / "analysis"
    report = json.loads((analysis_dir / "validation" / "validation_report.json").read_text(encoding="utf-8"))
    fdr = float(report["fdr"])
    candidates: list[dict[str, Any]] = []
    for candidate in report.get("candidate_ranking", []):
        if not candidate.get("passed"):
            continue
        result = next(
            (row for row in _rows(analysis_dir / "cohorts" / str(candidate["discovery"]) / "edgeR" / "edger_results.csv")
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
    available = {row["gene_id"] for row in _rows(analysis_dir / "cohorts" / cohort_id / "edgeR" / "logcpm.csv")}
    genes = [gene for gene in dict.fromkeys(gene_ids) if gene in available]
    if len(genes) < 2:
        return {"cohort": cohort_id, "status": "skipped", "reason": "at least two selected genes are required"}
    figures = output_dir / "report" / "figures" / cohort_id
    gene_file = figures / "selected_genes.txt"
    figures.mkdir(parents=True, exist_ok=True)
    gene_file.write_text("\n".join(genes) + "\n", encoding="utf-8")
    result = subprocess.run(
        [
            _rscript(), str(Path(__file__).with_name("visualize.R")),
            str(analysis_dir / "cohorts" / cohort_id / "edgeR" / "logcpm.csv"), str(output_dir / "cohorts" / cohort_id / "input" / "samples.csv"),
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
    external = json.loads((output_dir / "provenance" / "analysis_config.json").read_text(encoding="utf-8"))["external"]["name"]
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
    go_enrichment = create_go_enrichment(output_dir, primary, float(report["fdr"]))
    critic_path = output_dir / "analysis" / "validation" / "critic_report.json"
    critic = json.loads(critic_path.read_text(encoding="utf-8")) if critic_path.is_file() else {}
    return {
        "selection_basis": "Internally locked, FDR-significant candidates ranked by absolute logFC across development-cohort discoveries; external data are not used for ranking.",
        "primary_development_cohort": primary,
        "top_five": top_five,
        "visualizations": {"development": development_visual, "external": external_visual},
        "go_enrichment": go_enrichment,
        "critic": critic,
    }
