"""Small command-line surface for the GEO retrieval and edgeR workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ncbi_mcp import fetch_geo_sample_metadata, fetch_raw_count_matrix, find_bulk_rnaseq
from agent import run_summarizer_agent
from run import resume_analysis, run_disease
from summarizer import run_summarizer
from web import serve
from workflow import run_pipeline


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(prog="rnahero")
    commands = parser.add_subparsers(dest="command", required=True)
    search = commands.add_parser("search", help="find GEO RNA-seq studies with likely count matrices")
    search.add_argument("disease")
    search.add_argument("--organism", default="Homo sapiens")
    search.add_argument("--limit", type=int, default=4)
    metadata = commands.add_parser("metadata", help="get sample titles and characteristics")
    metadata.add_argument("gse_id")
    metadata.add_argument("--max-samples", type=int, default=200)
    download = commands.add_parser("download", help="download a verified raw-count candidate")
    download.add_argument("gse_id")
    download.add_argument("filename")
    download.add_argument("--output-dir", type=Path, default=Path("data/raw"))
    pipeline = commands.add_parser("pipeline", help="run 3 development cohorts and one held-out cohort")
    pipeline.add_argument("config", type=Path)
    run = commands.add_parser("run", help="search, save, select cohorts, and analyse in one command")
    run.add_argument("disease")
    run.add_argument("--output", type=Path, default=Path("output"))
    run.add_argument("--organism", default="Homo sapiens")
    run.add_argument("--limit", type=int, default=10)
    resume = commands.add_parser("resume", help="finish a retained analysis without repeating GEO collection or completed edgeR")
    resume.add_argument("--output", type=Path, default=Path("output"))
    summarize = commands.add_parser("summarize", help="write a top-logFC report and PCA/heatmap artifacts from a completed run")
    summarize.add_argument("--output", type=Path, default=Path("output"))
    web = commands.add_parser("web", help="open the local RNAHero browser interface")
    web.add_argument("--port", type=int, default=8000)

    args = parser.parse_args()
    if args.command == "search":
        result = find_bulk_rnaseq(args.disease, args.organism, args.limit)
    elif args.command == "metadata":
        result = fetch_geo_sample_metadata(args.gse_id, args.max_samples)
    elif args.command == "download":
        result = fetch_raw_count_matrix(args.gse_id, args.filename, args.output_dir)
    elif args.command == "pipeline":
        result = run_pipeline(args.config)
    elif args.command == "run":
        result = run_disease(args.disease, args.output, args.organism, args.limit)
    elif args.command == "resume":
        result = resume_analysis(args.output)
    elif args.command == "summarize":
        summary = run_summarizer(args.output)
        report_path = args.output / "report" / "summary_report.md"
        report_path.write_text(run_summarizer_agent(summary), encoding="utf-8")
        result = {"report": str(report_path), "figures": summary["visualizations"]}
    else:
        serve(args.port)
        return
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
