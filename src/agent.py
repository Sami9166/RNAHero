"""ADK agents used by RNAHero."""

from __future__ import annotations

import asyncio
import json
import os
import sys
from pathlib import Path
from typing import Any

from dotenv import load_dotenv


load_dotenv(Path(__file__).resolve().parent.parent / ".env")


INSTRUCTION = """You are RNAHero's study-selection orchestrator.
For a disease name in Korean or English, use the GEO MCP tools to find candidate bulk RNA-seq studies.
Call sample metadata before claiming case/control groups, bulk status, or subject independence.
Return evidence and uncertainty, never a clinical diagnosis. Select three development cohorts and one held-out cohort only when the metadata supports the comparison; otherwise report the missing condition.
Do not run edgeR, change a candidate list, or tune a cutoff from held-out external-cohort results."""

CRITIC_INSTRUCTION = """You are RNAHero's final read-only pipeline critic.
You receive a compact evidence packet assembled from the completed GEO/edgeR run. Critically assess both internal validation (the two development cohorts used after each discovery cohort) and external validation (the held-out cohort).

Do not change candidates, directions, FDR, cutoffs, cohort roles, or acceptance criteria. Do not invent missing evidence. Treat candidate biomarkers as research candidates only; never make clinical-use claims.

Return JSON only with this exact shape:
{
  "verdict": "interpretable|caution|not_interpretable",
  "summary": "short Korean summary",
  "checks": [{"area":"cohort_selection|internal_validation|external_validation|reproducibility", "status":"pass|caution|fail", "evidence":"short Korean evidence"}],
  "strengths": ["..."],
  "concerns": ["..."],
  "next_steps": ["..."]
}
Use caution when evidence is incomplete, a metric is below the locked threshold, an external result is absent, or ID mapping/skipped candidates could affect interpretation."""

SUMMARIZER_INSTRUCTION = """You are RNAHero's final research summarizer.
You receive a completed deterministic biomarker summary: the internally selected top-five genes, internal and external validation metrics, PCA/heatmap artifact paths, and the critic report.
Do not re-rank genes, relax thresholds, or use external data to select a candidate. Treat every gene as a research candidate, never a clinical diagnostic.
The PCA and heatmap files have already been created. Do not call a tool during this task.
Return Korean Markdown only, beginning with "# RNAHero 요약 보고서". Include a concise conclusion, a table for the supplied five genes and their metrics, a visualization section citing the supplied PCA/heatmap PNG paths, a GO enrichment section when supplied (with its cohort-specific limitation and dot-plot path), the critic's concerns under "## 비판점", and a next step. Do not wrap the Markdown in a code fence."""


def _local_mcp_tools():
    from google.adk.tools.mcp_tool.mcp_session_manager import StdioConnectionParams
    from google.adk.tools.mcp_tool.mcp_toolset import McpToolset
    from mcp import StdioServerParameters

    return [McpToolset(connection_params=StdioConnectionParams(server_params=StdioServerParameters(
        command=sys.executable, args=["-m", "ncbi_mcp"], cwd=str(Path(__file__).parent.parent)
    )))]


def build_agent():
    """Create the GEO MCP orchestrator."""
    if not os.getenv("GEMINI_API_KEY"):
        raise RuntimeError("GEMINI_API_KEY is required to run the Gemini agent")
    from google.adk.agents import Agent

    return Agent(
        name="rnahero_orchestrator",
        model=os.getenv("GEMINI_MODEL", "gemini-3.5-flash-lite"),
        instruction=INSTRUCTION,
        tools=_local_mcp_tools(),
    )


def build_critic():
    """Create the final critic. It has no tools and cannot mutate the run."""
    if not os.getenv("GEMINI_API_KEY"):
        raise RuntimeError("GEMINI_API_KEY is required to run the Gemini agent")
    from google.adk.agents import Agent

    return Agent(
        name="pipeline_critic",
        model=os.getenv("GEMINI_MODEL", "gemini-3.5-flash-lite"),
        instruction=CRITIC_INSTRUCTION,
    )


def build_summarizer():
    """Create the ADK summarizer with access to the local visualization MCP tool."""
    if not os.getenv("GEMINI_API_KEY"):
        raise RuntimeError("GEMINI_API_KEY is required to run the Gemini agent")
    from google.adk.agents import Agent

    return Agent(
        name="rnahero_summarizer",
        model=os.getenv("GEMINI_MODEL", "gemini-3.5-flash-lite"),
        instruction=SUMMARIZER_INSTRUCTION,
        tools=_local_mcp_tools(),
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _analysis_artifact(output_dir: Path, filename: str) -> Path:
    current = output_dir / "analysis" / "validation" / filename
    legacy = output_dir / "analysis" / filename
    return current if current.is_file() else legacy if legacy.is_file() else output_dir / filename


def _provenance_artifact(output_dir: Path, filename: str) -> Path:
    current = output_dir / "provenance" / filename
    return current if current.is_file() else output_dir / filename


def _critic_packet(output_dir: Path) -> dict[str, Any]:
    """Keep the ADK prompt small while preserving every decision-relevant metric."""
    report = _read_json(_analysis_artifact(output_dir, "validation_report.json"))
    ranking = report.get("candidate_ranking", [])
    locked = [row for row in ranking if row.get("passed")] or report.get("locked_candidates", [])
    cohorts_path = _provenance_artifact(output_dir, "cohorts.json")
    cohorts = _read_json(cohorts_path) if cohorts_path.is_file() else _read_json(_provenance_artifact(output_dir, "run_manifest.json")).get("cohorts", [])
    return {
        "cohorts": cohorts,
        "analysis_config": _read_json(_provenance_artifact(output_dir, "analysis_config.json")),
        "validation": {
            "thresholds": {key: report.get(key) for key in (
                "fdr", "min_validation_auc", "min_validation_sensitivity", "min_validation_specificity"
            )},
            "internal_validation_row_count": len(report.get("development_validation", [])),
            "ranked_candidate_count": len(ranking),
            "locked_candidates": locked,
            "top_ranked_candidates": ranking[:20],
            "external_validation": report.get("external_validation", []),
            "skipped_candidate_count": len(report.get("skipped_candidates", [])),
            "gene_id_standardization": report.get("gene_id_standardization", {}),
        },
    }


async def _run_agent_async(agent: Any, app_name: str, prompt: str) -> str:
    from google.adk.runners import Runner
    from google.adk.sessions import InMemorySessionService
    from google.genai import types

    sessions = InMemorySessionService()
    session = await sessions.create_session(app_name=app_name, user_id="local")
    runner = Runner(agent=agent, app_name=app_name, session_service=sessions)
    responses: list[str] = []
    message = types.Content(role="user", parts=[types.Part(text=prompt)])
    async for event in runner.run_async(user_id="local", session_id=session.id, new_message=message):
        if event.is_final_response() and event.content:
            responses.extend(part.text for part in event.content.parts if part.text)
    if not responses:
        raise RuntimeError("ADK critic returned no final text")
    return "".join(responses)


async def _run_critic_async(prompt: str) -> str:
    return await _run_agent_async(build_critic(), "rnahero_critic", prompt)


def run_critic(output_dir: Path) -> dict[str, Any]:
    """Run one ADK critic turn against completed analysis artifacts."""
    packet = _critic_packet(output_dir)
    prompt = "Completed RNAHero run evidence packet:\n" + json.dumps(packet, ensure_ascii=False)
    text = asyncio.run(_run_critic_async(prompt)).strip()
    if text.startswith("```"):
        text = text.split("\n", 1)[1].rsplit("```", 1)[0].strip()
    result = json.loads(text)
    if not isinstance(result, dict) or "verdict" not in result or "summary" not in result:
        raise RuntimeError("ADK critic did not return the required JSON report")
    result["engine"] = "google-adk"
    return result


def run_summarizer_agent(summary: dict[str, Any]) -> str:
    """Ask the ADK summarizer to turn completed artifacts into a Markdown report."""
    prompt = "Completed deterministic summary packet:\n" + json.dumps(summary, ensure_ascii=False)
    text = asyncio.run(_run_agent_async(build_summarizer(), "rnahero_summarizer", prompt)).strip()
    if text.startswith("```"):
        text = text.split("\n", 1)[1].rsplit("```", 1)[0].strip()
    if not text.startswith("# RNAHero 요약 보고서"):
        raise RuntimeError("ADK summarizer did not return the required Markdown report")
    return text + "\n"
