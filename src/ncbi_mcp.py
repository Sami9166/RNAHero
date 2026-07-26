"""Minimal MCP server for finding and downloading GEO bulk RNA-seq count files."""

from __future__ import annotations

import gzip
import hashlib
from http.client import IncompleteRead, RemoteDisconnected
import json
import os
import re
import threading
import time
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse
from urllib.request import Request, urlopen

from mcp.server.fastmcp import FastMCP


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series"
USER_AGENT = "RNAHero/0.1 (NCBI GEO retrieval)"
ACCESSION_RE = re.compile(r"^GSE\d+$", re.IGNORECASE)
COUNT_TERMS = ("count", "counts", "readcount", "read_count")
NON_COUNT_TERMS = ("tpm", "fpkm", "rpkm", "cpm", "normalized", "norm_", "log2", "circrna")
SINGLE_CELL_TERMS = ("single-cell", "single cell", "scrna", "sc-rna", "10x genomics", "10x")

mcp = FastMCP(
    "RNAHero NCBI GEO",
    instructions=(
        "Search public GEO bulk RNA-seq studies and retrieve submitter or "
        "NCBI raw-count matrices. Use sample metadata to classify cases and controls "
        "before running differential-expression analysis."
    ),
    json_response=True,
)

_request_lock = threading.Lock()
_last_request_started = 0.0


class NcbiError(RuntimeError):
    """Raised when NCBI cannot return a usable response."""


def normalize_accession(accession: str) -> str:
    accession = accession.strip().upper()
    if not ACCESSION_RE.fullmatch(accession):
        raise ValueError("accession must look like GSE12345")
    return accession


def geo_series_root(accession: str) -> str:
    """Return the GEO FTP directory that contains one GSE accession."""
    accession = normalize_accession(accession)
    digits = accession[3:]
    return f"GSE{digits[:-3]}nnn/{accession}"


def geo_soft_url(accession: str) -> str:
    root = geo_series_root(accession)
    return f"{GEO_FTP_URL}/{root}/soft/{accession}_family.soft.gz"


def geo_soft_http_url(accession: str) -> str:
    return "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + urlencode({"acc": normalize_accession(accession), "targ": "all", "view": "full", "form": "text"})


def geo_supplement_url(accession: str, filename: str) -> str:
    root = geo_series_root(accession)
    safe_name = Path(filename).name
    if safe_name != filename:
        raise ValueError("filename must not include a path")
    return f"{GEO_FTP_URL}/{root}/suppl/{safe_name}"


def geo_supplement_http_url(accession: str, filename: str) -> str:
    safe_name = Path(filename).name
    if safe_name != filename:
        raise ValueError("filename must not include a path")
    return "https://www.ncbi.nlm.nih.gov/geo/download/?" + urlencode({"acc": normalize_accession(accession), "format": "file", "file": safe_name})


def request_bytes(url: str) -> bytes:
    global _last_request_started
    interval = 0.1 if url.startswith(EUTILS_URL) and _setting("NCBI_API_KEY") else 0.35
    # ponytail: process-wide rate limit; use a shared limiter only for multi-process MCP traffic.
    with _request_lock:
        wait = interval - (time.monotonic() - _last_request_started)
        if wait > 0:
            time.sleep(wait)
        _last_request_started = time.monotonic()
    request = Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return response.read()
        except HTTPError as error:
            if error.code in {429, 500, 502, 503} and attempt < 2:
                time.sleep(2**attempt)
                continue
            raise NcbiError(f"NCBI returned HTTP {error.code} for {url}") from error
        except URLError as error:
            raise NcbiError(f"Could not reach NCBI: {error.reason}") from error
        except (IncompleteRead, RemoteDisconnected, ConnectionResetError):
            if attempt < 2:
                time.sleep(2**attempt)
                continue
            raise NcbiError(f"NCBI response was incomplete for {url}")
    raise NcbiError(f"Could not reach NCBI: retries exhausted for {url}")


def request_json(path: str, **params: str | int) -> dict[str, Any]:
    api_key = _setting("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    params.setdefault("tool", "RNAHero")
    email = _setting("NCBI_EMAIL")
    if email:
        params.setdefault("email", email)
    url = f"{EUTILS_URL}/{path}?{urlencode(params)}"
    try:
        return json.loads(request_bytes(url))
    except json.JSONDecodeError as error:
        raise NcbiError("NCBI returned invalid JSON") from error


def _setting(name: str) -> str:
    """Read an optional local setting without adding a dotenv dependency."""
    if value := os.getenv(name):
        return value
    env_file = Path.cwd() / ".env"
    if not env_file.is_file():
        return ""
    for line in env_file.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip() == name:
            return value.strip().strip('"').strip("'")
    return ""


@lru_cache(maxsize=128)
def read_family_soft(accession: str) -> str:
    try:
        return gzip.decompress(request_bytes(geo_soft_url(accession))).decode("utf-8", errors="replace")
    except NcbiError:
        try:
            return request_bytes(geo_soft_http_url(accession)).decode("utf-8", errors="replace")
        except NcbiError:
            raise
    except OSError as error:
        raise NcbiError(f"Could not decompress the GEO SOFT file for {accession}") from error


def raw_count_filenames(soft_text: str) -> list[str]:
    """Extract likely raw-count supplementary files from a GEO family SOFT record."""
    filenames: list[str] = []
    for line in soft_text.splitlines():
        if line.startswith("^SAMPLE"):
            break
        if line.startswith("!Series_supplementary_file ="):
            location = line.split("=", 1)[1].strip()
            filename = Path(urlparse(location).path).name
            lowered = filename.lower()
            if (
                any(term in lowered for term in COUNT_TERMS)
                and not any(term in lowered for term in NON_COUNT_TERMS)
                and not lowered.endswith((".tar", ".tar.gz", ".tgz"))
            ):
                filenames.append(filename)
    return sorted(set(filenames))


def parse_sample_metadata(soft_text: str, max_samples: int = 200) -> list[dict[str, Any]]:
    """Return complete submitter-supplied sample fields for AI review."""
    samples: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None

    for line in soft_text.splitlines():
        if line.startswith("^SAMPLE ="):
            if current:
                samples.append(current)
                if len(samples) >= max_samples:
                    break
            current = {"accession": line.split("=", 1)[1].strip(), "characteristics": [], "descriptions": [], "fields": []}
            continue
        if current is None or not line.startswith("!") or " =" not in line:
            continue

        key, value = (part.strip() for part in line[1:].split("=", 1))
        if key.startswith("Sample_"):
            current["fields"].append({"field": key, "value": value})
        if key == "Sample_title":
            current["title"] = value
        elif key.startswith("Sample_source_name"):
            current["source_name"] = value
        elif key.startswith("Sample_organism"):
            current["organism"] = value
        elif key.startswith("Sample_characteristics"):
            current["characteristics"].append(value)
        elif key.startswith("Sample_description"):
            current["descriptions"].append(value)

    if current and len(samples) < max_samples:
        samples.append(current)
    return samples


def match_count_columns_to_samples(raw_columns: list[str], samples: list[dict[str, Any]]) -> list[dict[str, str]]:
    """Match matrix headers only to explicit GEO sample identifiers or library names."""
    aliases: dict[str, list[tuple[dict[str, Any], str]]] = {}
    for sample in samples:
        values = [(sample["accession"], "GEO accession"), (sample.get("title", ""), "GEO sample title")]
        title = sample.get("title", "")
        if title.lower().endswith("(rna-seq)"):
            values.append((title[: -len("(rna-seq)")].strip(), "GEO sample title without RNA-seq suffix"))
        if title.lower().endswith("for rna-seq"):
            values.append((title[: -len("for rna-seq")].strip(), "GEO sample title without RNA-seq suffix"))
        for description in sample.get("descriptions", []):
            if ":" in description:
                values.append((description.split(":", 1)[1].strip(), f"GEO Sample_description: {description}"))
        for value, evidence in values:
            if value:
                aliases.setdefault(value.casefold(), []).append((sample, evidence))
    matches: list[dict[str, str]] = []
    for column in raw_columns:
        candidates = aliases.get(column.casefold(), [])
        evidence = ""
        if not candidates and re.match(r"^[A-Za-z][-_]", column):
            candidates = aliases.get(column[2:].casefold(), [])
            evidence = "GEO alias without count-column prefix" if candidates else ""
        if not candidates:
            partial = [
                pair
                for alias, pairs in aliases.items()
                if len(alias) >= 4 and (alias in column.casefold() or column.casefold() in alias)
                for pair in pairs
            ]
            unique_partial = {sample["accession"]: (sample, source) for sample, source in partial}
            if len(unique_partial) == 1:
                candidates = list(unique_partial.values())
                evidence = "unique partial GEO alias"
        unique = {sample["accession"]: (sample, evidence) for sample, evidence in candidates}
        if len(unique) == 1:
            sample, source = next(iter(unique.values()))
            matches.append({"raw_column": column, "accession": sample["accession"], "evidence": evidence or source})
        elif not unique:
            matches.append({"raw_column": column, "accession": "", "evidence": "no exact GEO alias"})
        else:
            matches.append({"raw_column": column, "accession": "", "evidence": "ambiguous GEO alias"})
    return matches


def search_query(disease: str, organism: str) -> str:
    return (
        f'({disease}) AND "expression profiling by high throughput sequencing"[DataSet Type] '
        f'AND {organism}[Organism] AND gse[Entry Type]'
    )


def likely_single_cell(study: dict[str, Any]) -> bool:
    text = f"{study.get('title', '')} {study.get('summary', '')}".lower()
    return any(term in text for term in SINGLE_CELL_TERMS)


def search_gse_summaries(disease: str, organism: str, limit: int) -> list[dict[str, Any]]:
    search = request_json(
        "esearch.fcgi",
        db="gds",
        term=search_query(disease, organism),
        retmode="json",
        retmax=max(limit * 4, 20),
    )
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    summary = request_json("esummary.fcgi", db="gds", id=",".join(ids), retmode="json")
    result = summary.get("result", {})
    studies: list[dict[str, Any]] = []
    for uid in result.get("uids", []):
        record = result[uid]
        accession = record.get("accession", "")
        if accession.upper().startswith("GSE"):
            studies.append(record)
    return studies


def find_bulk_rnaseq(
    disease: str,
    organism: str = "Homo sapiens",
    limit: int = 10,
) -> list[dict[str, Any]]:
    """Find likely bulk GEO RNA-seq GSE studies that expose a likely raw-count file.

    `disease` must be an English disease name or NCBI-ready synonym query. This tool
    verifies RNA-seq/GSE/organism/count-file conditions and rejects studies that are
    explicitly labelled single-cell. It deliberately does not infer case-control labels
    or guarantee bulk status; call `get_geo_sample_metadata` for that evidence.
    """
    if not disease.strip():
        raise ValueError("disease is required")
    if not 1 <= limit <= 25:
        raise ValueError("limit must be between 1 and 25")

    matches: list[dict[str, Any]] = []
    for study in search_gse_summaries(disease, organism, limit):
        if likely_single_cell(study):
            continue
        accession = normalize_accession(study["accession"])
        count_files = raw_count_filenames(read_family_soft(accession))
        if not count_files:
            continue
        matches.append(
            {
                "gse_id": accession,
                "title": study.get("title", ""),
                "summary": study.get("summary", ""),
                "organism": study.get("taxon", organism),
                "sample_count": study.get("n_samples"),
                "raw_count_files": [
                    {"filename": name, "url": geo_supplement_url(accession, name)} for name in count_files
                ],
                "retrieval_checks": ["rna_seq_gse", "not_explicitly_single_cell", "raw_count_filename"],
                "requires_metadata_check": ["bulk_status", "case_control_groups", "independent_subjects"],
            }
        )
        if len(matches) == limit:
            break
    return matches


@mcp.tool()
def search_bulk_rnaseq(
    disease: str,
    organism: str = "Homo sapiens",
    limit: int = 10,
) -> list[dict[str, Any]]:
    """Find likely bulk GEO RNA-seq GSE studies that expose a likely raw-count file."""
    return find_bulk_rnaseq(disease, organism, limit)


def fetch_geo_sample_metadata(gse_id: str, max_samples: int = 200) -> dict[str, Any]:
    """Get GEO sample titles and characteristics for AI-assisted case/control assignment."""
    gse_id = normalize_accession(gse_id)
    if not 1 <= max_samples <= 500:
        raise ValueError("max_samples must be between 1 and 500")
    soft_text = read_family_soft(gse_id)
    samples = parse_sample_metadata(soft_text, max_samples=max_samples)
    return {
        "gse_id": gse_id,
        "sample_metadata": samples,
        "truncated": len(parse_sample_metadata(soft_text, max_samples=max_samples + 1)) > max_samples,
    }


@mcp.tool()
def get_geo_sample_metadata(gse_id: str, max_samples: int = 200) -> dict[str, Any]:
    """Get GEO sample titles and characteristics for AI-assisted case/control assignment."""
    return fetch_geo_sample_metadata(gse_id, max_samples)


@mcp.tool()
def match_geo_sample_columns(gse_id: str, raw_columns: list[str]) -> list[dict[str, str]]:
    """Match raw-count column names to explicit GEO accessions or library names."""
    return match_count_columns_to_samples(raw_columns, fetch_geo_sample_metadata(gse_id)["sample_metadata"])


@mcp.tool()
def assess_geo_clinical_cohort(gse_id: str, max_samples: int = 200) -> dict[str, Any]:
    """Return conservative sample-level case/control evidence before raw-count download."""
    return fetch_geo_sample_metadata(gse_id, max_samples)


def fetch_raw_count_matrix(gse_id: str, filename: str, output_dir: Path | None = None) -> dict[str, Any]:
    """Download one previously discovered raw-count file into RNAHero/data/raw.

    The requested filename must be a likely raw-count supplementary file listed in the
    same GSE SOFT record; arbitrary URLs and paths are rejected.
    """
    gse_id = normalize_accession(gse_id)
    count_files = raw_count_filenames(read_family_soft(gse_id))
    if filename not in count_files:
        raise ValueError("filename is not a verified raw-count candidate for this GSE")

    destination = (output_dir or Path.cwd() / "data" / "raw") / Path(filename).name
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        data = request_bytes(geo_supplement_url(gse_id, filename))
    except NcbiError:
        data = request_bytes(geo_supplement_http_url(gse_id, filename))
    destination.write_bytes(data)
    return {
        "gse_id": gse_id,
        "filename": filename,
        "path": str(destination.resolve()),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


@mcp.tool()
def download_raw_count_matrix(gse_id: str, filename: str) -> dict[str, Any]:
    """Download one previously discovered raw-count file into RNAHero/data/raw."""
    return fetch_raw_count_matrix(gse_id, filename)


@mcp.tool()
def create_pca_heatmap(output_dir: str, cohort_id: str, gene_ids: list[str]) -> dict[str, Any]:
    """Create PCA and heatmap PNGs for selected genes from a completed local RNAHero run."""
    from summarizer import create_pca_heatmap as create_visuals

    return create_visuals(Path(output_dir), cohort_id, gene_ids)


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
