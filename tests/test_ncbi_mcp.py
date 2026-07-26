import gzip
import unittest
from unittest.mock import patch

from ncbi_mcp import (
    geo_series_root,
    geo_soft_http_url,
    geo_supplement_http_url,
    likely_single_cell,
    match_count_columns_to_samples,
    parse_sample_metadata,
    raw_count_filenames,
    read_family_soft,
    search_query,
)


class NcbiMcpTests(unittest.TestCase):
    @patch("ncbi_mcp.request_bytes", return_value=gzip.compress(b"^SERIES = GSE1\n"))
    def test_family_soft_is_cached_within_one_run(self, request) -> None:
        read_family_soft.cache_clear()
        self.assertEqual(read_family_soft("GSE1"), "^SERIES = GSE1\n")
        self.assertEqual(read_family_soft("GSE1"), "^SERIES = GSE1\n")
        self.assertEqual(request.call_count, 1)

    def test_geo_series_root(self) -> None:
        self.assertEqual(geo_series_root("GSE164073"), "GSE164nnn/GSE164073")
        self.assertIn("acc=GSE164073", geo_soft_http_url("GSE164073"))
        self.assertIn("file=matrix.txt", geo_supplement_http_url("GSE164073", "matrix.txt"))

    def test_raw_count_filter_excludes_tpm_and_archives(self) -> None:
        soft = """!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/GSE1_raw_counts.tsv.gz
!Series_supplementary_file = GSE1_TPM.tsv.gz
!Series_supplementary_file = GSE1_circRNA.count.txt.gz
!Series_supplementary_file = GSE1_RAW.tar
^SAMPLE = GSM1
"""
        self.assertEqual(raw_count_filenames(soft), ["GSE1_raw_counts.tsv.gz"])

    def test_sample_metadata_keeps_all_sample_fields_for_ai_review(self) -> None:
        soft = """^SAMPLE = GSM1
!Sample_title = Tumor 1
!Sample_source_name_ch1 = lung tissue
!Sample_characteristics_ch1 = disease state: adenocarcinoma
!Sample_treatment_protocol_ch1 = none
^SAMPLE = GSM2
!Sample_title = Control 1
!Sample_characteristics_ch1 = disease state: healthy
"""
        samples = parse_sample_metadata(soft)
        self.assertEqual(samples[0]["accession"], "GSM1")
        self.assertIn("disease state: adenocarcinoma", samples[0]["characteristics"])
        self.assertIn({"field": "Sample_treatment_protocol_ch1", "value": "none"}, samples[0]["fields"])
        self.assertEqual(samples[1]["title"], "Control 1")

    def test_matches_matrix_columns_to_explicit_library_name(self) -> None:
        soft = """^SAMPLE = GSM1
!Sample_title = DMSO control
!Sample_description = Library name: E_306_1
!Sample_description = Column name in counts: E_306_1
^SAMPLE = GSM2
!Sample_title = DEX treatment
!Sample_description = Library name: E_306_DEX_1
^SAMPLE = GSM3
!Sample_title = NYU1019-N (RNA-seq)
^SAMPLE = GSM4
!Sample_title = Patient-0001
^SAMPLE = GSM5
!Sample_title = FHTCA for RNA-seq
"""
        matches = match_count_columns_to_samples(["E_306_1", "E_306_DEX_1", "NYU1019-N", "raw_Patient-0001_count", "C-FHTCA", "unknown"], parse_sample_metadata(soft))
        self.assertEqual(matches[0]["accession"], "GSM1")
        self.assertEqual(matches[1]["accession"], "GSM2")
        self.assertEqual(matches[2]["accession"], "GSM3")
        self.assertEqual(matches[3]["accession"], "GSM4")
        self.assertEqual(matches[3]["evidence"], "unique partial GEO alias")
        self.assertEqual(matches[4]["accession"], "GSM5")
        self.assertEqual(matches[5]["evidence"], "no exact GEO alias")

    def test_search_query_keeps_required_geo_filters(self) -> None:
        query = search_query("lung adenocarcinoma", "Homo sapiens")
        self.assertIn("[DataSet Type]", query)
        self.assertIn("[Organism]", query)
        self.assertIn("gse[Entry Type]", query)

    def test_single_cell_screening_rejects_explicit_studies(self) -> None:
        self.assertTrue(likely_single_cell({"title": "10x single-cell lung atlas"}))
        self.assertFalse(likely_single_cell({"title": "Bulk RNA-seq of lung tumor"}))

if __name__ == "__main__":
    unittest.main()
