import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from run import _ai_sample_metadata, _assess_cohort, _columns, _groups_from_assessment, _prepare_matrix, _retry_seconds, clear_output, run_disease


class RunTests(unittest.TestCase):
    def test_clear_output_only_removes_a_project_subfolder(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            (output / "provenance" / "run_manifest.json").parent.mkdir(parents=True)
            (output / "provenance" / "run_manifest.json").write_text("{}", encoding="utf-8")
            with patch("run.Path.cwd", return_value=root):
                self.assertEqual(clear_output(output)["stage"], "cleared")
                with self.assertRaises(ValueError):
                    clear_output(root)
            self.assertFalse(output.exists())

    @patch("run._assess_cohorts", return_value={"GSE1": {"status": "insufficient_usable_groups", "case_count": 0, "control_count": 0, "samples": []}})
    @patch("run.fetch_geo_sample_metadata")
    @patch("run.find_bulk_rnaseq")
    def test_run_persists_search_metadata_before_safe_stop(self, search, metadata, _assessment) -> None:
        search.return_value = [{"gse_id": "GSE1", "title": "study", "sample_count": 4, "raw_count_files": []}]
        metadata.return_value = {"gse_id": "GSE1", "sample_metadata": [{"accession": "GSM1", "fields": []}]}
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            result = run_disease("disease", output, limit=1)
            self.assertEqual(result["stage"], "needs_cohort_review")
            self.assertTrue((output / "provenance" / "search_results.json").is_file())
            self.assertFalse((output / "cohorts").exists())
            self.assertFalse((output / "evidence").exists())
            self.assertTrue((output / "provenance" / "run_manifest.json").is_file())

    @patch("run._assess_cohorts", return_value={"GSE1": {"status": "assessment_failed", "reason": "API key not valid", "samples": []}})
    @patch("run.fetch_geo_sample_metadata")
    @patch("run.find_bulk_rnaseq")
    def test_run_reports_ai_configuration_failure(self, search, metadata, _assessment) -> None:
        search.return_value = [{"gse_id": "GSE1", "title": "study", "sample_count": 4, "raw_count_files": []}]
        metadata.return_value = {"gse_id": "GSE1", "sample_metadata": []}
        with tempfile.TemporaryDirectory() as temporary:
            result = run_disease("disease", Path(temporary), limit=1)
        self.assertEqual(result["stage"], "cohort_assessment_failed")
        self.assertEqual(result["reason"], "API key not valid")

    @patch("run.find_bulk_rnaseq")
    def test_clinical_candidates_prefers_strict_then_fills_without_duplicates(self, search) -> None:
        search.side_effect = [
            [{"gse_id": "GSE1"}],
            [{"gse_id": "GSE1"}, {"gse_id": "GSE2"}],
        ]
        from run import _clinical_candidates

        queries, studies = _clinical_candidates("disease", "Homo sapiens", 2)
        self.assertEqual(len(queries), 2)
        self.assertEqual([study["gse_id"] for study in studies], ["GSE1", "GSE2"])

    def test_prepare_matrix_handles_geo_matrix_without_gene_header(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw.tsv"
            raw.write_text("c1\tc2\tk1\tk2\nGENE1\t1\t2\t10\t11\n", encoding="utf-8")
            prepared = root / "prepared.csv"
            groups = [{"sample_id": sample, "group": "case"} for sample in ("c1", "c2", "k1", "k2")]
            self.assertEqual(_columns(raw), ["c1", "c2", "k1", "k2"])
            _prepare_matrix(raw, groups, prepared)
            self.assertEqual(prepared.read_text(encoding="utf-8").splitlines(), ["gene_id,c1,c2,k1,k2", "GENE1,1,2,10,11"])

    def test_ai_metadata_omits_repeated_contact_and_processing_fields(self) -> None:
        compact = _ai_sample_metadata([{"accession": "GSM1", "fields": [
            {"field": "Sample_title", "value": "tumor"},
            {"field": "Sample_characteristics_ch1", "value": "tissue: lung"},
            {"field": "Sample_contact_email", "value": "person@example.org"},
            {"field": "Sample_data_processing", "value": "x" * 1000},
        ]}])
        self.assertEqual(compact, [{"accession": "GSM1", "fields": [
            {"field": "Sample_title", "value": "tumor"},
            {"field": "Sample_characteristics_ch1", "value": "tissue: lung"},
        ]}])
        self.assertEqual(_retry_seconds("retryDelay': '45s'"), 45)

    @patch("run._ask_gemini")
    def test_ai_assessment_excludes_semantically_invalid_samples(self, ask) -> None:
        samples = [
            {"accession": f"GSM{i}", "fields": [{"field": "Sample_title", "value": value}]}
            for i, value in enumerate(("tumor 1", "tumor 2", "adjacent normal 1", "adjacent normal 2", "A549 tumor"), 1)
        ]
        ask.return_value = ({"samples": [
            {"accession": "GSM1", "group": "case", "disease_match": True, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "tumor 1"}]},
            {"accession": "GSM2", "group": "case", "disease_match": True, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "tumor 2"}]},
            {"accession": "GSM3", "group": "control", "disease_match": True, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "adjacent normal 1"}]},
            {"accession": "GSM4", "group": "control", "disease_match": True, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "adjacent normal 2"}]},
            {"accession": "GSM5", "group": "case", "disease_match": True, "clinical_tissue": False, "cell_line": True, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "A549 tumor"}]},
        ]}, None)
        assessment, error = _assess_cohort("lung adenocarcinoma", {"sample_metadata": samples})
        self.assertIsNone(error)
        self.assertEqual(assessment["status"], "eligible")
        self.assertEqual((assessment["case_count"], assessment["control_count"]), (2, 2))
        self.assertEqual(assessment["samples"][-1]["excluded_by"], ["cell_line", "not_clinical_tissue"])

    @patch("run._ask_gemini")
    def test_normal_control_is_not_excluded_for_missing_disease_label(self, ask) -> None:
        samples = [
            {"accession": "GSM1", "fields": [{"field": "Sample_title", "value": "lung tumor"}]},
            {"accession": "GSM2", "fields": [{"field": "Sample_title", "value": "adjacent normal lung"}]},
            {"accession": "GSM3", "fields": [{"field": "Sample_title", "value": "adjacent normal lung 2"}]},
            {"accession": "GSM4", "fields": [{"field": "Sample_title", "value": "lung tumor 2"}]},
        ]
        ask.return_value = ({"samples": [
            {"accession": "GSM1", "group": "case", "disease_match": True, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "lung tumor"}]},
            {"accession": "GSM2", "group": "control", "disease_match": False, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "adjacent normal lung"}]},
            {"accession": "GSM3", "group": "control", "disease_match": False, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "adjacent normal lung 2"}]},
            {"accession": "GSM4", "group": "case", "disease_match": True, "clinical_tissue": True, "cell_line": False, "treatment_experiment": False, "contradictions": [], "evidence": [{"field": "Sample_title", "value": "lung tumor 2"}]},
        ]}, None)
        assessment, error = _assess_cohort("lung adenocarcinoma", {"title": "lung adenocarcinoma", "summary": "", "sample_metadata": samples})
        self.assertIsNone(error)
        self.assertEqual(assessment["status"], "eligible")
        self.assertEqual(assessment["samples"][1]["excluded_by"], [])

    def test_group_mapping_requires_unique_gsm_matches(self) -> None:
        assessment = {"samples": [
            {"accession": "GSM1", "group": "case", "included": True, "evidence": []},
            {"accession": "GSM2", "group": "case", "included": True, "evidence": []},
            {"accession": "GSM3", "group": "control", "included": True, "evidence": []},
            {"accession": "GSM4", "group": "control", "included": True, "evidence": []},
        ]}
        crosswalk = [{"raw_column": f"s{i}", "accession": f"GSM{i}"} for i in range(1, 5)]
        groups, error = _groups_from_assessment(assessment, ["s1", "s2", "s3", "s4"], crosswalk)
        self.assertIsNone(error)
        self.assertEqual([item["group"] for item in groups], ["case", "case", "control", "control"])

    @patch("run.run_critic", return_value={"engine": "google-adk", "verdict": "caution", "summary": "reviewed"})
    @patch("run.run_pipeline", return_value={"candidate_ranking": [], "fdr": 0.05, "min_validation_auc": 0.8, "min_validation_sensitivity": 0.7, "min_validation_specificity": 0.7})
    @patch("run._select_cohorts")
    @patch("run._assess_cohorts")
    @patch("run.match_count_columns_to_samples")
    @patch("run.fetch_raw_count_matrix")
    @patch("run.fetch_geo_sample_metadata")
    @patch("run.find_bulk_rnaseq")
    def test_one_command_prepares_four_cohorts_and_runs_analysis(self, search, metadata, download, crosswalk, assessment, selection, pipeline, critic) -> None:
        studies = [{"gse_id": f"GSE{i}", "title": f"study {i}", "sample_count": 4, "raw_count_files": [{"filename": "matrix.csv"}]} for i in range(1, 5)]
        search.return_value = studies
        metadata.side_effect = lambda gse: {"gse_id": gse, "sample_metadata": []}
        assessment.side_effect = lambda _disease, items: {study["gse_id"]: {"status": "eligible", "case_count": 2, "control_count": 2, "samples": [
            {"accession": "GSM1", "group": "control", "included": True, "evidence": []},
            {"accession": "GSM2", "group": "control", "included": True, "evidence": []},
            {"accession": "GSM3", "group": "case", "included": True, "evidence": []},
            {"accession": "GSM4", "group": "case", "included": True, "evidence": []},
        ]} for study in items}
        crosswalk.return_value = [
            {"raw_column": "c1", "accession": "GSM1"}, {"raw_column": "c2", "accession": "GSM2"},
            {"raw_column": "k1", "accession": "GSM3"}, {"raw_column": "k2", "accession": "GSM4"},
        ]
        selection.return_value = ([{"gse_id": "GSE1", "role": "development", "reason": "ok"}, {"gse_id": "GSE2", "role": "development", "reason": "ok"}, {"gse_id": "GSE3", "role": "development", "reason": "ok"}, {"gse_id": "GSE4", "role": "external", "reason": "ok"}], None)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)

            def write_raw(gse: str, _filename: str, directory: Path):
                directory.mkdir(parents=True, exist_ok=True)
                path = directory / "matrix.csv"
                path.write_text("gene_id,c1,c2,k1,k2\nGENE1,1,2,10,11\n", encoding="utf-8")
                return {"path": str(path)}

            download.side_effect = write_raw
            result = run_disease("disease", root, limit=4)
            self.assertEqual(result["stage"], "complete")
            self.assertTrue((root / "provenance" / "analysis_config.json").is_file())
            self.assertTrue((root / "cohorts" / "GSE4" / "raw" / "matrix.csv").is_file())
            self.assertTrue((root / "cohorts" / "GSE4" / "input" / "crosswalk.json").is_file())
            self.assertTrue((root / "analysis" / "biomarkers.json").is_file())
            self.assertTrue((root / "analysis" / "biomarkers.csv").is_file())
            self.assertTrue((root / "analysis" / "candidate_scores.csv").is_file())
            self.assertTrue((root / "analysis" / "validation" / "internal_validation.json").is_file())
            self.assertTrue((root / "analysis" / "validation" / "internal_validation.csv").is_file())
            self.assertTrue((root / "analysis" / "validation" / "external_validation.json").is_file())
            self.assertTrue((root / "analysis" / "validation" / "external_validation.csv").is_file())
            self.assertEqual(json.loads((root / "analysis" / "validation" / "critic_report.json").read_text(encoding="utf-8"))["verdict"], "caution")
            self.assertEqual(result["critic"]["engine"], "google-adk")
            self.assertEqual(json.loads((root / "provenance" / "search_results.json").read_text(encoding="utf-8")), [{"gse_id": f"GSE{i}", "title": f"study {i}"} for i in range(1, 5)])
            self.assertEqual(pipeline.call_count, 1)
            self.assertEqual(pipeline.call_args.args[0], root / "provenance" / "analysis_config.json")
            self.assertEqual(critic.call_count, 1)


if __name__ == "__main__":
    unittest.main()
