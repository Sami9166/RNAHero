import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from workflow import Cohort, _candidate_ranking, _edger_complete, _ensembl_for_symbols, _validation_gene, run_edger, score_gene


class WorkflowTests(unittest.TestCase):
    def write_csv(self, path: Path, header: list[str], rows: list[list[object]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def test_edger_writes_results_and_locked_validation_scores(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cohorts = []
            for name, shift in (("dev1", 0), ("dev2", 2)):
                counts = root / f"{name}_counts.csv"
                samples = root / f"{name}_samples.csv"
                self.write_csv(counts, ["gene_id", "c1", "c2", "c3", "k1", "k2", "k3"], [["GENE1", 10 + shift, 12 + shift, 11 + shift, 100 + shift, 105 + shift, 98 + shift], ["GENE2", 50, 52, 49, 51, 50, 52], ["GENE3", 3, 4, 3, 5, 4, 3]])
                self.write_csv(samples, ["sample_id", "group"], [["c1", "control"], ["c2", "control"], ["c3", "control"], ["k1", "case"], ["k2", "case"], ["k3", "case"]])
                cohorts.append(Cohort(name, counts, samples))

            outputs = root / "outputs"
            for cohort in cohorts:
                run_edger(cohort, outputs / cohort.name)
            score = score_gene("GENE1", 1.0, cohorts[0], cohorts[1], outputs)
            self.assertEqual(score["cohort"], "dev2")
            self.assertGreaterEqual(float(score["auc"]), 0.5)
            self.assertTrue((outputs / "dev1" / "edger_results.csv").is_file())

    def test_validation_gene_maps_symbols_and_ensembl_ids(self) -> None:
        self.assertEqual(_validation_gene("A1BG", {"ENSG00000121410": {}}, {"A1BG": "ENSG00000121410"}, {}), "ENSG00000121410")
        self.assertEqual(_validation_gene("ENSG00000121410", {"A1BG": {}}, {}, {"ENSG00000121410": "A1BG"}), "A1BG")

    def test_ensembl_mapping_reuses_local_cache(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cache = root / ".rnahero" / "ensembl_gene_id_cache.json"
            cache.parent.mkdir()
            cache.write_text('{"symbol_to_ensembl":{"A1BG":"ENSG00000121410"},"ensembl_to_symbol":{}}', encoding="utf-8")
            with patch("workflow.Path.cwd", return_value=root), patch("workflow.urlopen") as request:
                self.assertEqual(_ensembl_for_symbols(["A1BG"]), {"A1BG": "ENSG00000121410"})
            request.assert_not_called()

    def test_ensembl_timeout_keeps_directly_matched_analysis_running(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with patch("workflow.Path.cwd", return_value=Path(temporary)), patch("workflow.urlopen", side_effect=TimeoutError):
                self.assertEqual(_ensembl_for_symbols(["NOT_IN_CACHE"]), {})

    def test_edger_complete_requires_both_result_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "edger_results.csv").write_text("", encoding="utf-8")
            self.assertFalse(_edger_complete(output))
            (output / "logcpm.csv").write_text("", encoding="utf-8")
            self.assertTrue(_edger_complete(output))

    def test_candidate_ranking_uses_lowest_validation_auc_as_score(self) -> None:
        grouped = {("dev", "GENE1"): [
            {"cohort": "test1", "direction": 1.0, "auc": 0.9, "sensitivity": 0.8, "specificity": 0.7},
            {"cohort": "test2", "direction": 1.0, "auc": 0.8, "sensitivity": 0.9, "specificity": 0.8},
        ]}
        ranking = _candidate_ranking(grouped, 0.8, 0.7, 0.7)
        self.assertEqual((ranking[0]["rank"], ranking[0]["min_auc"], ranking[0]["passed"]), (1, 0.8, True))


if __name__ == "__main__":
    unittest.main()
