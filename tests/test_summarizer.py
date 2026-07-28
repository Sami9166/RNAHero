import csv
import json
import tempfile
import unittest
from pathlib import Path

from summarizer import _go_gene_sets, _top_five


class SummarizerTests(unittest.TestCase):
    def test_go_gene_sets_use_fdr_and_direction(self) -> None:
        rows = [
            {"gene_id": "UP", "logFC": "1.2", "FDR": "0.01"},
            {"gene_id": "DOWN", "logFC": "-1.2", "FDR": "0.01"},
            {"gene_id": "NOT_FDR", "logFC": "2", "FDR": "0.2"},
        ]
        self.assertEqual(_go_gene_sets(rows, 0.05), {"up": ["UP"], "down": ["DOWN"]})

    def test_top_five_uses_internal_fdr_significant_logfc_from_one_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            analysis = root / "analysis"
            for cohort, rows in {
                "dev1": [["A", "2.0", "0.01"], ["B", "-4.0", "0.01"]],
                "dev2": [["C", "10.0", "0.2"]],
            }.items():
                path = analysis / "cohorts" / cohort / "edgeR"
                path.mkdir(parents=True)
                with (path / "edger_results.csv").open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(["gene_id", "logFC", "FDR"])
                    writer.writerows(rows)
            report = {"fdr": 0.05, "candidate_ranking": [
                {"gene_id": "A", "discovery": "dev1", "passed": True},
                {"gene_id": "B", "discovery": "dev1", "passed": True},
                {"gene_id": "C", "discovery": "dev2", "passed": True},
            ]}
            (analysis / "validation").mkdir(parents=True)
            (analysis / "validation" / "validation_report.json").write_text(json.dumps(report), encoding="utf-8")
            cohort, genes, _ = _top_five(root)
            self.assertEqual(cohort, "dev1")
            self.assertEqual([gene["gene_id"] for gene in genes], ["B", "A"])


if __name__ == "__main__":
    unittest.main()
