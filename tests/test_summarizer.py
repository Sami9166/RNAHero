import csv
import json
import tempfile
import unittest
from pathlib import Path

from summarizer import _top_five


class SummarizerTests(unittest.TestCase):
    def test_top_five_uses_internal_fdr_significant_logfc_from_one_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            analysis = root / "analysis"
            for cohort, rows in {
                "dev1": [["A", "2.0", "0.01"], ["B", "-4.0", "0.01"]],
                "dev2": [["C", "10.0", "0.2"]],
            }.items():
                path = analysis / cohort
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
            (analysis / "validation_report.json").write_text(json.dumps(report), encoding="utf-8")
            cohort, genes, _ = _top_five(root)
            self.assertEqual(cohort, "dev1")
            self.assertEqual([gene["gene_id"] for gene in genes], ["B", "A"])


if __name__ == "__main__":
    unittest.main()
