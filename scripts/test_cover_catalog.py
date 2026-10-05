"""Metadata joins only. Cover mathematics is tested in the Rust binary."""
import copy
import unittest

import build_lab_browser as browser


class CoverCatalogTests(unittest.TestCase):
    def test_complete_catalog_join(self):
        data = browser.build()
        self.assertEqual(data["counts"]["curves"], len(browser.load(browser.REGISTRY)["curves"]))
        for curve in data["curves"]:
            finding = curve["hyperelliptic_cover"]
            self.assertEqual(finding["slug"], curve["slug"])
            self.assertEqual(finding["status"], "verified_over_declared_field")
        self.assertIn(browser.COVERS, data["sources"])

    def test_stale_and_ambiguous_findings_rejected(self):
        curves = browser.curve_rows(browser.load(browser.REGISTRY), browser.load(browser.LEADERBOARD))
        original = browser.load(browser.COVERS)
        for mutation in ("registry", "model", "missing", "duplicate"):
            report = copy.deepcopy(original)
            if mutation == "registry":
                report["registry_sha256"] = "0" * 64
            elif mutation == "model":
                report["curves"][0]["model_sha256"] = "0" * 64
            elif mutation == "missing":
                report["curves"].pop()
            else:
                report["curves"].append(report["curves"][0])
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                browser.attach_covers(copy.deepcopy(curves), report)


if __name__ == "__main__":
    unittest.main()
