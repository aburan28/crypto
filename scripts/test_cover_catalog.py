"""Metadata joins only. Cover mathematics is tested in the Rust binary."""
import copy
import unittest
import update_curve_standards as standards

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
        self.assertIn(browser.COVER_LINKS, data["sources"])
        for curve in data["curves"]:
            self.assertTrue(curve["cover_links"])
            self.assertTrue(curve["cover_links"][0]["cover_uid"].startswith("urn:hc-model:1:sha256:"))

    def test_source_inventory_and_idempotent_merge(self):
        catalog = browser.load(browser.REGISTRY)
        imported = browser.load("docs/curves/standards/registry.json")
        self.assertEqual(standards.merge(copy.deepcopy(catalog), imported), catalog)
        coverage = browser.load(browser.STANDARDS_COVERAGE)
        records = coverage["records"]
        source_rows = sum((browser.load("docs/curves/standards/" + name)["curves"]
                           for name in ("parameters.json", "supplemental.json")), [])
        self.assertEqual({r["source_id"] for r in records}, {r["source_id"] for r in source_rows})
        self.assertEqual(len(records), len(source_rows))
        slugs = {c["slug"] for c in catalog["curves"]}
        for row in records:
            if row["status"] == "imported":
                self.assertIn(row["icv1_slug"], slugs)
            else:
                self.assertIsNone(row["exists"])
                self.assertTrue(row["reason"])

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
