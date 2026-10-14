#!/usr/bin/env python3
"""Regression for consuming explicit model forms without changing their meaning."""
import json
import unittest
from unittest.mock import patch

import build_curve_registry as builder


class NativePrimeModelTests(unittest.TestCase):
    def test_quadratic_records_use_the_native_identity_and_short_coefficients(self):
        source = builder.REPO / (
            "research/isogeny_fibers_20261010/round4/"
            "16_7_agm_isogeny_oracles/out_seed7.json")
        registry = builder.Registry()
        problems, fatal = [], []
        with patch.object(builder, "tracked", return_value=[source]):
            builder.harvest_json(registry, problems, fatal)
        self.assertEqual(problems, [])
        self.assertEqual(fatal, [])
        rows = {entry["slug"]: entry for entry in registry.curves.values()}
        self.assertEqual(set(rows), {
            "icv1-fp21-tm700-d62ba3a6", "icv1-fp21-tm700-07d957b2"})
        for slug, a, b, j in [
                ("icv1-fp21-tm700-d62ba3a6", "699056", "737891", "128"),
                ("icv1-fp21-tm700-07d957b2", "699046", "388354", "10976")]:
            row = rows[slug]
            self.assertEqual(row["params"], {"p": "1048583", "a": a, "b": b})
            self.assertEqual(row["order"], "1049284")
            self.assertEqual(row["j"], j)
            self.assertEqual(row["sources"], [str(source.relative_to(builder.REPO))])

    def test_missing_stale_duplicate_or_wrong_pointer_native_records_are_rejected(self):
        source = builder.REPO / (
            "research/isogeny_fibers_20261010/round4/"
            "16_7_agm_isogeny_oracles/out_seed7.json")
        cache = builder.REPO / "docs/curves/sources/prime-models.json"
        original_read = builder.Path.read_text
        original_doc = json.loads(cache.read_text())
        for failure in ("missing", "stale", "duplicate", "pointer"):
            doc = json.loads(json.dumps(original_doc))
            if failure == "missing":
                doc["records"] = []
            elif failure == "stale":
                for row in doc["records"]:
                    row["source_sha256"] = "0" * 64
            elif failure == "duplicate":
                doc["records"] *= 2
            else:
                for row in doc["records"]:
                    row["pointer"] = "/wrong"
            def read(path, *args, **kwargs):
                return json.dumps(doc) if path == cache else original_read(path, *args, **kwargs)
            registry = builder.Registry()
            problems, fatal = [], []
            with self.subTest(failure=failure), patch.object(builder, "tracked", return_value=[source]), \
                    patch.object(builder.Path, "read_text", autospec=True, side_effect=read):
                builder.harvest_json(registry, problems, fatal)
            self.assertEqual(registry.curves, {})
            self.assertEqual(len(fatal), 2)
            self.assertTrue(all("normalization is missing or stale" in error for error in fatal))

    def test_an_unknown_explicit_form_is_rejected_before_short_model_extraction(self):
        source = builder.REPO / (
            "research/isogeny_fibers_20261010/round4/"
            "16_7_agm_isogeny_oracles/out_seed7.json")
        original_read = builder.Path.read_text
        document = json.loads(source.read_text())
        document["curve_E"]["model"] = "unsupported coordinate equation"
        document["curve_E2"]["model"] = "unsupported coordinate equation"
        def read(path, *args, **kwargs):
            return json.dumps(document) if path == source else original_read(path, *args, **kwargs)
        registry = builder.Registry()
        problems, fatal = [], []
        with patch.object(builder, "tracked", return_value=[source]), \
                patch.object(builder.Path, "read_text", autospec=True, side_effect=read):
            builder.harvest_json(registry, problems, fatal)
        self.assertEqual(registry.curves, {})
        self.assertEqual(len(fatal), 2)
        self.assertTrue(all("unsupported explicit prime model" in error for error in fatal))


if __name__ == "__main__":
    unittest.main()
