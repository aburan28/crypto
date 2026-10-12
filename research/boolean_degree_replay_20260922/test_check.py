import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from check import InvalidEvidence, Ring, check_record, load_corpus


ROOT = Path(__file__).resolve().parent


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = load_corpus(ROOT)
        cls.sat = next(r for r in cls.rows if r["solution_assignments"])
        cls.unsat = next(r for r in cls.rows if not r["solution_assignments"])
        cls.growing = next(r for r in cls.rows
                           if r["truncated_closure_bound"] > r["input_degree"])

    def reject(self, row, update, message):
        changed = copy.deepcopy(row)
        update(changed)
        with self.assertRaisesRegex(InvalidEvidence, message):
            check_record(changed)

    def test_sat_and_unsat_fixtures(self):
        self.assertGreater(check_record(self.sat)["solutions"], 0)
        self.assertEqual(check_record(self.unsat)["solutions"], 0)

    def test_missing_solution(self):
        self.reject(self.sat, lambda r: r["solution_assignments"].pop(), "solution set")

    def test_invented_solution(self):
        self.reject(self.unsat, lambda r: r["solution_assignments"].append(0), "solution set")

    def test_wrong_basis(self):
        self.reject(self.sat, lambda r: r.update(boolean_gb_hex=["0x1"]), "basis fails")

    def test_missing_basis(self):
        self.reject(self.unsat, lambda r: r.update(boolean_gb_hex=[]), "dimension")

    def test_duplicate_basis_row(self):
        self.reject(self.sat, lambda r: r["boolean_gb_hex"].append(r["boolean_gb_hex"][0]),
                    "duplicate basis")

    def test_false_completion(self):
        self.reject(self.growing,
                    lambda r: r["profiles_through_bound"][0].update(contains_basis=True),
                    "profile mismatch")

    def test_wrong_rank(self):
        self.reject(self.sat, lambda r: r["profiles_through_bound"][0].update(rank=0),
                    "profile mismatch")

    def test_missing_lower_degree(self):
        self.reject(self.growing, lambda r: r["profiles_through_bound"].pop(0),
                    "missing degree")

    def test_wrong_input_degree(self):
        self.reject(self.sat, lambda r: r.update(input_degree=0), "input degree mismatch")

    def test_resource_and_type_limits(self):
        for variables in [0, 9, True, 1.0]:
            with self.subTest(variables=variables):
                self.reject(self.sat, lambda r: r.update(variables=variables), "variables")
        self.reject(self.sat, lambda r: r.update(equations_hex=["0x" + "f" * 1000]),
                    "encoding")
        self.reject(self.sat, lambda r: r.update(solution_assignments=[True]),
                    "invalid solution")

    def test_one_variable_decoder(self):
        self.assertEqual(Ring(1).decode("0x3"), 3)
        with self.assertRaises(InvalidEvidence):
            Ring(1).decode("0x4")

    def test_corpus_hash_missing_and_duplicate_records(self):
        original = (ROOT / "corpus.jsonl").read_bytes()
        manifest = json.loads((ROOT / "manifest.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for kind in ["hash", "missing", "duplicate"]:
                rows = copy.deepcopy(self.rows)
                if kind == "missing":
                    rows.pop()
                if kind == "duplicate":
                    rows[1] = rows[0]
                raw = (original + b" ") if kind == "hash" else (
                    "".join(json.dumps(r) + "\n" for r in rows).encode())
                metadata = dict(manifest)
                if kind != "hash":
                    metadata["corpus_sha256"] = hashlib.sha256(raw).hexdigest()
                (root / "manifest.json").write_text(json.dumps(metadata))
                (root / "corpus.jsonl").write_bytes(raw)
                with self.subTest(kind=kind), self.assertRaises(InvalidEvidence):
                    load_corpus(root)


if __name__ == "__main__":
    unittest.main()
