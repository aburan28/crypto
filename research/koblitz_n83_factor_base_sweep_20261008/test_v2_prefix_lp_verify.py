"""Controls for exact prefix acceptance and single-row corruption."""

import gzip
import json
from pathlib import Path
import tempfile
import unittest

import v2_prefix_lp_verify as verifier


def object_file(path: Path, columns: int, altered: bool = False) -> None:
    header = {
        "study": verifier.STUDY,
        "schema": "n83.factor-base/v1",
        "curve_a": 0,
        "curve_slug": "fixture-curve",
        "n": 7,
        "policy": "public_x_hash",
        "seed": 1,
        "closure": "signed_frobenius",
        "modulus_low_terms": [0, 1],
        "subgroup_order": "101",
        "cofactor": "2",
        "frobenius_lambda": "3",
        "generator": ["1", "2"],
        "orbit_columns": columns,
        "point_count": columns * 2,
        "representatives": [[column, 0] for column in range(columns)],
    }
    with gzip.open(path, "wb") as output:
        output.write(json.dumps(header).encode() + b"\n")
        for column in range(columns):
            for phase in range(2):
                point = [column, phase]
                if altered and column == 0 and phase == 1:
                    point = [column, 99]
                output.write(json.dumps({"column": column, "point": point}).encode() + b"\n")


class PrefixVerifierTests(unittest.TestCase):
    def test_exact_prefix_and_disjoint_remainder(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            small = Path(directory) / "small.gz"
            large = Path(directory) / "large.gz"
            object_file(small, 2)
            object_file(large, 3)
            result = verifier.check_prefix(small, large, 2, 3, 2)
            self.assertEqual(result["small_point_rows"], 4)
            self.assertEqual(result["residual_columns"], 1)
            self.assertEqual(result["residual_point_rows"], 2)
            self.assertTrue(result["point_row_prefix"])

    def test_changed_shared_point_row_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            small = Path(directory) / "small.gz"
            large = Path(directory) / "large.gz"
            object_file(small, 2)
            object_file(large, 3, altered=True)
            with self.assertRaisesRegex(ValueError, "point-row prefix differs"):
                verifier.check_prefix(small, large, 2, 3, 2)

    def test_changed_representative_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            small = Path(directory) / "small.gz"
            large = Path(directory) / "large.gz"
            object_file(small, 2)
            object_file(large, 3)
            with gzip.open(large, "rb") as source:
                lines = source.readlines()
            header = json.loads(lines[0])
            header["representatives"][0] = [77, 0]
            with gzip.open(large, "wb") as target:
                target.write(json.dumps(header).encode() + b"\n")
                target.writelines(lines[1:])
            with self.assertRaisesRegex(ValueError, "representative list"):
                verifier.check_prefix(small, large, 2, 3, 2)


if __name__ == "__main__":
    unittest.main()
