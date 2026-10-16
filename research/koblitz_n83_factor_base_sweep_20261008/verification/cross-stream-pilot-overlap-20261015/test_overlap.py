"""Corruption controls for the archived-orbit and ordered-pair auditors."""

import gzip
import json
from pathlib import Path
import tempfile
import unittest

import blake3

import audit
import overlap


class OrbitOverlapTests(unittest.TestCase):
    def make_archive(self, directory: Path, duplicate_phase: bool = False):
        header = {"curve_a": 0, "closure": "signed_frobenius", "orbit_columns": 1,
                  "point_count": 166, "policy": "public_x_hash", "seed": 1,
                  "frobenius_lambda": "2", "representatives": [[0, 1]]}
        lines = [json.dumps(header)]
        for phase in range(83):
            for negative in (False, True):
                used_phase = 0 if duplicate_phase and phase == 1 and negative else phase
                lines.append(json.dumps({"column": 0, "phase": used_phase,
                                         "negative": negative,
                                         "coefficient": str((-1 if negative else 1)
                                                            * pow(2, used_phase, overlap.ORDER)
                                                            % overlap.ORDER),
                                         "point": [2 * phase + int(negative), 1]}))
        plain = ("\n".join(lines) + "\n").encode()
        compressed = gzip.compress(plain, mtime=0)
        path = directory / "archive.jsonl.gz"
        path.write_bytes(compressed)
        expected = {"compressed_blake3": blake3.blake3(compressed).hexdigest(),
                    "plain_blake3": blake3.blake3(plain).hexdigest(),
                    "columns": 1, "policy": "public_x_hash", "seed": 1}
        return path, expected

    def test_two_independent_orbit_identifiers_and_duplicate_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path, expected = self.make_archive(Path(directory))
            hashes, first_sha = overlap.archive_orbits(path, expected)
            minima, second_sha = audit.minimum_point_orbits(path, expected)
            self.assertEqual((len(hashes), len(minima)), (1, 1))
            self.assertEqual(first_sha, second_sha)
            path, expected = self.make_archive(Path(directory), duplicate_phase=True)
            with self.assertRaises(ValueError):
                overlap.archive_orbits(path, expected)
            with self.assertRaises(ValueError):
                audit.minimum_point_orbits(path, expected)

    def test_auditor_rejects_changed_residual_count(self):
        populations = {(policy, seed, columns): set(range(columns))
                       for policy in ("a", "b", "c")
                       for seed in (1, 2, 3)
                       for columns in (64, 256, 600)}
        pairs = []
        for base, base_set in populations.items():
            for pool, pool_set in populations.items():
                if base[:2] == pool[:2]:
                    continue
                common = len(base_set & pool_set)
                pairs.append({"base": list(base), "pool": list(pool),
                              "shared_orbits": common,
                              "base_only_orbits": base[2] - common,
                              "residual_orbits": pool[2] - common,
                              "union_orbits": base[2] + pool[2] - common,
                              "residual_points": 166 * (pool[2] - common)})
        self.assertEqual(sum(audit.check_pairs(pairs, populations).values()), 648)
        pairs[0]["residual_orbits"] += 1
        with self.assertRaisesRegex(ValueError, "residual count changed"):
            audit.check_pairs(pairs, populations)


if __name__ == "__main__":
    unittest.main()
