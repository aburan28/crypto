"""Controls for exact residual-column and counting claims."""

import tempfile
import unittest
from pathlib import Path

import closure


def rank(rows: list[list[int]], prime: int) -> int:
    """Small independent dense modular echelon for the test only."""
    matrix = [row[:] for row in rows]
    pivot = 0
    for column in range(len(matrix[0])):
        selected = next((i for i in range(pivot, len(matrix))
                         if matrix[i][column] % prime), None)
        if selected is None:
            continue
        matrix[pivot], matrix[selected] = matrix[selected], matrix[pivot]
        inverse = pow(matrix[pivot][column], -1, prime)
        matrix[pivot] = [(value * inverse) % prime for value in matrix[pivot]]
        for i in range(pivot + 1, len(matrix)):
            factor = matrix[i][column]
            matrix[i] = [(a - factor * b) % prime
                         for a, b in zip(matrix[i], matrix[pivot])]
        pivot += 1
        if pivot == len(matrix):
            break
    return pivot


class ClosureTests(unittest.TestCase):
    def test_residual_nullity_guarantee_and_rank_limit(self):
        for rows in ([[1, 0], [0, 1], [1, 1], [1, 2]],
                     [[1, 1]] * 4,
                     [[0, 0]] * 4):
            self.assertGreaterEqual(len(rows) - rank(rows, 7), 2)
        # Residual-free combinations alone say nothing about independent
        # base columns: all four rows can have the same zero base projection.
        self.assertEqual(rank([[0, 0]] * 4, 7), 0)

    def test_selected_pair_and_saved_replay(self):
        result = closure.calculate()
        gray, hashed = result["selected_gray_hash_3_plus_2"]
        self.assertEqual((gray["residual_orbits"],
                          gray["residual_nullity_threshold_partials"],
                          gray["half_threshold_decomposition_query_floor"]),
                         (97, 697, 39468))
        self.assertEqual((hashed["residual_orbits"],
                          hashed["residual_nullity_threshold_partials"],
                          hashed["half_threshold_decomposition_query_floor"]),
                         (600, 1200, 1777))
        self.assertEqual(closure.OUTPUT.read_bytes(), closure.canonical(result))

    def test_modified_parent_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "screen.json"
            raw = closure.SUPPORT.read_bytes()
            altered.write_bytes(raw[:-2] + b" \n")
            with self.assertRaisesRegex(ValueError, "receipt changed"):
                closure.calculate(altered)


if __name__ == "__main__":
    unittest.main()
