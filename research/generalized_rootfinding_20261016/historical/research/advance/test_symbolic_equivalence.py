"""Exact tiny truth-table oracles and forged-export rejection checks."""

import copy
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest

from rootfinder.binaryfield import BinaryField
from research.advance.symbolic_equivalence import (
    BooleanPolynomialRing, SymbolicBudgetError, SymbolicEquivalenceError,
    coordinate_anf_sets, verify_file, verify_instance,
)


def direct_chain(data, assignment):
    n = data["extension_field"]["degree"]
    field = BinaryField(n, int(data["extension_field"]["modulus_hex"], 0))
    d, m = data["source_dimension"], data["summands"]
    basis = [int(value, 0) for value in data["source_basis_hex"]]
    sources = []
    for index in range(m):
        value = 0
        for bit, column in enumerate(basis):
            if assignment >> (index * d + bit) & 1:
                value ^= column
        sources.append(value)
    auxiliaries = [(assignment >> (m * d + index * n)) & ((1 << n) - 1)
                   for index in range(m - 2)]
    target = int(data["target"]["x_hex"], 0)
    b = int(data["curve"]["b_hex"], 0)
    output = 0
    for index in range(m - 1):
        x = sources[0] if index == 0 else auxiliaries[index - 1]
        y = sources[index + 1]
        z = target if index == m - 2 else auxiliaries[index]
        xy = field.mul(x, y)
        s3 = field.mul(field.square(x ^ y), field.square(z)) ^ field.mul(xy, z) ^ field.square(xy) ^ b
        output |= s3 << (index * n)
    return output


def refresh_declared_counts(data):
    data["equation_count"] = len(data["equations"])
    data["term_occurrences"] = sum(map(len, data["equations"]))
    data["max_degree"] = max((len(term) for row in data["equations"] for term in row), default=0)
    data["equations_sha256"] = hashlib.sha256(json.dumps(data["equations"], separators=(",", ":")).encode()).hexdigest()


def tiny_export(n=2, d=1, m=3, *, basis=None, target=3, curve_b=1):
    """Construct ANFs via exhaustive truth-table Mobius inversion, not symbolic algebra."""
    field = BinaryField(n)
    n_vars = m * d + (m - 2) * n
    if basis is None:
        basis = [1 << bit for bit in range(d)]
    if m == 2:
        links = [["source_x:0", "source_x:1", "target_x"]]
    else:
        links = [["source_x:0", "source_x:1", "auxiliary_x:0"]]
        links.extend([[f"auxiliary_x:{index - 1}", f"source_x:{index + 1}", f"auxiliary_x:{index}"]
                      for index in range(1, m - 2)])
        links.append([f"auxiliary_x:{m - 3}", f"source_x:{m - 1}", "target_x"])
    data = {
        "schema": "neutral_boolean_polynomial_system_v1", "instance_id": "tiny-truth-table-oracle",
        "field": {"kind": "prime", "p": 2}, "domain": [0, 1],
        "extension_field": {"characteristic": 2, "kind": "binary", "degree": n,
                            "modulus_hex": hex(field.modulus),
                            "low_terms": [bit for bit in range(n) if field.modulus >> bit & 1]},
        "source_dimension": d, "source_bits": m * d, "summands": m, "n_vars": n_vars,
        "source_basis_hex": [hex(value) for value in basis],
        "variable_layout": ([{"index": index, "kind": "source_x", "offset": index * d, "width": d}
                             for index in range(m)] +
                            [{"index": index, "kind": "auxiliary_x", "offset": m * d + index * n, "width": n}
                             for index in range(m - 2)]),
        "chain_links": links,
        "coordinate_equation_order": "S3 link index then polynomial-basis field coordinate, least significant bit first",
        "curve": {"model": "y^2+x*y=x^3+a*x^2+b", "a_hex": "0x0", "b_hex": hex(curve_b)},
        "target": {"x_hex": hex(target)},
    }
    truth_values = [direct_chain(data, assignment) for assignment in range(1 << n_vars)]
    coefficients = truth_values[:]
    for variable in range(n_vars):
        for mask in range(1 << n_vars):
            if mask >> variable & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << variable)]
    data["equations"] = [
        [[variable for variable in range(n_vars) if mask >> variable & 1]
         for mask, coefficient in enumerate(coefficients) if coefficient >> coordinate & 1]
        for coordinate in range((m - 1) * n)
    ]
    refresh_declared_counts(data)
    data["native_coordinate_evaluation_checks"] = [
        {"assignment_hex": hex(index % (1 << n_vars)),
         "coordinate_output_hex": hex(truth_values[index % (1 << n_vars)])}
        for index in range(8)
    ]
    return data


class SymbolicEquivalenceTests(unittest.TestCase):
    def test_coefficient_arithmetic_against_certified_field(self):
        rng = random.Random(324224)
        for n in (1, 2, 3, 4, 8, 11, 17, 83):
            field = BinaryField(n)
            ring = BooleanPolynomialRing(field)
            for _ in range(30):
                left, right = rng.randrange(field.order), rng.randrange(field.order)
                self.assertEqual(ring.coefficient_mul(left, right), field.mul(left, right))
                self.assertEqual(ring.coefficient_square(left), field.square(left))
            for i in range(n):
                j = n - 1 - i
                self.assertEqual(ring.coefficient_mul(1 << i, 1 << j), field.mul(1 << i, 1 << j))

    def test_boolean_multiplication_is_idempotent_and_cancels(self):
        ring = BooleanPolynomialRing(BinaryField(3))
        x, y = {1: 1}, {2: 1}
        self.assertEqual(ring.mul(x, x), x)
        self.assertEqual(ring.mul(x, y), {3: 1})
        self.assertEqual(ring.add(x, x), {})
        self.assertEqual(ring.mul(ring.add(x, y), ring.add(x, y)), ring.add(x, y))

    def test_frobenius_squares_coefficients_only(self):
        field = BinaryField(4)
        ring = BooleanPolynomialRing(field)
        polynomial = {0: 3, 1: 5, 3: 7, 4: 9}
        squared = ring.square(polynomial)
        self.assertEqual(squared, {mask: field.square(value) for mask, value in polynomial.items()})
        for assignment in range(8):
            self.assertEqual(ring.evaluate(squared, assignment), field.square(ring.evaluate(polynomial, assignment)))

    def test_tiny_full_truth_table_anf_equivalence(self):
        for n, d, m in ((1, 1, 2), (2, 1, 2), (2, 1, 3), (3, 1, 3), (2, 1, 5)):
            data = tiny_export(n, d, m, target=min(3, (1 << n) - 1))
            receipt = verify_instance(data)
            self.assertTrue(receipt["all_boolean_assignments_covered"])
            self.assertFalse(receipt["roots_solved"])
            self.assertEqual(receipt["term_occurrences"], data["term_occurrences"])

    def test_arbitrary_independent_basis_and_nontrivial_b(self):
        data = tiny_export(3, 2, 3, basis=[3, 5], target=6, curve_b=7)
        self.assertEqual(verify_instance(data)["equation_count"], 6)

    def test_forged_equation_with_valid_hash_and_matching_receipts_is_rejected(self):
        data = tiny_export(2, 1, 3)
        # All eight recorded assignments are zero, so toggling x_0 preserves
        # every receipt.  A valid recomputed digest still cannot prove ANF
        # equivalence, and the exact verifier must reject it.
        zero_output = direct_chain(data, 0)
        data["native_coordinate_evaluation_checks"] = [
            {"assignment_hex": "0x0", "coordinate_output_hex": hex(zero_output)} for _ in range(8)
        ]
        verify_instance(data)
        forged = copy.deepcopy(data)
        if [0] in forged["equations"][0]:
            forged["equations"][0].remove([0])
        else:
            forged["equations"][0].append([0])
        refresh_declared_counts(forged)
        with self.assertRaises(SymbolicEquivalenceError):
            verify_instance(forged)

    def test_wrong_semantic_metadata_with_valid_equation_hash_is_rejected(self):
        data = tiny_export(2, 1, 3)
        data["target"]["x_hex"] = "0x0"
        with self.assertRaises(SymbolicEquivalenceError):
            verify_instance(data)

    def test_frozen_file_hash_rejects_non_equation_metadata_change(self):
        data = tiny_export(2, 1, 3)
        raw = json.dumps(data).encode()
        digest = hashlib.sha256(raw).hexdigest()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.json"
            path.write_bytes(raw)
            self.assertTrue(verify_file(path, expected_file_sha256=digest)["frozen_file_hash_checked"])
            data["instance_id"] = "forged-label"
            path.write_text(json.dumps(data))
            with self.assertRaises(SymbolicEquivalenceError):
                verify_file(path, expected_file_sha256=digest)

    def test_repeated_or_noncanonical_monomials_are_rejected(self):
        data = tiny_export(2, 1, 3)
        for monomial in ([0, 0], [-1], [data["n_vars"]]):
            forged = copy.deepcopy(data)
            forged["equations"][0].append(monomial)
            refresh_declared_counts(forged)
            with self.assertRaises(SymbolicEquivalenceError):
                verify_instance(forged)

    def test_declared_counts_digest_layout_and_basis_are_checked(self):
        original = tiny_export(2, 1, 3)
        for key, value in (("term_occurrences", -1), ("equations_sha256", "0" * 64),
                           ("max_degree", 17), ("source_bits", 0), ("n_vars", 0),
                           ("chain_links", []), ("variable_layout", []),
                           ("source_basis_hex", ["0x0"])):
            forged = copy.deepcopy(original)
            forged[key] = value
            with self.assertRaises(SymbolicEquivalenceError):
                verify_instance(forged)

    def test_budget_refusals(self):
        ring = BooleanPolynomialRing(BinaryField(3), max_monomials=1)
        with self.assertRaises(SymbolicBudgetError):
            ring.input(0, [1, 2])
        with self.assertRaises(SymbolicBudgetError):
            coordinate_anf_sets({1: 7}, 3, max_coordinate_terms=2)
        with self.assertRaises(SymbolicBudgetError):
            verify_instance(tiny_export(2, 1, 3), max_monomials=1)


if __name__ == "__main__":
    unittest.main()
