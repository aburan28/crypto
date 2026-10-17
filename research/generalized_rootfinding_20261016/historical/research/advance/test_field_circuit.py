"""Independent arithmetic, gate semantics, and exported-coordinate checks."""

import json
from pathlib import Path
import random
import unittest

from rootfinder.binaryfield import BinaryField
from research.advance.field_circuit import (
    BinaryFieldCircuit, CircuitBudgetError, ONE, ZERO, compile_s3_chain,
)


def clause_holds(clause, wires):
    return any(wires[abs(literal)] ^ (literal < 0) for literal in clause)


def xor_holds(row, wires):
    variables, parity = row
    actual = False
    for variable in variables:
        actual ^= wires[variable]
    return actual == parity


def direct_s3(field, a, b, u, curve_b):
    product = field.mul(a, b)
    inner = field.mul(a ^ b, u) ^ product
    return field.square(inner) ^ field.mul(u, product) ^ curve_b


def direct_chain(instance, assignment):
    n = instance["extension_field"]["degree"]
    field = BinaryField(n, int(instance["extension_field"]["modulus_hex"], 0))
    d = instance["source_dimension"]
    m = instance["summands"]
    basis = [int(value, 0) for value in instance["source_basis_hex"]]
    sources = []
    for index in range(m):
        value = 0
        for bit, column in enumerate(basis):
            if assignment >> (index * d + bit) & 1:
                value ^= column
        sources.append(value)
    auxiliaries = [(assignment >> (m * d + index * n)) & ((1 << n) - 1)
                   for index in range(m - 2)]
    target = int(instance["target"]["x_hex"], 0)
    curve_b = int(instance["curve"]["b_hex"], 0)
    result = 0
    for index in range(m - 1):
        a = sources[0] if index == 0 else auxiliaries[index - 1]
        b = sources[index + 1]
        u = target if index == m - 2 else auxiliaries[index]
        result |= direct_s3(field, a, b, u, curve_b) << (index * n)
    return result


def instance_for(field, dimension=2, summands=5, target=3):
    return {
        "instance_id": "test-independent-s3",
        "extension_field": {"characteristic": 2, "degree": field.degree,
                            "modulus_hex": hex(field.modulus)},
        "source_dimension": dimension, "summands": summands,
        "n_vars": summands * dimension + (summands - 2) * field.degree,
        "source_basis_hex": [hex(1 << bit) for bit in range(dimension)],
        "target": {"x_hex": hex(target)}, "curve": {"b_hex": "0x1"},
    }


class FieldCircuitTests(unittest.TestCase):
    def assert_definitions(self, circuit, wires, cnf_count=None, xor_count=None):
        clauses = circuit.cnf_clauses if cnf_count is None else circuit.cnf_clauses[:cnf_count]
        rows = circuit.xor_rows if xor_count is None else circuit.xor_rows[:xor_count]
        self.assertTrue(all(clause_holds(clause, wires) for clause in clauses))
        self.assertTrue(all(xor_holds(row, wires) for row in rows))

    def test_random_karatsuba_operations_against_field(self):
        rng = random.Random(919193)
        for degree in (1, 2, 3, 4, 5, 8, 11):
            field = BinaryField(degree)
            circuit = BinaryFieldCircuit(field, 2 * degree)
            a = circuit.input_element(0)
            b = circuit.input_element(degree)
            product = circuit.mul(a, b)
            added = circuit.add(a, b)
            squared = circuit.square(a)
            constant = field.validate(rng.randrange(field.order))
            scaled = circuit.constant_mul(a, constant)
            for _ in range(60):
                av = rng.randrange(field.order)
                bv = rng.randrange(field.order)
                wires = circuit.evaluate_wires(av | bv << degree)
                self.assertEqual(circuit.evaluate_element(product, wires), field.mul(av, bv))
                self.assertEqual(circuit.evaluate_element(added, wires), av ^ bv)
                self.assertEqual(circuit.evaluate_element(squared, wires), field.square(av))
                self.assertEqual(circuit.evaluate_element(scaled, wires), field.mul(av, constant))
                self.assert_definitions(circuit, wires)

    def test_square_and_constant_multiplication_have_no_and_gates(self):
        field = BinaryField(11, 0x805)
        circuit = BinaryFieldCircuit(field, 11)
        a = circuit.input_element(0)
        square = circuit.square(a)
        self.assertEqual(circuit.mul(a, a), square)
        circuit.constant_mul(square, 0x427)
        circuit.mul(a, circuit.element(0x427))
        self.assertEqual(circuit.statistics()["and_gates"], 0)
        self.assertEqual(circuit.n_variables, 11)
        self.assertEqual(circuit.xor_rows, [])

    def test_affine_canonicalization_and_lazy_materialization(self):
        circuit = BinaryFieldCircuit(BinaryField(4), 4)
        a = circuit.affine([1, 2, 2, 3], 1)
        b = circuit.xor(circuit.primary(1), circuit.primary(3))
        self.assertEqual(a, circuit.complement(b))
        self.assertIs(circuit.affine([1, 3]), b)
        self.assertEqual(circuit.n_variables, 4)
        positive = circuit.materialize(b)
        negative = circuit.materialize(a)
        self.assertEqual(positive, -negative)
        self.assertEqual(circuit.statistics()["xor_gates"], 1)
        for assignment in range(16):
            wires = circuit.evaluate_wires(assignment)
            self.assertEqual(wires[positive], bool((assignment & 1) ^ ((assignment >> 2) & 1)))
            self.assert_definitions(circuit, wires)

    def test_signed_and_truth_table_and_hash_consing(self):
        circuit = BinaryFieldCircuit(BinaryField(2), 2)
        a = circuit.primary(1)
        b = circuit.complement(circuit.primary(2))
        result = circuit.and_bit(a, b)
        self.assertEqual(result, circuit.and_bit(b, a))
        self.assertEqual(circuit.statistics()["and_gates"], 1)
        for assignment in range(4):
            wires = circuit.evaluate_wires(assignment)
            expected = bool(assignment & 1) and not bool(assignment >> 1 & 1)
            self.assertEqual(circuit.evaluate_bit(result, wires), expected)
            self.assert_definitions(circuit, wires)
            # The three clauses admit exactly the expected output bit.
            for output in (False, True):
                changed = list(wires)
                changed[circuit.gates[0].output] = output
                self.assertEqual(all(clause_holds(c, changed) for c in circuit.cnf_clauses), output == expected)

    def test_and_simplifications_do_not_generate_gates(self):
        circuit = BinaryFieldCircuit(BinaryField(2), 2)
        a = circuit.xor(circuit.primary(1), circuit.primary(2))
        self.assertEqual(circuit.and_bit(a, ZERO), ZERO)
        self.assertEqual(circuit.and_bit(a, ONE), a)
        self.assertEqual(circuit.and_bit(a, a), a)
        self.assertEqual(circuit.and_bit(a, circuit.complement(a)), ZERO)
        self.assertEqual(circuit.n_variables, 2)
        self.assertEqual(circuit.gates, [])

    def test_prefix_product_has_no_high_coefficients(self):
        field = BinaryField(11, 0x805)
        circuit = BinaryFieldCircuit(field, 6)
        product = circuit.mul(circuit.input_element(0, 3), circuit.input_element(3, 3))
        self.assertEqual(product[5:], (ZERO,) * 6)
        self.assertLessEqual(circuit.statistics()["and_gates"], 9)

    def test_chain_matches_direct_field_equations(self):
        rng = random.Random(21071)
        for summands in (2, 3, 4, 5):
            instance = instance_for(BinaryField(4), summands=summands)
            compiled = compile_s3_chain(instance)
            for _ in range(40):
                assignment = rng.getrandbits(compiled.primary_count)
                expected = direct_chain(instance, assignment)
                self.assertEqual(compiled.evaluate(assignment), expected)
                wires = compiled.circuit.evaluate_wires(assignment)
                self.assert_definitions(compiled.circuit, wires, compiled.definition_cnf_count,
                                        compiled.definition_xor_count)
                all_constraints = (all(clause_holds(c, wires) for c in compiled.circuit.cnf_clauses) and
                                   all(xor_holds(r, wires) for r in compiled.circuit.xor_rows))
                self.assertEqual(all_constraints, expected == 0)

    def test_arbitrary_source_embedding_and_curve_b(self):
        instance = instance_for(BinaryField(5), dimension=3, target=7)
        instance["source_basis_hex"] = ["0x1", "0xa", "0x17"]
        instance["curve"]["b_hex"] = "0x19"
        compiled = compile_s3_chain(instance)
        rng = random.Random(141)
        for _ in range(40):
            assignment = rng.getrandbits(compiled.primary_count)
            self.assertEqual(compiled.evaluate(assignment), direct_chain(instance, assignment))

    def test_independent_export_receipts_and_expanded_polynomials(self):
        path = Path(__file__).parent / "instances" / "n11-d3-m5-target17.json"
        instance = json.loads(path.read_text())
        compiled = compile_s3_chain(instance)
        for receipt in instance["native_coordinate_evaluation_checks"]:
            assignment = int(receipt["assignment_hex"], 0)
            expected = int(receipt["coordinate_output_hex"], 0)
            self.assertEqual(compiled.evaluate(assignment), expected)
            self.assertEqual(direct_chain(instance, assignment), expected)
            output = 0
            for index, equation in enumerate(instance["equations"]):
                bit = 0
                for term in equation:
                    bit ^= int(all(assignment >> variable & 1 for variable in term))
                output |= bit << index
            self.assertEqual(output, expected)

    def test_constraints_can_be_left_unasserted(self):
        compiled = compile_s3_chain(instance_for(BinaryField(4)), assert_constraints=False)
        self.assertEqual(len(compiled.circuit.cnf_clauses), compiled.definition_cnf_count)
        self.assertEqual(len(compiled.circuit.xor_rows), compiled.definition_xor_count)

    def test_budget_refusals_are_explicit(self):
        field = BinaryField(4)
        with self.assertRaises(CircuitBudgetError):
            BinaryFieldCircuit(field, 4, max_variables=3)
        circuit = BinaryFieldCircuit(field, 4, max_gates=0)
        with self.assertRaises(CircuitBudgetError):
            circuit.and_bit(circuit.primary(1), circuit.primary(2))
        circuit = BinaryFieldCircuit(field, 4, max_xor_literals=2)
        with self.assertRaises(CircuitBudgetError):
            circuit.materialize(circuit.xor(circuit.primary(1), circuit.primary(2)))
        self.assertEqual(circuit.n_variables, 4)
        with self.assertRaises(CircuitBudgetError):
            compile_s3_chain(instance_for(field), max_gates=1)

    def test_invalid_metadata_is_rejected(self):
        instance = instance_for(BinaryField(4))
        instance["n_vars"] += 1
        with self.assertRaises(ValueError):
            compile_s3_chain(instance)
        instance = instance_for(BinaryField(4))
        instance["extension_field"]["characteristic"] = 3
        with self.assertRaises(ValueError):
            compile_s3_chain(instance)


if __name__ == "__main__":
    unittest.main()
