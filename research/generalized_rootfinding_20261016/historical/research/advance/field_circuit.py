"""Exact binary-field circuits with lazy affine XOR expressions.

This module changes the representation of a finite-field equation; it does
not weaken the equation or enumerate a field.  Nonlinear gates are ANDs of
signed Boolean literals and affine operations stay as XOR expressions until
an AND needs a literal.  A SAT/XOR backend can consume ``cnf_clauses`` and
``xor_rows`` directly.  The latter are ``(positive_variable_ids, parity)``.

There are no solver dependencies.  The original input coordinates retain
their one-based variable IDs, so projected models can be checked against the
independent polynomial exports without a coordinate conversion.
"""

from dataclasses import dataclass
from functools import lru_cache

from rootfinder.binaryfield import BinaryField


class CircuitBudgetError(RuntimeError):
    """An exact compilation exceeded an explicit engineering budget."""


@dataclass(frozen=True)
class AffineBit:
    variables: frozenset
    constant: bool = False

    @property
    def is_constant(self):
        return not self.variables


ZERO = AffineBit(frozenset(), False)
ONE = AffineBit(frozenset(), True)


@dataclass(frozen=True)
class Gate:
    kind: str
    output: int
    inputs: tuple


def _positive_int(value, name, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        raise ValueError(f"{name} must be positive")
    return value


class BinaryFieldCircuit:
    """Hash-consed Boolean circuit over a certified polynomial-basis field."""

    def __init__(self, field, n_primary, *, max_gates=1_000_000,
                 max_variables=1_000_000, max_xor_literals=50_000_000,
                 affine_cache_limit=65_536):
        if not isinstance(field, BinaryField):
            raise TypeError("field must be a certified BinaryField")
        self.field = field
        self.n_primary = _positive_int(n_primary, "n_primary", True)
        self.n_variables = self.n_primary
        self.max_gates = _positive_int(max_gates, "max_gates", True)
        self.max_variables = _positive_int(max_variables, "max_variables", True)
        self.max_xor_literals = _positive_int(max_xor_literals, "max_xor_literals", True)
        self.affine_cache_limit = _positive_int(affine_cache_limit, "affine_cache_limit", True)
        if self.n_primary > self.max_variables:
            raise CircuitBudgetError("primary variables exceed max_variables")
        self.cnf_clauses = []
        self.xor_rows = []
        self.gates = []
        self._affine_cache = {(frozenset(), False): ZERO, (frozenset(), True): ONE}
        self._literal_cache = {}
        self._and_cache = {}
        self._xor_literals = 0
        self._and_gate_count = 0
        self._xor_gate_count = 0
        self._constant_and_simplifications = 0
        self._and_cache_hits = 0
        self._literal_cache_hits = 0

    def affine(self, variables=(), constant=False):
        """Canonicalize an affine bit; repeated IDs cancel in characteristic two."""
        if not isinstance(constant, (bool, int)) or constant not in (0, 1):
            raise ValueError("constant must be a bit")
        if isinstance(variables, frozenset):
            support = variables
            for variable in support:
                self._check_variable(variable)
        else:
            support = set()
            for variable in variables:
                self._check_variable(variable)
                if variable in support:
                    support.remove(variable)
                else:
                    support.add(variable)
            support = frozenset(support)
        key = (support, bool(constant))
        cached = self._affine_cache.get(key)
        if cached is not None:
            return cached
        result = AffineBit(*key)
        # Dense expressions remain canonical by value, without retaining all
        # intermediate large sets in the interning cache.
        if len(support) <= 64 and len(self._affine_cache) < self.affine_cache_limit:
            self._affine_cache[key] = result
        return result

    def _check_variable(self, variable):
        _positive_int(variable, "variable ID")
        if variable > self.n_variables:
            raise ValueError("variable ID has not been allocated")

    def primary(self, variable):
        self._check_variable(variable)
        if variable > self.n_primary:
            raise ValueError("not a primary variable")
        return self.affine(frozenset((variable,)))

    def xor(self, a, b):
        if not isinstance(a, AffineBit) or not isinstance(b, AffineBit):
            raise TypeError("XOR inputs must be AffineBit values")
        if a == ZERO:
            return b
        if b == ZERO:
            return a
        return self.affine(a.variables ^ b.variables, a.constant ^ b.constant)

    def xor_many(self, bits):
        support = set()
        constant = False
        for bit in bits:
            if not isinstance(bit, AffineBit):
                raise TypeError("XOR inputs must be AffineBit values")
            support.symmetric_difference_update(bit.variables)
            constant ^= bit.constant
        return self.affine(frozenset(support), constant)

    def complement(self, bit):
        return self.affine(bit.variables, not bit.constant)

    def _new_variable(self):
        if len(self.gates) >= self.max_gates:
            raise CircuitBudgetError("generated gates exceed max_gates")
        if self.n_variables >= self.max_variables:
            raise CircuitBudgetError("variables exceed max_variables")
        self.n_variables += 1
        return self.n_variables

    def _append_xor(self, variables, parity):
        variables = tuple(sorted(variables))
        if self._xor_literals + len(variables) > self.max_xor_literals:
            raise CircuitBudgetError("XOR row literals exceed max_xor_literals")
        self.xor_rows.append((variables, bool(parity)))
        self._xor_literals += len(variables)

    def materialize(self, bit):
        """Return a signed literal; constant bits are deliberately not literals."""
        if bit.is_constant:
            raise ValueError("a constant has no signed variable literal")
        if len(bit.variables) == 1:
            result = next(iter(bit.variables))
        else:
            result = self._literal_cache.get(bit.variables)
            if result is None:
                # A wire for the uncomplemented support is shared by both
                # affine signs.  Definition: out XOR support = 0.
                if self._xor_literals + len(bit.variables) + 1 > self.max_xor_literals:
                    raise CircuitBudgetError("XOR row literals exceed max_xor_literals")
                result = self._new_variable()
                inputs = tuple(sorted(bit.variables))
                self._append_xor((*inputs, result), False)
                self.gates.append(Gate("xor", result, inputs))
                self._xor_gate_count += 1
                self._literal_cache[bit.variables] = result
            else:
                self._literal_cache_hits += 1
        return -result if bit.constant else result

    def and_bit(self, a, b):
        if not isinstance(a, AffineBit) or not isinstance(b, AffineBit):
            raise TypeError("AND inputs must be AffineBit values")
        if a.is_constant or b.is_constant:
            self._constant_and_simplifications += 1
            if a == ZERO or b == ZERO:
                return ZERO
            return b if a == ONE else a
        if a == b:
            return a
        if a.variables == b.variables:
            return ZERO  # a AND NOT a
        left, right = sorted((self.materialize(a), self.materialize(b)))
        key = (left, right)
        output = self._and_cache.get(key)
        if output is not None:
            self._and_cache_hits += 1
            return self.affine(frozenset((output,)))
        output = self._new_variable()
        self.cnf_clauses.extend(((-left, -right, output), (left, -output), (right, -output)))
        self.gates.append(Gate("and", output, key))
        self._and_gate_count += 1
        self._and_cache[key] = output
        return self.affine(frozenset((output,)))

    def element(self, value):
        value = self.field.validate(value)
        return tuple(ONE if value >> bit & 1 else ZERO for bit in range(self.field.degree))

    def input_element(self, offset, width=None, basis=None):
        """Embed consecutive primary bits; ``offset`` is zero based."""
        _positive_int(offset, "offset", True)
        if width is None:
            width = self.field.degree
        _positive_int(width, "width", True)
        if offset + width > self.n_primary:
            raise ValueError("input element exceeds the primary variable range")
        if basis is None:
            if width > self.field.degree:
                raise ValueError("polynomial prefix exceeds field degree")
            return tuple(self.primary(offset + bit + 1) if bit < width else ZERO
                         for bit in range(self.field.degree))
        if len(basis) != width:
            raise ValueError("basis length does not match input width")
        columns = [self.field.validate(value) for value in basis]
        return tuple(self.affine(frozenset(offset + index + 1 for index, value in enumerate(columns)
                                          if value >> bit & 1))
                     for bit in range(self.field.degree))

    def _check_element(self, a):
        if len(a) != self.field.degree or any(not isinstance(bit, AffineBit) for bit in a):
            raise ValueError("element must have exactly degree affine bits")

    def add(self, a, b):
        self._check_element(a)
        self._check_element(b)
        return tuple(self.xor(x, y) for x, y in zip(a, b))

    def _linear_map(self, a, columns):
        self._check_element(a)
        output = [[] for _ in range(self.field.degree)]
        for bit, column in zip(a, columns):
            if bit == ZERO:
                continue
            while column:
                low = column & -column
                output[low.bit_length() - 1].append(bit)
                column ^= low
        return tuple(self.xor_many(bits) for bits in output)

    def square(self, a):
        return self._linear_map(a, _square_columns(self.field.degree, self.field.modulus))

    def constant_mul(self, a, value):
        value = self.field.validate(value)
        return self._linear_map(a, _constant_columns(self.field.degree, self.field.modulus, value))

    @staticmethod
    def _trim_polynomial(a):
        length = len(a)
        while length and a[length - 1] == ZERO:
            length -= 1
        return a[:length]

    def _poly_add(self, a, b):
        return self._trim_polynomial(tuple(self.xor(a[i] if i < len(a) else ZERO,
                                                   b[i] if i < len(b) else ZERO)
                                           for i in range(max(len(a), len(b)))))

    def _poly_mul(self, a, b):
        """Recursive Karatsuba multiplication in GF(2)[z], pruning zero tails."""
        a = self._trim_polynomial(a)
        b = self._trim_polynomial(b)
        if not a or not b:
            return ()
        width = max(len(a), len(b))
        if width == 1:
            return (self.and_bit(a[0], b[0]),)
        split = (width + 1) // 2
        a0, a1 = a[:split], a[split:]
        b0, b1 = b[:split], b[split:]
        low = self._poly_mul(a0, b0)
        high = self._poly_mul(a1, b1)
        middle = self._poly_mul(self._poly_add(a0, a1), self._poly_add(b0, b1))
        middle = self._poly_add(self._poly_add(middle, low), high)
        # The three intermediate polynomials can extend beyond the final
        # degree for unequal input lengths.  Combine them before truncating:
        # those coefficients vanish algebraically, though an opaque AND
        # wire representation need not recognize distributive cancellation.
        output = [ZERO] * max(len(low), split + len(middle), 2 * split + len(high))
        for terms, offset in ((low, 0), (middle, split), (high, 2 * split)):
            for index, bit in enumerate(terms):
                position = offset + index
                output[position] = self.xor(output[position], bit)
        return self._trim_polynomial(tuple(output[:len(a) + len(b) - 1]))

    def mul(self, a, b):
        self._check_element(a)
        self._check_element(b)
        if a == b:
            return self.square(a)
        if all(bit.is_constant for bit in a):
            constant = sum(int(bit.constant) << index for index, bit in enumerate(a))
            return self.constant_mul(b, constant)
        if all(bit.is_constant for bit in b):
            constant = sum(int(bit.constant) << index for index, bit in enumerate(b))
            return self.constant_mul(a, constant)
        product = list(self._poly_mul(a, b))
        n = self.field.degree
        # Exact symbolic reduction by the field's monic defining polynomial.
        low_terms = tuple(bit for bit in range(n) if self.field.modulus >> bit & 1)
        for position in range(len(product) - 1, n - 1, -1):
            coefficient = product[position]
            if coefficient == ZERO:
                continue
            shift = position - n
            for bit in low_terms:
                product[shift + bit] = self.xor(product[shift + bit], coefficient)
        return tuple(product[bit] if bit < len(product) else ZERO for bit in range(n))

    def assert_zero(self, bit):
        if bit.is_constant:
            if bit.constant:
                self.cnf_clauses.append(())
            return
        self._append_xor(bit.variables, bit.constant)

    def evaluate_wires(self, assignment):
        """Evaluate every definition in creation order from primary bits."""
        if isinstance(assignment, int) and not isinstance(assignment, bool):
            if not 0 <= assignment < 1 << self.n_primary:
                raise ValueError("assignment integer is outside primary bit range")
            values = [False] + [bool(assignment >> bit & 1) for bit in range(self.n_primary)]
        else:
            if len(assignment) != self.n_primary or any(bit not in (0, 1, False, True) for bit in assignment):
                raise ValueError("assignment must contain exactly n_primary Boolean bits")
            values = [False] + [bool(bit) for bit in assignment]
        values.extend([False] * (self.n_variables - self.n_primary))
        for gate in self.gates:
            if gate.kind == "xor":
                value = False
                for variable in gate.inputs:
                    value ^= values[variable]
            else:
                value = all(values[abs(literal)] ^ (literal < 0) for literal in gate.inputs)
            values[gate.output] = value
        return tuple(values)

    @staticmethod
    def evaluate_bit(bit, wires):
        result = bit.constant
        for variable in bit.variables:
            result ^= wires[variable]
        return result

    def evaluate_element(self, a, wires):
        self._check_element(a)
        return sum(int(self.evaluate_bit(bit, wires)) << index for index, bit in enumerate(a))

    def statistics(self):
        return {
            "primary_variables": self.n_primary, "variables": self.n_variables,
            "generated_gates": len(self.gates), "and_gates": self._and_gate_count,
            "xor_gates": self._xor_gate_count, "cnf_clauses": len(self.cnf_clauses),
            "xor_rows": len(self.xor_rows), "xor_literals": self._xor_literals,
            "largest_xor_row": max((len(row) for row, _ in self.xor_rows), default=0),
            "and_cache_hits": self._and_cache_hits,
            "literal_cache_hits": self._literal_cache_hits,
            "constant_and_simplifications": self._constant_and_simplifications,
        }


@lru_cache(maxsize=32)
def _square_columns(degree, modulus):
    field = BinaryField(degree, modulus)
    return tuple(field.square(1 << bit) for bit in range(degree))


@lru_cache(maxsize=256)
def _constant_columns(degree, modulus, constant):
    field = BinaryField(degree, modulus)
    return tuple(field.mul(1 << bit, constant) for bit in range(degree))


@dataclass
class CompiledChain:
    circuit: BinaryFieldCircuit
    residual_bits: tuple
    metadata: dict
    definition_cnf_count: int
    definition_xor_count: int

    @property
    def primary_count(self):
        return self.circuit.n_primary

    def evaluate(self, assignment):
        wires = self.circuit.evaluate_wires(assignment)
        return sum(int(self.circuit.evaluate_bit(bit, wires)) << index
                   for index, bit in enumerate(self.residual_bits))

    def statistics(self):
        result = self.circuit.statistics()
        result.update({"extension_degree": self.circuit.field.degree,
                       "constraint_coordinates": len(self.residual_bits),
                       "definition_cnf_clauses": self.definition_cnf_count,
                       "definition_xor_rows": self.definition_xor_count})
        return result


def _parse_integer(value, name):
    if isinstance(value, str):
        return int(value, 0)
    return _positive_int(value, name, True)


def compile_s3_chain(instance, *, assert_constraints=True, **budgets):
    """Compile an exported binary-curve summation chain without expanded cubics.

    For each link, the exact equation is
    ``((a+b)*u + a*b)^2 + u*(a*b) + B = 0``.
    The final fixed target uses only one nonlinear product ``a*b``.  Source
    coordinates precede unrestricted auxiliary field coordinates exactly as
    in the exporter's layout.  All points are algebraic coordinate solutions;
    curve lifting and sign consistency are separate downstream checks.
    """
    extension = instance.get("extension_field", {})
    if extension.get("characteristic", 2) != 2:
        raise ValueError("S3 circuit compiler requires characteristic two")
    degree = _positive_int(extension.get("degree"), "extension degree")
    modulus = _parse_integer(extension.get("modulus_hex", extension.get("modulus")), "modulus")
    field = BinaryField(degree, modulus)
    dimension = _positive_int(instance.get("source_dimension"), "source dimension")
    if dimension > degree:
        raise ValueError("source dimension exceeds the extension degree")
    summands = _positive_int(instance.get("summands", 5), "summands")
    if summands < 2:
        raise ValueError("at least two summands are required")
    primary_count = summands * dimension + (summands - 2) * degree
    if instance.get("n_vars", primary_count) != primary_count:
        raise ValueError("primary variable count differs from source/auxiliary layout")
    expected_layout = ([{"index": index, "kind": "source_x", "offset": index * dimension,
                         "width": dimension} for index in range(summands)] +
                       [{"index": index, "kind": "auxiliary_x", "offset": summands * dimension + index * degree,
                         "width": degree} for index in range(summands - 2)])
    if "variable_layout" in instance and instance["variable_layout"] != expected_layout:
        raise ValueError("unsupported primary variable layout")
    expected_links = []
    if summands == 2:
        expected_links.append(["source_x:0", "source_x:1", "target_x"])
    else:
        expected_links.append(["source_x:0", "source_x:1", "auxiliary_x:0"])
        for index in range(1, summands - 2):
            expected_links.append([f"auxiliary_x:{index - 1}", f"source_x:{index + 1}", f"auxiliary_x:{index}"])
        expected_links.append([f"auxiliary_x:{summands - 3}", f"source_x:{summands - 1}", "target_x"])
    if "chain_links" in instance and instance["chain_links"] != expected_links:
        raise ValueError("unsupported S3 chain link order")
    basis = tuple(_parse_integer(value, "source basis") for value in
                  instance.get("source_basis_hex", [1 << bit for bit in range(dimension)]))
    if len(basis) != dimension:
        raise ValueError("source basis length differs from source dimension")
    target = field.validate(_parse_integer(instance.get("target", {}).get("x_hex"), "target x"))
    curve_b = field.validate(_parse_integer(instance.get("curve", {}).get("b_hex", "0x1"), "curve B"))
    circuit = BinaryFieldCircuit(field, primary_count, **budgets)
    source = [circuit.input_element(index * dimension, dimension, basis) for index in range(summands)]
    auxiliary = [circuit.input_element(summands * dimension + index * degree) for index in range(summands - 2)]
    bits = []
    for index in range(summands - 1):
        a = source[0] if index == 0 else auxiliary[index - 1]
        b = source[index + 1]
        product = circuit.mul(a, b)
        if index == summands - 2:
            # ((a+b)*T+ab)^2 + T*ab + B.  Frobenius and
            # fixed-target multiplication introduce no nonlinear gates.
            squared_sum = circuit.square(circuit.add(a, b))
            first = circuit.constant_mul(squared_sum, field.square(target))
            residual = circuit.add(circuit.add(first, circuit.square(product)),
                                   circuit.constant_mul(product, target))
        else:
            u = auxiliary[index]
            inner = circuit.add(circuit.mul(circuit.add(a, b), u), product)
            residual = circuit.add(circuit.square(inner), circuit.mul(u, product))
        residual = circuit.add(residual, circuit.element(curve_b))
        bits.extend(residual)
    definition_cnf_count = len(circuit.cnf_clauses)
    definition_xor_count = len(circuit.xor_rows)
    if assert_constraints:
        for bit in bits:
            circuit.assert_zero(bit)
    return CompiledChain(circuit, tuple(bits), {
        "instance_id": instance.get("instance_id"), "summands": summands,
        "source_dimension": dimension, "target_x": target, "curve_b": curve_b,
        "field_modulus": modulus, "constraints_asserted": bool(assert_constraints),
        "encoding": "lazy-affine-xor-karatsuba-and-v1",
    }, definition_cnf_count, definition_xor_count)
