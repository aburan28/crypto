"""Prove exported Boolean ANFs equal the stated binary-field S3 chain.

This is an independent input-equivalence check, not a root solver or timing
measurement.  It does not import the field-circuit compiler.  A field-valued
Boolean polynomial is a sparse ``variable_mask -> field_coefficient`` map.
Multiplication ORs variable masks; equal monomials combine by XOR.  Squaring
squares coefficients only because Boolean monomials are idempotent.

All coordinate monomial sets are compared exactly.  Random evaluation
receipts are supplementary checks and cannot substitute for this proof.
For preregistered input integrity, pass its frozen SHA-256 with
``verify_file(..., expected_file_sha256=...)``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rootfinder.binaryfield import BinaryField


class SymbolicEquivalenceError(ValueError):
    """An export or its declared S3 semantics is inconsistent."""


class SymbolicBudgetError(RuntimeError):
    """An explicit symbolic-verification resource cap was reached."""


def _integer(value, label, *, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise SymbolicEquivalenceError(f"{label} must be an integer >= {minimum}")
    return value


def _encoded_integer(value, label):
    if isinstance(value, str):
        try:
            result = int(value, 0)
        except ValueError as error:
            raise SymbolicEquivalenceError(f"invalid {label}") from error
        return _integer(result, label)
    return _integer(value, label)


class BooleanPolynomialRing:
    """Exact sparse polynomials in GF(2^n)[x]/(x_i^2-x_i)."""

    def __init__(self, field, *, max_monomials=1_000_000,
                 coefficient_cache_limit=500_000):
        if not isinstance(field, BinaryField):
            raise TypeError("field must be a certified BinaryField")
        self.field = field
        self.max_monomials = _integer(max_monomials, "max_monomials", minimum=1)
        self.coefficient_cache_limit = _integer(coefficient_cache_limit, "coefficient_cache_limit")
        self._products = {}
        self._squares = {}
        self.coefficient_products = 0
        self.coefficient_cache_hits = 0
        self.largest_polynomial = 0
        # Basis powers cover every unreduced product of three field basis
        # columns.  Exact polynomial remainders justify this optimization.
        self._basis_powers = [1]
        for _ in range(3 * (field.degree - 1)):
            self._basis_powers.append(self._reduce(self._basis_powers[-1] << 1))

    def _reduce(self, value):
        degree = self.field.degree
        while value.bit_length() > degree:
            value ^= self.field.modulus << (value.bit_length() - degree - 1)
        return value

    def coefficient_mul(self, left, right):
        """Independent carryless product/remainder, cached by exact operands."""
        if not left or not right:
            return 0
        if left == 1:
            return right
        if right == 1:
            return left
        if left & (left - 1) == 0 and right & (right - 1) == 0:
            exponent = left.bit_length() + right.bit_length() - 2
            return self._basis_powers[exponent]
        if left > right:
            left, right = right, left
        key = (left, right)
        cached = self._products.get(key)
        if cached is not None:
            self.coefficient_cache_hits += 1
            return cached
        self.coefficient_products += 1
        if left & (left - 1) == 0:
            result = self._reduce(right << (left.bit_length() - 1))
        elif right & (right - 1) == 0:
            result = self._reduce(left << (right.bit_length() - 1))
        else:
            # Iterating set bits rather than every bit keeps sparse defining
            # polynomials and basis-derived coefficients inexpensive.
            if left.bit_count() > right.bit_count():
                left, right = right, left
            product = 0
            remaining = left
            while remaining:
                low = remaining & -remaining
                product ^= right << (low.bit_length() - 1)
                remaining ^= low
            result = self._reduce(product)
        if len(self._products) < self.coefficient_cache_limit:
            self._products[key] = result
        return result

    def coefficient_square(self, value):
        cached = self._squares.get(value)
        if cached is not None:
            return cached
        product = 0
        remaining = value
        while remaining:
            low = remaining & -remaining
            product ^= 1 << (2 * (low.bit_length() - 1))
            remaining ^= low
        result = self._reduce(product)
        self._squares[value] = result
        return result

    def _budget(self, polynomial):
        size = len(polynomial)
        if size > self.max_monomials:
            raise SymbolicBudgetError("field-valued polynomial exceeds max_monomials")
        self.largest_polynomial = max(self.largest_polynomial, size)
        return polynomial

    @staticmethod
    def _merge(polynomial, mask, coefficient):
        coefficient ^= polynomial.get(mask, 0)
        if coefficient:
            polynomial[mask] = coefficient
        else:
            polynomial.pop(mask, None)

    def constant(self, coefficient):
        coefficient = self.field.validate(coefficient)
        return {0: coefficient} if coefficient else {}

    def input(self, offset, basis):
        _integer(offset, "input offset")
        result = {1 << (offset + index): self.field.validate(coefficient)
                  for index, coefficient in enumerate(basis) if coefficient}
        return self._budget(result)

    def add(self, *polynomials):
        result = {}
        for polynomial in polynomials:
            for mask, coefficient in polynomial.items():
                self._merge(result, mask, coefficient)
        return self._budget(result)

    def square(self, polynomial):
        return self._budget({mask: self.coefficient_square(coefficient)
                             for mask, coefficient in polynomial.items()})

    def mul(self, left, right):
        result = {}
        for left_mask, left_coefficient in left.items():
            for right_mask, right_coefficient in right.items():
                coefficient = self.coefficient_mul(left_coefficient, right_coefficient)
                self._merge(result, left_mask | right_mask, coefficient)
            if len(result) > self.max_monomials:
                raise SymbolicBudgetError("field-valued polynomial exceeds max_monomials")
        return self._budget(result)

    def s3(self, a, b, c, curve_b):
        # This expansion is deliberately different from the compiler's
        # ((a+b)*c+ab)^2+c*ab+B circuit construction.
        product = self.mul(a, b)
        quadratic = self.mul(self.square(self.add(a, b)), self.square(c))
        cubic = self.mul(product, c)
        return self.add(quadratic, cubic, self.square(product), self.constant(curve_b))

    @staticmethod
    def evaluate(polynomial, assignment):
        value = 0
        for mask, coefficient in polynomial.items():
            if assignment & mask == mask:
                value ^= coefficient
        return value


def _metadata(instance):
    if instance.get("schema") != "neutral_boolean_polynomial_system_v1":
        raise SymbolicEquivalenceError("unsupported exported schema")
    if instance.get("field") != {"kind": "prime", "p": 2}:
        raise SymbolicEquivalenceError("exported coordinate field must be GF(2)")
    if instance.get("domain") != [0, 1] or any(type(bit) is not int for bit in instance["domain"]):
        raise SymbolicEquivalenceError("exported coordinate domain must be the Boolean bits")
    extension = instance.get("extension_field", {})
    if extension.get("characteristic") != 2 or extension.get("kind") != "binary":
        raise SymbolicEquivalenceError("extension metadata must specify a binary field")
    n = _integer(extension.get("degree"), "extension degree", minimum=1)
    modulus = _encoded_integer(extension.get("modulus_hex"), "field modulus")
    try:
        field = BinaryField(n, modulus)
    except (TypeError, ValueError) as error:
        raise SymbolicEquivalenceError(str(error)) from error
    if "low_terms" in extension:
        expected = [bit for bit in range(n) if modulus >> bit & 1]
        if extension["low_terms"] != expected:
            raise SymbolicEquivalenceError("low_terms disagree with the exact field modulus")
    d = _integer(instance.get("source_dimension"), "source dimension", minimum=1)
    m = _integer(instance.get("summands"), "summands", minimum=2)
    if d > n:
        raise SymbolicEquivalenceError("source dimension exceeds extension degree")
    count = m * d + (m - 2) * n
    if _integer(instance.get("n_vars"), "n_vars") != count:
        raise SymbolicEquivalenceError("n_vars disagrees with source/auxiliary dimensions")
    if "source_bits" in instance and instance["source_bits"] != m * d:
        raise SymbolicEquivalenceError("source_bits disagrees with source dimensions")
    values = instance.get("source_basis_hex")
    if not isinstance(values, list) or len(values) != d:
        raise SymbolicEquivalenceError("source basis length disagrees with source dimension")
    basis = tuple(field.validate(_encoded_integer(value, "source basis column")) for value in values)
    echelon = {}
    for column in basis:
        remaining = column
        while remaining:
            pivot = remaining.bit_length() - 1
            if pivot not in echelon:
                echelon[pivot] = remaining
                break
            remaining ^= echelon[pivot]
        if not remaining:
            raise SymbolicEquivalenceError("source basis columns are linearly dependent")
    layout = ([{"index": index, "kind": "source_x", "offset": index * d, "width": d}
               for index in range(m)] +
              [{"index": index, "kind": "auxiliary_x", "offset": m * d + index * n, "width": n}
               for index in range(m - 2)])
    if instance.get("variable_layout") != layout:
        raise SymbolicEquivalenceError("variable layout disagrees with source/auxiliary coordinates")
    if m == 2:
        links = [["source_x:0", "source_x:1", "target_x"]]
    else:
        links = [["source_x:0", "source_x:1", "auxiliary_x:0"]]
        links.extend([[f"auxiliary_x:{index - 1}", f"source_x:{index + 1}", f"auxiliary_x:{index}"]
                      for index in range(1, m - 2)])
        links.append([f"auxiliary_x:{m - 3}", f"source_x:{m - 1}", "target_x"])
    if instance.get("chain_links") != links:
        raise SymbolicEquivalenceError("chain link order disagrees with S3 semantics")
    if "coordinate_equation_order" in instance and instance["coordinate_equation_order"] != (
            "S3 link index then polynomial-basis field coordinate, least significant bit first"):
        raise SymbolicEquivalenceError("unsupported coordinate equation order")
    target = field.validate(_encoded_integer(instance.get("target", {}).get("x_hex"), "target x"))
    curve = instance.get("curve", {})
    if "model" in curve and curve["model"] != "y^2+x*y=x^3+a*x^2+b":
        raise SymbolicEquivalenceError("unsupported curve model")
    if "a_hex" in curve:
        field.validate(_encoded_integer(curve["a_hex"], "curve A"))
    b = field.validate(_encoded_integer(curve.get("b_hex"), "curve B"))
    if b == 0:
        raise SymbolicEquivalenceError("binary elliptic-curve B must be nonzero")
    return field, d, m, count, basis, target, b


def symbolic_links(instance, *, max_monomials=1_000_000, coefficient_cache_limit=500_000):
    """Yield exact field-valued link polynomials, sharing coefficient caches."""
    field, d, m, _, basis, target, b = _metadata(instance)
    ring = BooleanPolynomialRing(field, max_monomials=max_monomials,
                                 coefficient_cache_limit=coefficient_cache_limit)
    source = [ring.input(index * d, basis) for index in range(m)]
    auxiliary = [ring.input(m * d + index * field.degree,
                            [1 << bit for bit in range(field.degree)])
                 for index in range(m - 2)]
    for index in range(m - 1):
        a = source[0] if index == 0 else auxiliary[index - 1]
        c = ring.constant(target) if index == m - 2 else auxiliary[index]
        yield index, ring.s3(a, source[index + 1], c, b), ring


def coordinate_anf_sets(polynomial, degree, *, max_coordinate_terms=20_000_000):
    """Split field coefficients into canonical coordinate monomial sets."""
    _integer(max_coordinate_terms, "max_coordinate_terms", minimum=1)
    rows = [set() for _ in range(degree)]
    count = 0
    for mask, coefficient in polynomial.items():
        while coefficient:
            low = coefficient & -coefficient
            index = low.bit_length() - 1
            if index >= degree:
                raise SymbolicEquivalenceError("noncanonical field coefficient")
            rows[index].add(mask)
            count += 1
            if count > max_coordinate_terms:
                raise SymbolicBudgetError("coordinate ANF terms exceed max_coordinate_terms")
            coefficient ^= low
    return rows


def _canonical_row_digest(digest, row, byte_width):
    digest.update(len(row).to_bytes(8, "little"))
    for mask in sorted(row):
        digest.update(mask.to_bytes(byte_width, "little"))


def verify_instance(instance, *, max_monomials=1_000_000,
                    max_coordinate_terms=20_000_000, coefficient_cache_limit=500_000):
    """Compare every original coordinate polynomial to exact S3 semantics."""
    started = time.perf_counter()
    field, _, m, n_vars, _, _, _ = _metadata(instance)
    equations = instance.get("equations")
    if not isinstance(equations, list) or len(equations) != (m - 1) * field.degree:
        raise SymbolicEquivalenceError("equation list does not contain every ordered S3 coordinate")
    if instance.get("equation_count") != len(equations):
        raise SymbolicEquivalenceError("equation_count disagrees with equation list")
    equations_digest = hashlib.sha256(json.dumps(equations, separators=(",", ":")).encode()).hexdigest()
    if instance.get("equations_sha256") != equations_digest:
        raise SymbolicEquivalenceError("declared equation SHA-256 does not match equation list")
    receipts = instance.get("native_coordinate_evaluation_checks", [])
    if not isinstance(receipts, list):
        raise SymbolicEquivalenceError("coordinate receipts must be a list")
    assignments = []
    for receipt in receipts:
        assignment = _encoded_integer(receipt.get("assignment_hex"), "receipt assignment")
        expected = _encoded_integer(receipt.get("coordinate_output_hex"), "receipt output")
        if assignment >= 1 << n_vars or expected >= 1 << len(equations):
            raise SymbolicEquivalenceError("coordinate receipt is outside the declared bit range")
        assignments.append((assignment, expected))
    predicted = [0] * len(assignments)
    canonical_digest = hashlib.sha256(b"s3-coordinate-anf-masks-v1\0")
    byte_width = (n_vars + 7) // 8
    canonical_digest.update(n_vars.to_bytes(8, "little"))
    canonical_digest.update(len(equations).to_bytes(8, "little"))
    total_terms = 0
    maximum_degree = 0
    per_link = []
    ring = None
    for link, polynomial, ring in symbolic_links(instance, max_monomials=max_monomials,
                                                coefficient_cache_limit=coefficient_cache_limit):
        rows = coordinate_anf_sets(polynomial, field.degree, max_coordinate_terms=max_coordinate_terms)
        link_terms = sum(map(len, rows))
        for local_coordinate, expected_row in enumerate(rows):
            coordinate = link * field.degree + local_coordinate
            _canonical_row_digest(canonical_digest, expected_row, byte_width)
            equation = equations[coordinate]
            if not isinstance(equation, list):
                raise SymbolicEquivalenceError(f"coordinate {coordinate} is not a monomial list")
            if len(equation) != len(expected_row):
                raise SymbolicEquivalenceError(f"coordinate {coordinate} has {len(equation)} terms; exact S3 has {len(expected_row)}")
            for term in equation:
                if not isinstance(term, list) or any(type(variable) is not int or not 0 <= variable < n_vars
                                                     for variable in term):
                    raise SymbolicEquivalenceError(f"coordinate {coordinate} has an invalid monomial")
                if term != sorted(set(term)):
                    raise SymbolicEquivalenceError(f"coordinate {coordinate} monomial is not squarefree/canonical")
                mask = sum(1 << variable for variable in term)
                if mask not in expected_row:
                    raise SymbolicEquivalenceError(f"coordinate {coordinate} has an unexpected or repeated monomial {term}")
                expected_row.remove(mask)
                maximum_degree = max(maximum_degree, len(term))
            if expected_row:
                raise SymbolicEquivalenceError(f"coordinate {coordinate} omits exact S3 monomials")
            total_terms += len(equation)
        for index, (assignment, _) in enumerate(assignments):
            predicted[index] |= ring.evaluate(polynomial, assignment) << (link * field.degree)
        per_link.append({"link": link, "field_monomials": len(polynomial), "coordinate_terms": link_terms})
    if instance.get("term_occurrences") != total_terms:
        raise SymbolicEquivalenceError("term_occurrences disagrees with exact equation term count")
    if instance.get("max_degree") != maximum_degree or maximum_degree > 3:
        raise SymbolicEquivalenceError("max_degree disagrees with exact S3 coordinates")
    for index, (_, expected) in enumerate(assignments):
        if predicted[index] != expected:
            raise SymbolicEquivalenceError(f"recorded coordinate receipt {index} disagrees with exact S3")
    return {
        "schema": "exact_s3_symbolic_anf_equivalence_v1",
        "instance_id": instance.get("instance_id"), "n_vars": n_vars,
        "extension_degree": field.degree, "equation_count": len(equations),
        "term_occurrences": total_terms, "max_degree": maximum_degree,
        "equations_sha256": equations_digest,
        "canonical_anf_sha256": canonical_digest.hexdigest(),
        "canonical_digest_encoding": "v1 header; n_vars/u64le; n_rows/u64le; each row count/u64le then sorted fixed-width little-endian variable masks",
        "equivalence": "every coordinate monomial set matches exactly",
        "field_irreducibility_verified": True,
        "all_boolean_assignments_covered": True,
        "coordinate_receipts_checked": len(assignments), "links": per_link,
        "largest_field_polynomial": ring.largest_polynomial if ring else 0,
        "coefficient_products_computed": ring.coefficient_products if ring else 0,
        "coefficient_product_cache_hits": ring.coefficient_cache_hits if ring else 0,
        "verification_seconds": time.perf_counter() - started,
        "roots_solved": False, "performance_claim": False,
    }


def verify_file(path, *, expected_file_sha256=None, **budgets):
    raw = Path(path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if expected_file_sha256 is not None and digest != expected_file_sha256:
        raise SymbolicEquivalenceError("input file SHA-256 differs from the frozen preregistered input")
    result = verify_instance(json.loads(raw), **budgets)
    result.update({"file": str(path), "file_sha256": digest,
                   "frozen_file_hash_checked": expected_file_sha256 is not None,
                   "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--output")
    parser.add_argument("--max-monomials", type=int, default=1_000_000)
    parser.add_argument("--max-coordinate-terms", type=int, default=20_000_000)
    parser.add_argument("--coefficient-cache-limit", type=int, default=500_000)
    args = parser.parse_args()
    results = [verify_file(path, max_monomials=args.max_monomials,
                           max_coordinate_terms=args.max_coordinate_terms,
                           coefficient_cache_limit=args.coefficient_cache_limit)
               for path in args.files]
    output = json.dumps({"schema": "exact_s3_symbolic_anf_panel_v1", "files": results}, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(output)
    print(output)
    return 0


if __name__ == "__main__":
    # Support execution by pathname as well as ``python -m``.
    sys.exit(main())
