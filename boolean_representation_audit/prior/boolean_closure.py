"""Small exact Boolean-algebra audit; synthetic systems only, n <= 10.

Degree prioritizes work; it never deletes a term. This is not a scalable
solver, a novel basis algorithm, or an application-specific integration.
"""
from heapq import heappop, heappush

from baseline.bounded_sparse_gf2 import (
    BoundedMonomials, ResourceLimitError, SparseEchelon,
)


def normalize(n, generators):
    if type(n) is not int or not 1 <= n <= 10:
        raise ValueError("the audit supports only 1 <= n <= 10")
    output = []
    for generator in generators:
        terms = set()
        for mask in generator:
            if type(mask) is not int or not 0 <= mask < (1 << n):
                raise ValueError("invalid monomial mask")
            terms.symmetric_difference_update((mask,))
        output.append(sorted(terms))
    return output


def compute(n, generators, *, strategy="frontier", max_submissions=100_000,
            max_pending_terms=1_000_000, max_row_terms=100_000,
            max_total_terms=1_000_000):
    """Return a partial span on budget exhaustion, never a false completion.

    `frontier` closes accepted rows under variable multiplication.
    `exhaustive` enumerates all monomial multiples of input generators.
    Both use exactly the same recovered representation and reducer.
    """
    generators = normalize(n, generators)
    if strategy not in ("frontier", "exhaustive"):
        raise ValueError("unknown strategy")
    if max_submissions < 1 or max_pending_terms < 1:
        raise ValueError("budgets must be positive")
    space = BoundedMonomials(n, n)
    echelon = SparseEchelon(space, max_row_terms=max_row_terms,
                            max_total_terms=max_total_terms)
    inputs = [space.row_from_masks(g) for g in generators]
    certificate = {"n": n, "rows": [], "origins": []}
    stats = {"submitted": 0, "products_formed": 0,
             "peak_pending_terms": 0, "max_scheduled_degree": 0}
    queue = []
    pending_terms = serial = 0

    def multiply(row, multiplier):
        stats["products_formed"] += 1
        product, dropped = space.multiply_row(row, multiplier)
        if dropped:
            raise AssertionError("exact mode must never discard terms")
        return product

    def enqueue(row, origin):
        nonlocal pending_terms, serial
        if not row:
            return
        if pending_terms + len(row) > max_pending_terms:
            raise ResourceLimitError("pending term budget exceeded")
        degree = space.unrank(row[-1]).bit_count()
        heappush(queue, (degree, serial, row, origin))
        serial += 1
        pending_terms += len(row)
        stats["peak_pending_terms"] = max(stats["peak_pending_terms"], pending_terms)
        stats["max_scheduled_degree"] = max(stats["max_scheduled_degree"], degree)

    def submit(row, origin):
        if stats["submitted"] >= max_submissions:
            raise ResourceLimitError("submission budget exceeded")
        stats["submitted"] += 1
        if not echelon.add(row):
            return None
        pivot = echelon.pivots[next(reversed(echelon.pivots))]
        certificate["rows"].append([space.unrank(col) for col in pivot])
        certificate["origins"].append(origin)
        return pivot

    complete, reason = False, None
    try:
        if strategy == "exhaustive":
            for i, row in enumerate(inputs):
                for mask in range(1 << n):
                    product = multiply(row, mask)
                    submit(product, ["input_multiple", i, mask])
        else:
            for i, row in enumerate(inputs):
                enqueue(row, ["input", i])
            while queue:
                _, _, row, origin = heappop(queue)
                pending_terms -= len(row)
                pivot = submit(row, origin)
                if pivot is None:
                    continue
                parent = len(certificate["rows"]) - 1
                for variable in range(n):
                    product = multiply(pivot, 1 << variable)
                    # An unchanged product is already represented exactly.
                    if product != pivot:
                        enqueue(product, ["row_variable", parent, variable])
        complete = True
    except ResourceLimitError as error:
        reason = str(error)
    stats.update(rank=echelon.rank, stored_terms=echelon.stored_terms,
                 packed_payload_bytes=echelon.packed_payload_bytes,
                 xor_steps=echelon.xor_steps,
                 peak_row_terms=echelon.peak_row_terms,
                 pending_terms_at_stop=pending_terms)
    return {"complete": complete, "reason": reason,
            "stats": stats, "certificate": certificate}


# Independent verifier: columns use raw monomial masks, represented by bits
# in a Python integer. It does not use the recovered indexer or sparse reducer.
def _vector(masks):
    value = 0
    for mask in masks:
        value ^= 1 << mask
    return value


def _multiply_vector(value, mask):
    result = 0
    while value:
        bit = value & -value
        monomial = bit.bit_length() - 1
        result ^= 1 << (monomial | mask)
        value ^= bit
    return result


def _reduce(value, pivots):
    while value:
        column = value.bit_length() - 1
        if column not in pivots:
            break
        value ^= pivots[column]
    return value


def verify(n, generators, result, *, require_complete=True):
    """Check containment in both directions and multiplicative closure.

    The externally supplied n and generators are the statement being proved;
    the certificate cannot replace them. Raises ValueError on a failed proof.
    """
    generators = normalize(n, generators)
    cert = result["certificate"]
    if cert["n"] != n or len(cert["rows"]) != len(cert["origins"]):
        raise ValueError("malformed certificate")
    if require_complete and not result["complete"]:
        raise ValueError("incomplete result")
    inputs = [_vector(g) for g in generators]
    rows, pivots = [], {}
    for masks, origin in zip(cert["rows"], cert["origins"]):
        canonical = normalize(n, [masks])[0]
        if len(canonical) != len(masks):
            raise ValueError("noncanonical row")
        value = _vector(masks)
        kind = origin[0]
        if kind == "input" and len(origin) == 2:
            i = origin[1]
            if type(i) is not int or not 0 <= i < len(inputs):
                raise ValueError("invalid input origin")
            source = inputs[i]
        elif kind == "input_multiple" and len(origin) == 3:
            i, mask = origin[1:]
            if (type(i) is not int or not 0 <= i < len(inputs)
                    or type(mask) is not int or not 0 <= mask < (1 << n)):
                raise ValueError("invalid input multiple")
            source = _multiply_vector(inputs[i], mask)
        elif kind == "row_variable" and len(origin) == 3:
            parent, variable = origin[1:]
            if (type(parent) is not int or not 0 <= parent < len(rows)
                    or type(variable) is not int or not 0 <= variable < n):
                raise ValueError("invalid earlier-row origin")
            source = _multiply_vector(rows[parent], 1 << variable)
        else:
            raise ValueError("invalid origin")
        if _reduce(value ^ source, pivots):
            raise ValueError("provenance failed")
        remainder = _reduce(value, pivots)
        if not remainder:
            raise ValueError("zero or dependent certified row")
        pivots[remainder.bit_length() - 1] = remainder
        rows.append(value)
    # Partial results certify only that their rows belong to the input ideal.
    if not require_complete:
        return {"rank": len(rows), "verified": "containment_only"}
    if any(_reduce(g, pivots) for g in inputs):
        raise ValueError("generator containment failed")
    for row in rows:
        for variable in range(n):
            if _reduce(_multiply_vector(row, 1 << variable), pivots):
                raise ValueError("variable closure failed")
    return {"rank": len(rows), "verified": "exact_ideal"}


def solutions(n, generators):
    """Exhaustive independent oracle; includes assignments with zeros."""
    generators = normalize(n, generators)
    return [assignment for assignment in range(1 << n)
            if all(sum((mask & assignment) == mask for mask in g) % 2 == 0
                   for g in generators)]


def equal_spans(left, right):
    def basis(result):
        pivots = {}
        for masks in result["certificate"]["rows"]:
            row = _reduce(_vector(masks), pivots)
            if row:
                pivots[row.bit_length() - 1] = row
        return pivots
    a, b = basis(left), basis(right)
    return len(a) == len(b) and all(not _reduce(row, b) for row in a.values())
