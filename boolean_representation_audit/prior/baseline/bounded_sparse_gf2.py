"""Degree-bounded monomial indices and packed sparse GF(2) rows.

Standalone algebra code; no external solver or application integration.
Monomials are squarefree Boolean masks. A degree cap is explicit and lossy.
"""

from array import array
from bisect import bisect_right
from math import comb


class ResourceLimitError(MemoryError):
    """A configured row or matrix budget was exceeded."""


class BoundedMonomials:
    """Zero-based graded-colex combinatorial index, built without a lookup table.

    The first C(n,0) indices have degree zero, the next C(n,1) degree one,
    and so on. Within degree k, positions i_1 < ... < i_k have colex rank
    sum_{j=1}^k C(i_j,j). Variables are numbered 0 through n-1.
    """

    def __init__(self, n: int, max_degree: int):
        if not isinstance(n, int) or not isinstance(max_degree, int):
            raise TypeError("n and max_degree must be integers")
        if n < 1 or not 0 <= max_degree <= n:
            raise ValueError("require n >= 1 and 0 <= max_degree <= n")
        self.n, self.max_degree = n, max_degree
        starts = [0]
        for k in range(max_degree + 1):
            starts.append(starts[-1] + comb(n, k))
        self.starts = tuple(starts)
        self.count = starts[-1]
        if self.count > 1 << 64:
            raise ValueError("packed 64-bit indices cannot represent this space")
        self.typecode = 'I' if self.count <= 1 << 32 else 'Q'

    def rank(self, mask: int) -> int:
        if not isinstance(mask, int) or mask < 0 or mask.bit_length() > self.n:
            raise ValueError("monomial outside variable universe")
        degree = mask.bit_count()
        if degree > self.max_degree:
            raise ValueError("monomial exceeds degree cap")
        result, j = self.starts[degree], 1
        while mask:
            lowest = mask & -mask
            result += comb(lowest.bit_length() - 1, j)
            mask ^= lowest
            j += 1
        return result

    def unrank(self, index: int) -> int:
        if not isinstance(index, int) or not 0 <= index < self.count:
            raise IndexError("monomial index out of range")
        degree = bisect_right(self.starts, index) - 1
        remainder = index - self.starts[degree]
        mask, upper = 0, self.n - 1
        for j in range(degree, 0, -1):
            lo, hi = j - 1, upper
            while lo < hi:
                middle = (lo + hi + 1) // 2
                if comb(middle, j) <= remainder:
                    lo = middle
                else:
                    hi = middle - 1
            mask |= 1 << lo
            remainder -= comb(lo, j)
            upper = lo - 1
        assert remainder == 0
        return mask

    def row_from_masks(self, masks) -> array:
        return self.row_from_indices(self.rank(mask) for mask in masks)

    def row_from_indices(self, indices) -> array:
        ordered = sorted(indices)
        result = array(self.typecode)
        position = 0
        while position < len(ordered):
            index = ordered[position]
            if not isinstance(index, int) or not 0 <= index < self.count:
                raise IndexError("row index out of range")
            end = position + 1
            while end < len(ordered) and ordered[end] == index:
                end += 1
            if (end - position) & 1:
                result.append(index)
            position = end
        return result

    def multiply_row(self, row: array, monomial_mask: int) -> tuple[array, int]:
        """Boolean multiplication with an explicit count of dropped high terms.

        No claim about rank or ideal membership may ignore a nonzero drop count.
        """
        if not isinstance(monomial_mask, int) or monomial_mask < 0 or monomial_mask.bit_length() > self.n:
            raise ValueError("multiplier outside variable universe")
        if row.typecode != self.typecode:
            raise TypeError("row uses incompatible packed index width")
        output, dropped = [], 0
        for index in row:
            combined = self.unrank(index) | monomial_mask
            if combined.bit_count() > self.max_degree:
                dropped += 1
            else:
                output.append(self.rank(combined))
        return self.row_from_indices(output), dropped


def xor_sorted(left: array, right: array, typecode: str) -> array:
    """Symmetric difference of two sorted, duplicate-free packed rows."""
    if left.typecode != typecode or right.typecode != typecode:
        raise TypeError("incompatible packed rows")
    result = array(typecode)
    i = j = 0
    while i < len(left) and j < len(right):
        if left[i] < right[j]:
            result.append(left[i]); i += 1
        elif right[j] < left[i]:
            result.append(right[j]); j += 1
        else:
            i += 1; j += 1
    result.extend(left[i:])
    result.extend(right[j:])
    return result


class SparseEchelon:
    """Incremental GF(2) echelon rank with explicit fill-in budgets."""

    def __init__(self, space: BoundedMonomials, *, max_row_terms=100_000,
                 max_total_terms=1_000_000):
        if max_row_terms < 1 or max_total_terms < 1:
            raise ValueError("budgets must be positive")
        self.space = space
        self.max_row_terms = max_row_terms
        self.max_total_terms = max_total_terms
        self.pivots: dict[int, array] = {}
        self.stored_terms = 0
        self.xor_steps = 0
        self.peak_row_terms = 0

    @property
    def rank(self): return len(self.pivots)

    @property
    def packed_payload_bytes(self):
        return self.stored_terms * array(self.space.typecode).itemsize

    def add(self, input_row: array) -> bool:
        if input_row.typecode != self.space.typecode:
            raise TypeError("incompatible packed row")
        if len(input_row) > self.max_row_terms:
            raise ResourceLimitError("input row exceeds term budget")
        if any(index >= self.space.count for index in input_row):
            raise IndexError("row index outside degree-bounded space")
        if any(input_row[i] >= input_row[i + 1] for i in range(len(input_row) - 1)):
            raise ValueError("row must be strictly increasing")
        row = input_row[:]
        while row and row[-1] in self.pivots:
            row = xor_sorted(row, self.pivots[row[-1]], self.space.typecode)
            self.xor_steps += 1
            self.peak_row_terms = max(self.peak_row_terms, len(row))
            if len(row) > self.max_row_terms:
                raise ResourceLimitError("elimination fill-in exceeded row budget")
        if not row:
            return False
        if self.stored_terms + len(row) > self.max_total_terms:
            raise ResourceLimitError("matrix term budget exceeded")
        self.peak_row_terms = max(self.peak_row_terms, len(row))
        self.pivots[row[-1]] = row
        self.stored_terms += len(row)
        return True
