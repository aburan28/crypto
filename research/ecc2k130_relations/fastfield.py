#!/usr/bin/env python3
"""Fast `F_2^131` arithmetic for the ECC2K-130 relation sweeps.

`scripts/ecc2k130_point_decomposition.py` carries the reference
implementation of this field.  It is bit-at-a-time, which is the right
choice for the few thousand operations its cost model needs, and the
wrong one for a sweep over every pair of a 3.6k-point factor base:
measured on this container it costs 17.8us per multiply, 3.59ms per
inversion and 1.50ms per half-trace, which puts 6.7M pairs at about
nine hours.

This module keeps the same field -- same reduction polynomial, same
element encoding (an `int` whose bit `i` is the `z^i` coefficient) --
and makes three changes:

  * squaring, half-trace and every Frobenius power are `F_2`-linear, so
    each becomes a table of 17 byte-indexed lookups XORed together;
  * multiplication uses a 4-bit comb window and reduces once, rather
    than reducing after every bit;
  * inversion is amortised over a batch by Montgomery's trick, which
    turns `k` inversions into one inversion and `3k` multiplications.

Every routine is checked against the reference implementation in
`test_fastfield.py`; the two agree element for element, which is what
makes the speed admissible rather than a new and unverified field.
"""

from __future__ import annotations

N131 = 131
# x^131 + x^13 + x^2 + x + 1, the challenge reduction polynomial.
IRR131 = (1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1


class FastGF2m:
    """`F_2[z]/(irr)`, with byte tables for the linear maps."""

    def __init__(self, deg: int = N131, irr: int = IRR131):
        self.deg, self.irr = deg, irr
        self.mask = (1 << deg) - 1
        self.tail = irr ^ (1 << deg)
        self.taps = tuple(i for i in range(deg) if (self.tail >> i) & 1)
        self.nbytes = (deg + 7) // 8

        # z^i for i >= deg, reduced -- the basis of every reduction.
        self._pow = [1 << i for i in range(deg)]
        cur = self._fold(1 << deg)
        for _ in range(deg, 2 * deg):
            self._pow.append(cur)
            cur = self._fold(cur << 1)

        self.trace_bits = self._trace_bits()
        self._sqr_tbl = self._linear_table(self._sqr_slow)
        self._ht_tbl = self._linear_table(self._half_trace_slow)
        self._frob_tbl: dict[int, list] = {}

    # -- reduction -------------------------------------------------------

    def _fold(self, a: int) -> int:
        """Reduce `a` of any degree, folding the high part through the tail."""
        while a > self.mask:
            hi = a >> self.deg
            low = a & self.mask
            r = 0
            for t in self.taps:
                r ^= hi << t
            a = low ^ r
        return a

    def reduce(self, a: int) -> int:
        return self._fold(a)

    # -- table construction ----------------------------------------------

    def _linear_table(self, f):
        """Tabulate an `F_2`-linear `f` as `nbytes` x 256 byte lookups."""
        img = [f(1 << i) for i in range(self.deg)]
        tbl = []
        for byte in range(self.nbytes):
            row = [0] * 256
            base = byte * 8
            for v in range(1, 256):
                low = v & (v - 1)          # v with its lowest bit cleared
                bit = (v ^ low).bit_length() - 1
                idx = base + bit
                row[v] = row[low] ^ (img[idx] if idx < self.deg else 0)
            tbl.append(row)
        return tbl

    def _apply(self, tbl, a: int) -> int:
        r = 0
        for row in tbl:
            r ^= row[a & 0xFF]
            a >>= 8
        return r

    # -- slow reference maps, used only to build the tables --------------

    def _sqr_slow(self, a: int) -> int:
        s = 0
        i = 0
        while a:
            if a & 1:
                s ^= self._pow[2 * i]
            a >>= 1
            i += 1
        return s

    def _half_trace_slow(self, c: int) -> int:
        assert self.deg % 2 == 1
        t, acc = c, c
        for _ in range((self.deg - 1) // 2):
            acc = self._sqr_slow(self._sqr_slow(acc))
            t ^= acc
        return t

    def _trace_bits(self) -> int:
        bits = 0
        for i in range(self.deg):
            v, s = 1 << i, 0
            for _ in range(self.deg):
                s ^= v
                v = self._sqr_slow(v)
            assert s in (0, 1), "trace must land in F_2"
            bits |= (s & 1) << i
        return bits

    # -- the field operations --------------------------------------------

    def mul(self, a: int, b: int) -> int:
        # 4-bit comb: build the 16 multiples of `a`, then one reduction.
        t0 = 0
        t1 = a
        t2 = a << 1
        t3 = t2 ^ a
        t4 = a << 2
        t5 = t4 ^ a
        t6 = t4 ^ t2
        t7 = t4 ^ t3
        t8 = a << 3
        tbl = (t0, t1, t2, t3, t4, t5, t6, t7, t8,
               t8 ^ a, t8 ^ t2, t8 ^ t3, t8 ^ t4,
               t8 ^ t5, t8 ^ t6, t8 ^ t7)
        acc = 0
        shift = 0
        while b:
            acc ^= tbl[b & 0xF] << shift
            b >>= 4
            shift += 4
        return self._fold(acc)

    def sqr(self, a: int) -> int:
        return self._apply(self._sqr_tbl, a)

    def half_trace(self, c: int) -> int:
        """`t` with `t^2 + t = c`, valid for odd `deg` and `Tr(c) = 0`."""
        return self._apply(self._ht_tbl, c)

    def trace(self, a: int) -> int:
        return (a & self.trace_bits).bit_count() & 1

    def solve_artin_schreier(self, c: int):
        if self.trace(c):
            return None
        return self.half_trace(c)

    def frobenius(self, a: int, k: int = 1) -> int:
        k %= self.deg
        if k == 0:
            return a
        tbl = self._frob_tbl.get(k)
        if tbl is None:
            def f(v, k=k):
                for _ in range(k):
                    v = self._sqr_slow(v)
                return v
            tbl = self._linear_table(f)
            self._frob_tbl[k] = tbl
        return self._apply(tbl, a)

    def inv(self, a: int) -> int:
        assert a, "zero has no inverse"
        # a^(2^deg - 2) by square-and-multiply on the exponent's bit pattern
        r, base, e = 1, a, (1 << self.deg) - 2
        while e:
            if e & 1:
                r = self.mul(r, base)
            base = self.sqr(base)
            e >>= 1
        return r

    def batch_inv(self, xs):
        """Montgomery's trick: one inversion for the whole list.

        Zero entries are passed through as zero rather than raising, so
        callers can invert a column that legitimately contains one.
        """
        n = len(xs)
        out = [0] * n
        prefix = [0] * n
        run = 1
        for i, x in enumerate(xs):
            prefix[i] = run
            if x:
                run = self.mul(run, x)
        if run == 0:
            return out
        acc = self.inv(run)
        for i in range(n - 1, -1, -1):
            x = xs[i]
            if x:
                out[i] = self.mul(acc, prefix[i])
                acc = self.mul(acc, x)
        return out

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))
