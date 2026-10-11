#!/usr/bin/env python3
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
base = ROOT / "build/g7-xonly-doubling-src"
bridge = ROOT / "build/g7-xonly-23-bridge-src"
dst = ROOT / "build/g7-xonly-sparse-bridge-unified-src"
if dst.exists():
    raise SystemExit(f"destination exists: {dst}")
shutil.copytree(base, dst)
for name in ("ref.h", "solver.h", "packedengine.cuh"):
    shutil.copy2(bridge / "include" / name, dst / "include" / name)

p = dst / "include/ref.h"
s = p.read_text()
old = """#if ECC_PACKED_XONLY_BRIDGE3\n        if ((hw & 31) == 14) return addPt(p, frob(p, 3));\n#endif\n        const Point twice = dbl(p);\n        return ((hw >> 1) & 1) ? addPt(twice, p) : twice;\n"""
new = """#if ECC_PACKED_XONLY_BRIDGE3\n        if ((hw & 31) == 14) return addPt(p, frob(p, 3));\n#endif\n#if ECC_PACKED_XONLY_DOUBLE_ONLY\n        return dbl(p);\n#else\n        const Point twice = dbl(p);\n        return ((hw >> 1) & 1) ? addPt(twice, p) : twice;\n#endif\n"""
assert s.count(old) == 1
p.write_text(s.replace(old, new))

p = dst / "include/solver.h"
s = p.read_text()
old = """#if ECC_PACKED_XONLY_BRIDGE3\n            if ((hw & 31) == 14) out.counts[2]++;\n            else\n#endif\n            out.counts[(hw >> 1) & 1]++;\n"""
new = """#if ECC_PACKED_XONLY_BRIDGE3\n            if ((hw & 31) == 14) out.counts[2]++;\n            else\n#endif\n#if ECC_PACKED_XONLY_DOUBLE_ONLY\n            out.counts[0]++;\n#else\n            out.counts[(hw >> 1) & 1]++;\n#endif\n"""
assert s.count(old) == 1
p.write_text(s.replace(old, new))

p = dst / "include/packedxonly23.cuh"
s = p.read_text()
needle = "#if !ECC_PACKED_BLOCK_INVERSE || !ECC_PACKED_POLY_STATE || \\\n"
helper = r"""static __device__ __noinline__ PolynomialPair sparseBridge3Rational131(P131 x) {
    // x(P + sigma^3(P)) = A(x)^2 / (x B(x)^2).
    const P131 x2 = squarePolynomial131(x);
    const P131 x4 = squarePolynomial131(x2);
    const P131 x3 = mulPolynomial131(x, x2);
    const P131 x5 = mulPolynomial131(x, x4);
    const P131 x6 = mulPolynomial131(x2, x4);
    const P131 x7 = mulPolynomial131(x, x6);
    const P131 one{{1, 0, 0, 0, 0}};
    const P131 a = add131(add131(add131(add131(add131(x7, x6), x4), x3), x), one);
    const P131 b = add131(add131(add131(add131(add131(add131(x6, x5), x4), x3), x2), x), one);
    return PolynomialPair{squarePolynomial131(a),
        mulPolynomial131(x, squarePolynomial131(b))};
}

"""
assert s.count(needle) == 1
s = s.replace(needle, helper + needle)

old = r"""                // [2]P: x'=(x+1/x)^2. Store each denominator and the
                // product before it; backward recovery needs no numerator.
                const P131 denominator = x;
                if (localSlot) {
                    store(p.pchain, slot, tid, p.threads, prod);
                    prod = mulPolynomial131(prod, denominator);
                } else {
                    prod = denominator;
                }
#if ECC_PACKED_LAST_SLOT_CACHE
                if (localSlot == lastLocalSlot131) lastDenominator = denominator;
                else
#endif
                store(denominators, slot, tid, p.threads, denominator);
"""
new = r"""                // Common [2]P path: x'=(x+1/x)^2. The sparse bridge
                // keeps its rational numerator in the unused x-only y field.
                const bool sparseBridge = ECC_PACKED_XONLY_BRIDGE3 &&
                    (hw & 31) == 14;
                P131 denominator;
                if (sparseBridge) {
                    const PolynomialPair rational = sparseBridge3Rational131(x);
                    store(p.y, slot, tid, p.threads, rational.first);
                    denominator = rational.second;
                } else {
                    denominator = x;
                }
                if (localSlot) {
                    store(p.pchain, slot, tid, p.threads, prod);
                    prod = mulPolynomial131(prod, denominator);
                } else {
                    prod = denominator;
                }
                // Compact field tails use only bits 0..2. Bit 3 is metadata
                // added after multiplication and cleared before arithmetic.
                P131 storedDenominator = denominator;
                if (sparseBridge) storedDenominator.v[4] |= 8u;
#if ECC_PACKED_LAST_SLOT_CACHE
                if (localSlot == lastLocalSlot131)
                    lastDenominator = storedDenominator;
                else
#endif
                store(denominators, slot, tid, p.threads, storedDenominator);
"""
assert s.count(old) == 1
s = s.replace(old, new)

old = r"""#if ECC_PACKED_LAST_SLOT_CACHE
                const P131 denominator = localSlot == lastLocalSlot131
                    ? lastDenominator : load(denominators, slot, tid, p.threads);
#else
                const P131 denominator = load(denominators, slot, tid, p.threads);
#endif
#if ECC_PACKED_XONLY_DOUBLE_ONLY
#if ECC_PACKED_SHARED_X_SLOTS
                const P131 x = localSlot >= firstSharedXLocalSlot131
                    ? loadSharedX131(sharedXLow, sharedXTail,
                        (localSlot - firstSharedXLocalSlot131) * ECC_THREADS + threadIdx.x)
                    : load(p.x, slot, tid, p.threads);
#else
                const P131 x = load(p.x, slot, tid, p.threads);
#endif
                P131 inverseX;
                if (localSlot) {
                    const PolynomialPair pair = mulPolynomialPair131(inv,
                        load(p.pchain, slot, tid, p.threads), denominator);
                    inverseX = pair.first;
                    inv = pair.second;
                } else {
                    inverseX = inv;
                }
                const P131 nx = squarePolynomial131(add131(x, inverseX));
"""
new = r"""#if ECC_PACKED_LAST_SLOT_CACHE
                P131 denominator = localSlot == lastLocalSlot131
                    ? lastDenominator : load(denominators, slot, tid, p.threads);
#else
                P131 denominator = load(denominators, slot, tid, p.threads);
#endif
#if ECC_PACKED_XONLY_DOUBLE_ONLY
                const bool sparseBridge = (denominator.v[4] & 8u) != 0;
                denominator.v[4] &= 7u;
                P131 inverseDenominator;
                if (localSlot) {
                    const PolynomialPair pair = mulPolynomialPair131(inv,
                        load(p.pchain, slot, tid, p.threads), denominator);
                    inverseDenominator = pair.first;
                    inv = pair.second;
                } else {
                    inverseDenominator = inv;
                }
                P131 nx;
                if (sparseBridge) {
                    nx = mulPolynomial131(inverseDenominator,
                        load(p.y, slot, tid, p.threads));
                } else {
#if ECC_PACKED_SHARED_X_SLOTS
                    const P131 x = localSlot >= firstSharedXLocalSlot131
                        ? loadSharedX131(sharedXLow, sharedXTail,
                            (localSlot - firstSharedXLocalSlot131) * ECC_THREADS +
                                threadIdx.x)
                        : load(p.x, slot, tid, p.threads);
#else
                    const P131 x = load(p.x, slot, tid, p.threads);
#endif
                    nx = squarePolynomial131(add131(x, inverseDenominator));
                }
"""
assert s.count(old) == 1
p.write_text(s.replace(old, new))
