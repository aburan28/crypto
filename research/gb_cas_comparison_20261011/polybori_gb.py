"""External reference only: run PolyBoRi (BRiAL, via passagemath) on one
exported Boolean system and print its wall time and whether the ideal is
the unit ideal.  It reads an msolve input written by
`examples/gb_cas_export.rs bool`, drops the field equations (implicit in
the Boolean ring) and calls PolyBoRi's `groebner_basis` with its default
strategy.  No result of this repository is computed here; this stands in
for Magma's Boolean-ring F4, which is not available in the container.

    python polybori_gb.py <root.ms>
"""
import sys
import time

from sage.all__sagemath_brial import BooleanPolynomialRing, ideal

lines = open(sys.argv[1]).read().split("\n")
names = lines[0].split(",")
assert lines[1].strip() == "2"
polys = [p.strip().rstrip(",") for p in "\n".join(lines[2:]).split(",\n")]
field_eqs = {f"{v}^2+{v}" for v in names}
polys = [p for p in polys if p and p not in field_eqs]
B = BooleanPolynomialRing(len(names), names=names, order="degrevlex")
gens = [B(p) for p in polys]
t0 = time.perf_counter()
G = ideal(gens).groebner_basis()
dt = (time.perf_counter() - t0) * 1000.0
print(f"polybori_ms {dt:.1f} basis_size {len(G)} unit_ideal {int(list(G) == [B(1)])}")
