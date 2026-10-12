#!/usr/bin/env python3
"""A zero-compute readout of the m = 4 per-target solve-cost growth already in the repository.

Source: research/chain_split_order_20260924/tables.md (the frozen evidence of
research/notes/ecc2k130/RESEARCH_CHAIN_SPLIT_ORDER.md).  Its unit is word operations of the
Groebner decomposition oracle on the chained S3 system (m*l + (m-2)*n unknowns), summed over the
cell's targets, satisfiable and refuted mixed.  The 'candidate' column is the best engine in the
tree (interleaved split order + linear elimination); the 'reference' column is the engine before it.

The survey's exponent audit (decomp_exponent_model.py, Part B) says an m = 4 subspace arm can
have an exponent below rho's 1/2 only if the per-trial solve grows as 2^(c n) with c < 1/4 at
l ~ n/4.  This script reads the only m = 4 costs on record at two field sizes and prints the
implied per-n growth.  It is a two-size readout with confounds listed below, not a fit.
"""
import math
import statistics

# (cell, n, l, unknowns, targets, reference_total, candidate_total) -- copied from tables.md
ROWS = [
    ("K_0/2^9 m=4",            9,  6, 42, 8,     2_868_312,     917_450),
    ("K_1/2^9 m=4",            9,  6, 42, 8,     3_199_188,   1_168_831),
    ("K_1/2^11 m=4",          11, 10, 62, 4,    99_422_940,   2_983_641),
    ("K_0/2^15 m=4",          15,  4, 46, 8,    21_561_148,   5_648_706),
    ("K_1/2^15 m=4",          15,  4, 46, 8,   345_384_853,  39_537_587),
    ("K_0/2^15 m=4 div1",     15,  4, 46, 8,   149_499_912,  29_609_302),
    ("K_0/2^15 m=4 div2",     15,  4, 46, 8,   436_103_847,  64_781_197),
    ("K_1/2^15 m=4 div1",     15,  4, 46, 8,    28_328_427,   5_878_761),
    ("K_1/2^15 m=4 div2",     15,  4, 46, 8,   205_293_663,  33_838_553),
]

print(f"{'cell':<20} {'n':>3} {'l':>3} {'unk':>4} {'ref/target':>12} {'cand/target':>12}")
per = {}
for cell, n, l, unk, t, ref, cand in ROWS:
    per.setdefault(n, []).append((ref / t, cand / t))
    print(f"{cell:<20} {n:>3} {l:>3} {unk:>4} {ref / t:>12.3g} {cand / t:>12.3g}")

print()
for arm, k in (("reference", 0), ("candidate", 1)):
    m9 = statistics.median(v[k] for v in per[9])
    m15 = statistics.median(v[k] for v in per[15])
    c = math.log2(m15 / m9) / (15 - 9)
    print(f"{arm:>9}: median per-target word ops n=9 {m9:.3g}, n=15 {m15:.3g};"
          f"  implied growth {c:.2f} bits per unit n  (threshold for an exponent below 1/2: < 0.25)")

rho15 = 0.5 * (math.log2(math.pi) + 14 - math.log2(60))
print(f"\nContext: matched rho on a cofactor-2 subgroup at n = 15 is about 2^{rho15:.1f} = "
      f"{2 ** rho15:.0f} group operations for the WHOLE logarithm.")
print("Confounds: the n = 9 cells run l = 6 (x = l/n = 0.67) and the n = 15 cells l = 4 (x = 0.27);")
print("at a fixed x = 1/4 the n = 9 systems would be smaller and cheaper, so the per-n growth at")
print("fixed x is, if anything, steeper than printed.  Word operations are not group operations;")
print("the conversion is a constant and does not move a growth rate.  Two sizes fix one rate and")
print("no curvature, and the targets mix satisfiable and refuted systems in different proportions.")
