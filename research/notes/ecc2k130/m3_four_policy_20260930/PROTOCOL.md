# Frozen degree-7 m3 four-policy held-out PDP gate

Status: preregistered before any target hit or rank outcome in this panel.
This closes a specific gap between the earlier degree-7 four-policy m2
replication and the source/transported-only m3 chain. It is a small same-Q
diagnostic for factor-base policy and a complete three-summand pair/residual
solver; it is not a measured ECC2K-130 attack. [CONFIG.json](CONFIG.json)
freezes all source hashes, toy parameters, seeds, stream labels and caps.

Use the pinned `GF(2^21)` source `y²+xy=x³+1`, order-421 subgroup and
cofactor 4,988. Recompute all eight order-7 twist-kernel lines; select the
same lexicographically first nonhorizontal line and codomain as the saved
pilot. Construct the normalized forward and dual maps and independently
check their subgroup identities on the generator and Q. Derive one nonzero
planted challenge scalar from SHA-256 of the frozen secret label modulo 421;
the scalar is audit-only and must never enter the relation solver.

Build two source and two **independently sampled native leaf** bases. For
each seed/curve, hash the frozen base label with `seed`, literal role
`source` or `leaf`, and trial to a 21-bit x; skip repeated x, inspect both
rational lifts, apply the cofactor, discard O and duplicate projected points.
Classify candidates by the source signed-Frobenius orbit; for leaf points use
the charged normalized dual pullback. Take exactly two useful points from
each of the first four sorted nonzero source orbit IDs, giving eight
distinct points per base with identical occupancy. The transported base is
the forward image of the source base; the pullback is the normalized dual
image of the native base. Stop a scan at 4,096 x trials and preserve an
incomplete-base failure. No base selection may inspect target outcomes.

Use the ten complete source signed-Frobenius orbit IDs to split targets:
holdout A accepts the first five sorted IDs, B the other five. Generate
512 accepted targets per holdout from the frozen new SHA-256 `(u,v)` labels,
with `u mod 421`, `v` in 1..420, and `T=[u]G+[v]Q`; reject infinity and
targets outside the assigned holdout, charging rejected draws to the next
accepted target. Cap each holdout at 8,192 draws. Evaluate the **same**
accepted coefficients on the leaf generator/Q. Check forward target
covariance separately as an audit, not as candidate work.

For each of two seeds, two holdouts and four policies—original,
transported, descendant native, exact pullback—build the complete 36-entry
unordered pair table on the eight points. For each target, scan the eight
third factors in base order, look up the residual in the pair table, and
take the first lexicographic pair in the first successful residual. Record
every miss, the first witness, pair-entry multiplicity, tested third
factors and group-operation counts for all 512 targets, including the tail
after rank. A successful relation contributes the eight base-log
coefficients and `−v` challenge coefficient modulo 421; track row rank to
width nine, solve at the first full rank, and verify **every** recovered
base log and `[d]G=Q`. Independently replay all base scans, targets,
first witnesses, ranks, scalars, map covariance and phase accounting in a
fresh process with separately written enumeration and rank logic. Use a
bit-polynomial point-law reference on G, Q, every selected base point and
the first/last two targets of each holdout. Source/transported and
native/pullback hit, first-witness and rank trajectories must agree exactly.

The eight-point base has only `C(8+2,3)=120` unordered triples, so the
uniform-target one-shot support ceiling is `120/421 < 0.286`; ordered tuple
counts must not be presented as hit probability. The primary cold ledger
is the vector of field multiplications (including inversion internals),
squarings, inversion calls, group additions, modular row operations and
process CPU. Charge field setup, kernel search, forward/dual setup actually
used, generator/Q, source/leaf base scans, cofactor work, leaf dual
classification, target draws/rejections and codomain targets through rank,
pair construction, every failed residual scan through rank, row reduction,
solve and scalar verification. Audit covariance and the post-rank outcome
tail are separate. Report a field-multiplication ratio only as a toy
operation diagnostic, not an end-to-end attack speedup.

All 16 cells must either reach independently verified rank nine or retain
their exact failure/rank receipts; no seed, base size, target prefix or
policy may be adjusted from outcomes. A native operation-count advantage
requires lower fully charged field multiplications to verified rank in all
four paired seed/holdout cells, no rank or hit regression, no higher
modular-row burden, and no failed control. Otherwise retain the mixed or
negative result. A cap exceedance or missing archive is `CENSORED`, a
semantic mismatch is `FAIL`, and a complete unsupported or rank-deficient
panel is a measured negative. The producer has 120 seconds and 512 MiB;
the independent verifier has 180 seconds and 512 MiB. Preserve every raw
case, source hash, failure, host/cost receipt and decision in a follow-on
PR; update the canonical scoreboard. No n131 PDP-yield, full-DLP or rho
crossover claim follows from this toy panel.
