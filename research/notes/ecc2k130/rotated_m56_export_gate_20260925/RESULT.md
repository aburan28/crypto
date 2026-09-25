# Rotated m5/m6 exporter semantic gate: exact rational reference, no polynomial export

**Decision: PASS for a bounded correctness/capacity gate.** The corrected,
preregistered source enumerated every one of #767's 3,125 and 117,649
labelled rational point tuples and an equal-count repeated-`F0` control on
each rung. All 64 exact fixed-torsion counts for the 16 frozen `Q` points
match #767's independently produced complete point histograms. Every
rational tuple passed the Kosters–Yeo trace morphism check; every affine
recursive-S3 prefix and exact-target terminal check vanished, and each
identity/inverse step was classified separately. The output preserves every
accepted x-mask's rational lift multiplicity and a checked full point-tuple
witness. A fresh deterministic replay passed after the local run and the
focused CI replays the committed evidence. This admits a later direct/chain
**exporter correctness comparison**, not a solver timing, full index-calculus
cost, or ECC2K-130 speed claim.

The first frozen attempt at head `572ebd1` failed before any tuple loop:
its preflight assumed all four formal x-values in every two-dimensional
slice lift to rational points. At n13, one per slice does not. The exact
failure is in `evidence/failure_0/`; source and protocol were corrected and
re-frozen at `98c83a0` with `FROZEN.json` SHA-256
`be1cd1db1556a4a36cb79040e576adcbb1b6829e3d34934876d4a427d90c3d8a`
before a new output path was used. No failed attempt was relabelled or
combined with the successful run.

| Exact semantic arm | Frozen rotated tuples | Rational x-mask tuples / formal masks | Projected support: rotated / repeated F0 | Fixed coset count match | Direct x bits | Raw S3-chain bits / rows before exception auxiliaries | Gate child wall / CPU / peak RSS |
|:--|--:|--:|:--|:--|--:|:--|:--|
| n13, m5, d2 | 3,125 | 243 / 1,024 | 1,591 / 61 | 32/32 | 10 | 49 bits / 52 rows | 0.273 / 0.268 s / 28.0 MiB |
| n19, m6, d2 | 117,649 | 4,096 / 4,096 | 62,389 / 377 | 32/32 | 12 | 88 bits / 95 rows | 20.887 / 20.313 s / 95.8 MiB |

The direct x-bit counts describe the *layout*, not a generated S6/S7 ANF.
The 49/88-bit chain counts include `(m−2)n` intermediate x bits but exclude
infinity, sign, denominator, field-equation and SAT auxiliaries. The current
`u64` Boolean monomial representation admits the raw n13 49-bit layout and
rejects the raw n19 88-bit one. No direct resultant, recursive polynomial
system or solver input was emitted in this PR; the 15-second/2-GiB portfolio
pilot of [#763](https://github.com/aburan28/crypto/pull/763) has **not**
started. FES also cannot take an unquadraticized cubic S3 descent as a
quadratic system.

The trace check is exact for this rational-point PDP. On this curve
`φ(P)=Tr(x(P))`, with `φ(O)=0`; normal beta has trace one, so the XOR of
the two-bit masks' parities equals `φ(sum)`. The n13 rotated tuple parity
classes were 1,441 even and 1,684 odd; the n19 classes were 58,825 and
58,824. `Q∈[4]E` has trace zero, and the four fixed torsion branches have
parities `[0,0,1,1]`. All 64 `Q+T` branches obeyed that rule. A fixed-T
branch can discard the opposite parity, while the all-T projected PDP still
needs both classes. This is **not** a measured 2× overall reduction. The
linear trace condition filters rational candidates and is not necessarily
implied by higher Semaev ideals with nonrational algebraic lifts, as
[Kosters–Yeo, Prop. 4.2 and Remark 4.8](https://arxiv.org/html/1503.08001v3)
make explicit.

The chain necessity audit checked 12,500 n13 and 588,245 n19 prefix
additions. Their classifications were respectively `ordinary=12,240,
inverse-to-infinity=130, identity-prefix=130` and `ordinary=583,344,
inverse-to-infinity=2,451, identity-prefix=2,450`; no doubling occurred on
these frozen tuples. Every affine S3 transition vanished; an infinity step
has no x-coordinate and must remain an explicit branch in a future chain
exporter. The nine n13 and ten n19 exact point hits across all positive
`Q+T` branches passed terminal S3. These checks show necessity for rational
witnesses, not sufficiency of an algebraic S3 chain. Every raw S6/S7 or chain
root must still be independently lifted to signed rational points and
re-added to its exact `Q+T` target.

The rotated nonzero summand spaces have only x=0 in common, so exchanging
points between slots is not a symmetry of the factor-domain constraints.
The preregistered structural control found full projected labels, outside
the 16 frozen solver targets, that disappear if a naive first-slot order is
imposed across **all** four torsion cosets. Excluding the eight frozen target
labels on each rung, polynomial-basis integer `x0≤x1` loses 422 n13 and
31,236 n19 supported projected labels; a slot-local two-bit-mask order loses
257 and 15,624. The first eligible exact controls are `R=(154,4452)`,
`Q=(7821,5622)` at n13 and `R=(56,271694)`, `Q=(511055,446447)` at n19.
The committed receipts include each entire point tuple, full sum, torsion
shift and swapped tuple. Each swap preserves the group sum but places both
nonzero points outside their new `F_i` slots; no ordered witness exists for
that projected label. This is a concrete **invalid-admission control** for
copying WDSat's same-base `m!` ordering into the rotated system, not a solver
speed measurement. [Trimoska–Ionica–Dequen, §4](https://eprint.iacr.org/2019/313.pdf)
assumes interchangeable slots.

The repeated-`F0` control yielded **no witness in any torsion branch for any
of the 16 frozen Q targets**, including all eight rotated-planted positives.
Thus these fixed targets cannot support a same-Q positive solver-time
comparison against repeated `F0`. A later repeated-vs-rotated experiment must
preselect targets in shared support (or report the repeated arm as UNSAT),
and charge work per useful verified relation and independent rank row.
Support opportunity alone cannot be multiplied by the lost same-base
symmetry reduction to infer a PDP or whole-attack speedup.

The successful n13/n19 children ran sequentially at 12:43:58–12:44:19 UTC.
Their runner UTC intervals, including interpreter startup and final JSON
serialization, were 0.307 and 20.983 seconds; the child wall values in the
table stop just before serialization and are used only as diagnostics.
The frozen caps were 600 seconds and 512 MiB **per arm**. Their combined
native counters in each child's receipt include both rotated and repeated
semantic controls: n13 had 37,655 curve additions, 247,905 field
multiplications, 111,409 squares and 2,647 inversions; n19 had 1,848,716
curve additions, 13,863,410 multiplications, 6,238,885 squares and 134,939
inversions. Square calls include a multiplication inside the bit-serial
implementation, so these are native event counters, not independent field
operation-equivalents or calibrated group-addition S. Full raw files, exact
hashes and replay commands are in [evidence/README.md](evidence/README.md).

The next gate is a separately frozen direct/chain polynomial or circuit
exporter, with complete denominator/infinity/rational-lift semantics and a
model verifier against these exact mask sets on all 16 point targets.
Before importing the public ECC2K-130 challenge points, verify the public
field/model representation and exact P/Q curve and q-order membership, then
archive the source coordinates and import hashes.
Nothing in this toy corpus measures a target decomposition at n131.
