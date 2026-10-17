# Four follow-on gates for isogeny and torsion fibers

Frozen on 2026-10-16 before any run under this protocol. This protocol
extends [NEXT_TESTS.md](NEXT_TESTS.md), the immutable inverse-Vélu V2
receipt, and the different-`j` order-229 public fixture. New implementations
and failed runs receive versioned records; neither prior fixture nor prior
measurement is rewritten.

## Boundary and common rules

The cold one-target ECDLP question is recovery of the unknown scalar in
`Q=[k]G`, with a separate checker replaying `[k]G=Q`. The generic
two-class collision floor in counted group operations is
`sqrt(pi/4)=0.886226925` in `S=work/sqrt(r)`. The references are the
same-target source negation rho and source BSGS, with rho A/A for wall
noise. All five arms include validation, construction, every rejected
attempt, relation collection, rank and row solve, map transfer, target
recovery, and final verification. No answer scalar enters a solver input.
One-target cold and 64-target setup reuse are distinct questions.

The existing order-1201 public file has SHA-256
`24f2c56a9ba56f91c954f38b7752edf24fee5ed027b49b1ac0747771945099c4`.
Use its exact 64 targets and five arms in three rotated rounds, at
prefixes 1, 4, 16, 64, with `8r` relation trials per target, 32 rho
restarts of at most `20r` steps, 300 seconds and 512 MiB per child.
The V2 result and its independent checker are the equality reference.

## Gate 1: complete native price

Hypothesis: the direct inverse construction has a smaller complete cold
cost than the full source census, while a comparison against direct rho
and BSGS remains open. Instrument field operations separately inside
and outside counted group calls. The ledger must cover addition,
multiplication, inversion, exponentiation and square-root branches,
polynomial gcd/split work, row arithmetic, table/hash probes and
inserts, allocations, process setup and I/O, rejected candidates,
relation and final scalar replay. Count a primitive exactly once.

Calibrate each primitive with pinned native kernels against a source
group addition on Darwin and an independent Linux host. Freeze the
factors and conservative intervals before ranking children. Rerun the
60-child panel and compare its witness and answer stream bit for bit
with V2. An independent Rust checker must reconstruct every phase,
sum, target, witness and scalar. Publish cold N1 and N64 `S` intervals
for all five arms, plus physical CPU, wall and RSS. If any category is
missing or the two host factors conflict, leave complete portable `S`
unset and retain the vector and host-specific intervals.

## Gate 2: nonfactorable χ decomposition

Hypothesis: an eliminated χ-coordinate relation equation can change
the number of verified, independent rows per complete unit of work.
Use the smallest committed different-`j` cell,
`F_(211^2)`, `r=229`, source `y²=x³-21x+37`, degree-three kernel
`{O,(3,1),(3,-1)}`, `j=77` and `160`. The committed public fixture
`../../different_j_order_batch_20261015/PUBLIC_V1.json` has SHA-256
`f6c3126b4d95e3c457c9e552471bc2cdff322de0e0a9de9fbba701d7c05897ba`.

Define `z=(y+1)/(y-1)` where finite, `d=z-1`, `n=z+1`. On this source
the point constraint is
`F(X,z)=d²(X³-21X+37)-n²=0`. Form the cleared
chord or tangent constraint for `P1+P2=R`, then eliminate both source
abscissae by resultants. Intersect the eliminated polynomial with the
frozen signed factor-base χ polynomial and solve its roots without a
pair-sum table or a loop over every signed factor-base pair. Verify
every candidate by full source addition. Treat infinity, vertical
chords, `y=1`, repeated summands, torsion shifts, repeated resultant
factors and false roots as explicit branches and counters.

Before any DLP timing, compare the complete unordered pair multiset
and multiplicities for every point of the order-229 subgroup against
independent direct pair enumeration; include zero-result cases. Any
mismatch stops the solve panel. Then select four public targets outside
the prior 64 by the same domain-separated point-hash ordering rule,
store their scalars only in a checker file, and freeze both hashes.
Run complete relation collection, rank-nine row recovery and scalar
replay for χ, direct S3, quotient-invariant, inverse-Vélu, full image,
source BSGS, rho and rho A/A. Relation arms use the same eight columns,
random streams, `8r` trial cap and 512 MiB memory ceiling. Charge
every resultant, root, false candidate and verification. If a control
cannot use exactly the same base or equation, report the mismatch and
do not make a paired-yield claim for that row. The decisive output is
useful verified rank or relations per complete work, then a verified
same-target cold solve; polynomial root count alone is diagnostic.

## Gate 3: host replay and isolation

Build the inverse runner and checker from the committed source on a
Linux x86-64 host. Execute the 60-child five-arm panel on the identical
order-1201 public file; retain source, binary, input, stdout, stderr,
checker, host and artifact hashes. Compare exact answers, witnesses,
operation counters and costs with Darwin. A hosted VM can supply an
independent L1 performance replay, not the L2 wall admission.

For L2, run at least five interleaved paired rounds on a reserved whole
physical core, one NUMA node when applicable, with affinity readback,
thermal, IRQ, host load, frequency, RSS and A/A receipts. Quote a wall
ratio only if the candidate/control paired 95% interval excludes one
and the A/A interval, under the ecbench isolation and noise rules.
Record unavailable runner access as an access or resource blocker, not
as a mathematical outcome. No credential or quota bypass is authorized.

## Gate 4: direct inverse order trend

Hypothesis: the direct inverse construction's complete cold cost changes
with subgroup order, beyond the existing one-order measurement.
Screen ascending primes `1238 <= p <= 10007`, `p=1 mod 3`, in the exact
`y²=x³-21x+37` family. Require nonsingular different-`j` source and
destination, exact order-three kernel, equal independently counted
base-field orders, and a prime-to-`p` factor `r` of multiplicity one in
`#E(F_(p²))`. Select the first eligible cell, by ascending `p`, in
each disjoint order bin `[2400,4800)`, `[4800,9600)`, `[9600,19200)`,
and `[19200,38400)`. Retain every rejection and a bin with no eligible
cell. The screen ceiling is 3600 seconds and 512 MiB; a cap yields a
censored screen, not an empty bin.

Audit the selected field, `j`, kernel, quotient, subgroup generator,
sampled map homomorphism and full three-point fibers where the 512 MiB
and 300-second child caps permit. Freeze a fresh public point by the
domain-separated point-hash rule, then 64 distinct points with known
scalars kept only in the checker. At each cell run direct inverse,
full census where feasible, source BSGS, signed rho and rho A/A in
three interleaved rounds at N1 and N64. Apply `8r`, 32 by `20r`, 300
seconds and 512 MiB uniformly. Preserve timeout and memory censorship
for every arm; do not silently remove a difficult cell. Replay every
returned scalar and relation witness independently. Four checked cells
permit a descriptive trend; an exponent change requires complete
portable cost, uncertainty and admitted isolation over a credible size
range.

## Decision and records

Each gate gets its own run ID, immutable raw receipts, source and binary
hashes, independent replay, report, SVG/PNG graph and PDF. Performance
runs also copy the visual record into `benchmarks/` with a README row.
The existing problem record remains open until all complete-price and
proof obligations it names are met. This protocol is a plan, not a
receipt for any of the four runs.
