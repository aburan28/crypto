# Graded Boolean Macaulay reuse across affine-tail batches

## Hypothesis and exact invariant

The [sibling pivot screen](https://github.com/aburan28/crypto/pull/1202) is an
unmerged, source-pinned negative discovery at the time this protocol is
written: whole pivot traces change across branch restrictions, while the
matrix difference nearly fills the low-degree tail. Its observations are not
treated as an accepted downstream performance result. The next hypothesis
rests on a separate exact identity. Write each quadratic Boolean generator
as `f_j=q_j+a_j`, where `q_j` contains the degree-two terms and `a_j` is
affine. At Macaulay degree three, each multiplier has degree at most one.
Therefore the degree-three projection of `t*f_j` depends only on `q_j`:
`pi_3(t*(q_j+a_j))=pi_3(t*q_j)`. In the ordered matrix
`M_a=[H(q) | L(q,a)]`, all cubic columns `H(q)` are identical across systems
with the same quadratic core. The affine/constant columns may differ, cancel,
or change their rank.

Test whether compiling and eliminating `H(q)` once per batch, replaying its
exact row-operation schedule on each changing low block, and completing the
low-column elimination beats **the pointwise fastest correct matched full
matrix constructor and elimination**. The previous indexed-envelope 0/32
construction result cannot serve as the only reference: the merged
packed-direct and scalable ranked constructors are included. This is a
finite generic Boolean-matrix experiment, not a full solver or DLP test.

## Matrix and algorithm contract

Generate n-variable public quadratic systems at n=8/12/16/20/24 with m=n.
For each fixed seed and n, initialize the SplitMix64 transition specified in
`boolean_pivot_screen_20261003/PROTOCOL.md`, then give each generator 2n
distinct quadratic monomials using successive pairs of unequal variable
indices. Continue the PRNG stream to generate affine coefficients. No planted
solution, curve input, scalar, target, or external equation is admitted.

The degree bound is D=3. Every generator retains its nonempty quadratic core,
so its multiplier set is exactly all squarefree monomials of degree at most
one in ascending numeric order. Row labels are generator index then
multiplier. Keep labelled zero rows internally; remove zero output rows in the
same deterministic final step in every arm. Columns are degree-three
monomials in descending numeric order, then monomials of degree at most two
in descending numeric order. Boolean products use bitwise OR and XOR parity.
The independent oracle uses set-parity products and compares the exact
returned row space, rank, row labels and canonical output digest.

The candidate first constructs both fixed `H(q)` and fixed contributions to
the low block, then echelonizes H and records every row swap and XOR. Its
retained schedule, packed high echelon, fixed low contribution, coordinates,
and allocations are charged to **cold batch setup**. For each affine
assignment, it builds only changing low contributions, combines them with
the fixed low block, applies the exact high-block schedule, completes low
column elimination, materializes the full output and validates it. The
candidate must produce the same deterministic echelon as the fresh full
elimination, not merely the same rank. It falls back to a complete fresh
constructor if n, D, generator order, quadratic support, row labels, caps or
any schedule invariant differs. A repeat family and quadratic-support
escapes are retained as guardrails, with no timing selected from them.

Fresh controls in the same executable are: packed direct with dense lookup
where its `2^n` table fits the declared 64 MiB context cap; combinatorial
ranked lookup with dense intermediate rows; combinatorial ranked lookup with
sparse intermediate rows; and exact completed-matrix cache for the repeat
guardrail. Every control pays cold setup, construction, elimination, output
validation and destruction. At each primary cell, the denominator is the
pointwise fastest applicable correct control for each paired repetition.
No historical absolute time enters a new ratio. The dense lookup is
inapplicable if its context exceeds the cap; that condition is recorded, not
counted as a failure or omitted silently.

## Frozen workloads and gates

`protocol.json` fixes discovery seeds 20261004 and 3659031, untouched
holdout seeds 20261011 and 6931207, n=8/12/16/20/24, batches 1/2/8/32/64,
ten balanced arm-position repetitions and 4,000 paired bootstrap resamples.
Each batch starts with the all-zero affine tail. In `independent_affine`, each
later system draws a fresh constant and n linear coefficient bits per
generator from successive low bits of SplitMix64 outputs. In
`walk_affine`, each later system toggles one affine slot per generator,
selected by `next()%(n+1)`; slot 0 is the constant and slots 1..n are the
linear variables. Neither family changes any quadratic coefficient.
`repeat` keeps the first system unchanged. In `support_escape`, every fourth
later system toggles one selected quadratic monomial in generator zero; the
candidate must take a measured full-build fallback and match the oracle.

Discovery correctness requires exact outputs on every cell, exact fallback
and hit counts, unchanged source/protocol hashes, zero censored results, and
all resource receipts qualified. Discovery advances to fresh holdouts only
if the candidate's paired 95% lower timing bound exceeds **1.5** against the
pointwise fastest control in all eight n=12/16/20/24, batch-64 cells of the
two changing-affine families. Otherwise reject before holdouts and retain the
negative result. This 1.5 screen is a resource gate, not the dramatic claim.

The primary holdout claim requires a 95% paired-bootstrap lower bound above
**2.0 in all eight** of those same batch-64 cells, each above its own A/A
97.5th-percentile noise floor. No size, family, failed case or reference arm
may be dropped after timing. Batches 1/2/8/32, repeat and support-escape are
fully reported guardrails. Report setup, per-assignment construction,
high-schedule application, low reduction, output/validation, deterministic
GF(2) word operations, retained context bytes, process peak RSS, fallbacks,
and complete cold batch time. A one-size or construction-only win does not
pass the universal claim.

Run each phase under a single recorded Linux ARM64 core reservation through
the repository isolation controller. Admit a campaign only when other-process
CPU use is at most 10% and CPU PSI some avg10 is at most 5.0; retain every
preflight refusal. Arm order rotates; every arm receives a fresh cold batch,
and A/A pairs use the same reference path. Fixture and independent oracle
generation are common supplied-input work outside arm clocks but inside
whole-process receipts. No setup or failed attempt is subtracted from the
candidate. Caps: 4,096 rows, 8,192 columns, 64 MiB retained context and
20-minute campaign. Caps return UNKNOWN/censored, never a success or UNSAT.

All code, orchestration, tests, verification and evidence replay for this new
study must be native Rust except the repository's required CPU-isolation
controller and thin shell build/run commands. Commit the source before
timing, seal raw outputs and failed attempts with hashes, and use a separate
fresh holdout run only after discovery passes. The protocol creates no
speedup claim now. Full polynomial-solving cost, natural relation yield,
calibrated IC operations and rho ratio remain null until measured directly.
