# Seven-arm held-out relation control protocol (V3)

The public source curve is `y²=x³−21x+37` over `F_211²`, the degree-three
quotient has `j=160`, and the prime-order subgroup has order 229. The
public target positions 64–67, signed eight-column factor base, private
checker scalars, random-stream seeds, 1,832-trial cap, three rounds, and
same-target rho controls are frozen by `FOUR_TARGET_PUBLIC_V1.json`,
`FOUR_TARGET_CHECKER_V1.json`, and `HELDOUT_FIVE_ARM_SPEC_V2.json`.

**Hypothesis.** The quadratic source and quotient `S_3` decomposition
oracles have the same exact pair support as direct enumeration and the
resultant oracle on this cell; their complete per-target costs may differ.
The exhaustive prerequisite is an independent 229-case check on both curves,
including zero-pair cases. If any pair set differs, stop the panel and
investigate that case. A matching set permits the panel, but does not imply
matching performance.

**Arms.** Interleave `chi-resultant`, `direct-pair`, `source-s3`,
`quotient-s3`, `source-bsgs`, `rho-negation`, and `rho-negation-aa`, with
round-dependent cyclic order. Each child runs from process start through
setup, query, row solve where applicable, and source and quotient scalar
verification. The native process meter records wall time, CPU time, peak
RSS, exit status, and timeout. The 512 MiB Darwin memory ceiling is
observed; the 300-second child wall cap is enforced. Preserve all failed
children and retry only under a versioned protocol.

**Acceptance.** A separate Rust checker must replay all 84 scalars on both
curves, every relation random stream and witness, row rank, BSGS result,
and rho-returned scalar; all children must exit successfully. The primary
cost table is whole-child metered CPU and counted-group operations divided
by `sqrt(229)`, with exact arm/control ratios. These are implementation
diagnostics until a nonoverlapping field-operation conversion and isolated
wall admission establish portable `S`. A/A must be reported beside rho.

The comparison remains incomplete with respect to the direct inverse-Vélu
and full-image relation arms; their matching factor-base construction on
this specific cell and complete pricing are separate prerequisites.
