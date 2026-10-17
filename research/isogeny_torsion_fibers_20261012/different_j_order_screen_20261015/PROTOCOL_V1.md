# Frozen larger-order different-j candidate screen

**Status: preregistered before implementation or execution.** This is
Gate D from the [60-target reuse result](../different_j_reuse_20261015/NEXT_GATE.md).
It is a bounded construction test, not a DLP timing or speed claim.

## Hypothesis and selector

The rational-kernel family used for the verified order-61 different-j
pair may contain a larger subgroup. Enumerate **every ascending prime**
`p` with `61 <= p <= 271` and `p = 1 (mod 3)`. Use the canonical
quadratic field `Fp2::new(p)`, source
`E: y²=x³-21x+37`, and `T=(3,1)`. Require a nonsingular source,
`T` on the source, `T != O`, and `[3]T=O`. Use the explicit degree-three
Vélu construction with kernel `{O,T,-T}` and destination coefficients
`a'=-9a-30x_T²`, `b'=-27b-70x_T³-42a*x_T`. Require a nonsingular
destination and distinct source/destination `j` values. Count both
base-field group orders independently and require equality. Compute
`n2=(p+1)²-t²` from `t=p+1-#E(Fp)` and factor it exactly. A candidate
qualifies if `n2` has a prime factor `r>=211`, `r!=3`, with exponent
one. Select the **first qualifying p**, then its largest qualifying r.
Do not stop scanning after selection: retain an output row for every
eligible prime, including each rejection and its first failed check.

For every row that passes the algebraic predicates, check the kernel
maps to identity, maps of 64 seeded source-point pairs land on the
destination, and `φ(P+Q)=φ(P)+φ(Q)` for those pairs. A failed map
check rejects the row. Record field identity, both models and j values,
`n1`, `n2`, complete factorization, selected r, cofactor and operation
counters. The deterministic sample seed is `0x20261015 XOR p`. Check
the output with a separate native Rust reader that recomputes primality,
selection, group-order factor arithmetic, j distinction and all stated
map witnesses. A successful screen only admits a pair to a subsequent
full solve; it does not establish a relation-yield or cost advantage.

## Limits, receipts and decision

Build the native Rust example in release mode, then commit source and
binary SHA-256, exact command, host identity and this protocol before
the screen runs. Send only public synthetic curve parameters and points
to the screen. Set a 120-second whole-screen wall ceiling and 512 MiB
RSS ceiling. Preserve stdout, stderr, exit status and any timeout; a
failed or interrupted screen is censored, never a mathematical miss.
Do not rewrite a run; corrections get versioned spec and output paths.

If a qualifying pair appears, commit its full identity and a follow-on
protocol for at least 64 distinct public subgroup targets with source,
projected-image, destination-direct, full-image, direct source BSGS,
negation rho and A/A controls. Freeze public target order, caps and all
binary hashes before running it. Independent scalar and relation replay,
complete cold/process accounting, portable counted
`S=total operations/sqrt(r)`, isolated measurement and independent-host
replay remain separate acceptance gates. The generic two-class rho
floor for that later comparison is `sqrt(pi/4)=0.886226925`; matched
direct BSGS and rho are references, not assumed wins.

If no pair qualifies, report only a bounded miss for this exact family
and `p<=271`; publish the rejection census and the next bound before
another search. Neither result decides the direct χ-coordinate
resultant route, which has its own [pending gate](../different_j_reuse_20261015/NEXT_GATE.md).
