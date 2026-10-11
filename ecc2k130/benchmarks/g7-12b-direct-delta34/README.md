# Direct polynomial delta-3/4 toward 12 B/s on G7

Preregistered on the local AWS g7.2xlarge / RTX PRO 4500, GPU UUID
GPU-827e82b5-4739-d339-d95e-b7214337553e, at 165 W. The goal is 12 billion
complete ECC2K-130 scalar updates/s on one GPU. The immutable timing reference
is build/ecc2k130-local-packed, SHA-256
8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02.
Generic work is sqrt(n/262); full-DLP S is null.

## Hypothesis

Retain the validated two-jump map P -> P + sigma^j(P), j in {3,4}, selected by
bit one of the normal-basis X weight. Its multipliers generate the full scalar
group, and the frozen 6,000-target planted-DLP cohort measured a 1.0130803718
collision-work ratio to the selected eight-jump map.

Keep one polynomial-to-normal X conversion for the invariant weight and DP
record. Form both coordinate deltas directly in polynomial basis:

    D = x + x^(2^j)
    E = y + y^(2^j).

Compute powers by three polynomial squarings, plus a fourth result selected
when j=4. This removes normal Y conversion, paired normal-basis Frobenius
routing, and both delta conversions back to polynomial form. All products,
batch inversion, affine equations, state, DP, restart, checkpoint and launch
geometry remain fixed. Repeated polynomial squaring may cost more than the
removed linear networks; compilation and complete timing decide.

## Gates

Use an isolated copy of the validated two-jump source and a default-off mode 3.
Before timing require: all 117 built-in checks; independent delta equality on
all 131 basis vectors, edge and dense inputs for both coordinates and j values;
complete GPU state and DP equality against CPU replay; bidirectional resume;
partial populations; zero drops; and walk memcheck, initcheck and synccheck.
The disabled/indexed client remains the semantic reference.

Screen selected, indexed two-jump and direct-delta two-jump in three interleaved
paired repetitions, 524,288 workers, B16, 1,024 steps and four launches:
34,359,738,368 complete updates/sample. Advance only if the paired Student-t
95% log-ratio interval versus selected is wholly above one after dividing by
the frozen 1.0130803718 collision-work ratio. Promotion requires five fresh
benchmark and DP34 pairs, both intervals above one, exact records, zero drops
and the full validation suite. Preserve regressions. No hardware, power,
driver, service, cloud-resource or publication changes.
