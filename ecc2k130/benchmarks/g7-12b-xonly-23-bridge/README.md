# Sparse-bridge x-only rho toward 12 B/s on G7

Preregistered successor experiment on the local AWS g7.2xlarge / RTX PRO 4500,
GPU UUID GPU-827e82b5-4739-d339-d95e-b7214337553e, at 165 W. The immutable
timing reference remains build/ecc2k130-local-packed, SHA-256
8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02. The goal is
12 billion complete scalar iterations/s on one GPU.

## Map and hypothesis

The rejected [2]/[3] map generated only an index-two scalar subgroup. Retain its
x-only common path, but when `weight(x) mod 32 == 14`, apply the bridge
`P -> P + sigma^3(P)`. Otherwise bit one of the weight selects [2] or [3].
All selectors are Frobenius and negation invariant.

For `y^2+xy=x^3+1`, independent symbolic reduction gives

    x(P+sigma^3(P)) = A(x)^2 / (x B(x)^2)
    A=x^7+x^6+x^4+x^3+x+1
    B=x^6+x^5+x^4+x^3+x^2+x+1.

The bridge needs five extra field products only on its sparse selector. Its
scalar multiplier `1+s^3` is a non-square; lcm of its order with ord(2) is
`ell-1`, so the transition generators cover the full scalar group. Freeze
the selector and formulas before implementation.

## Gates

Use an isolated copy of the rejected candidate and a new checkpoint version.
Before timing require: exhaustive P13/P19 and dense F23/F131 independent formula
proofs; GPU polynomial-form oracle for all three maps; built-in tests; complete
GPU x state against CPU replay; at least 64 DP endpoint replays; bidirectional
resume; partial populations; zero drops; and walk memcheck, initcheck and
synccheck. Preserve all failures.

Screen against the selected binary in three interleaved paired repetitions at
524,288 workers, B16, 1,024 steps and four launches: 34,359,738,368 complete
updates per sample. Advance only if the paired log-ratio Student-t 95% interval
is wholly above one. Report raw and collision-adjusted throughput. Promotion
requires five fresh timing and DP34 pairs, matching canonical-x corpora, both
intervals above one, and a larger matched collision cohort. No hardware, power,
driver, service, cloud-resource or publication changes.

## Preregistered spill-recovery follow-up

The first bridge build spilled 72 store / 60 load bytes per thread and regressed
to 5.249 B/s versus 6.311 B/s selected. Test exactly one structural change in a
fresh source copy: move bridge numerator/denominator construction into an
ECC_BIG noinline device helper. The bridge selector, formulas, common path,
batch layout and all build flags remain fixed. Rebuild and inspect resources
before timing. Repeat built-in, partial DP replay and sanitizers if the device
code changes. Screen only if the walk's spill traffic falls. Retain the initial
negative screen.

## Preregistered low-live-range bridge follow-up

The outlined helper worsened walk spills to 80 store / 64 load bytes and is
rejected without timing. In another fresh copy of the initial bridge, keep the
bridge inlined but construct B in one accumulator, emit its denominator, then
derive A from B while recomputing x^5 and x^6/x^7. This raises only the rare
bridge from five to seven products, while removing four simultaneously live
P131 temporaries. All selectors, formulas and common-path operations remain
fixed. Advance to validation only if spill traffic falls below the initial
72/60 bytes.

## Preregistered non-promotable diagnostic

Both valid bridge implementations regress near 0.83x. Run one paired diagnostic
sample of the already rejected zero-spill [2]/[3] binary versus selected. This
binary fails scalar-group generation and can never be promoted; its timing is
used only to attribute the regression to the base rational map or bridge code.
No confidence interval or goal claim may use this diagnostic.
