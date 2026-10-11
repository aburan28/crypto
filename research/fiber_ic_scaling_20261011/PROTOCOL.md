# Protocol: does fiber-aware index calculus ever cross rho end to end?

Frozen 2026-10-11, before any session of this study ran.  This follows the
conductor-gap rounds
([1](../isogeny_conductor_gap_ic_20261010/README.md),
[2](../isogeny_conductor_gap_round2_20261011/README.md),
[3](../isogeny_conductor_gap_round3_20261011/README.md)).

Their levers each cut index-calculus cost on cofactor curves:
- the fiber-closed pair table (round 1);
- the corrected closure and partial fibers (round 2);
- the 2-torsion quotient (round 3).

None was compared with rho across sizes.  This study asks the end-to-end
question: **on a fixed-cofactor family, does the best fiber-aware IC
approach rho as the curve grows, or recede from it?**

Result class: whole-method cold cost against matched rho, one target per
run (AGENTS.md §2, §8).  Synthetic public curves, no deployed-curve claim.

## Derivation, stated before measuring

A meet-in-the-middle IC with `m = 2` uses:
- a factor base of `B` abscissae;
- a negation-folded pair table of `π ≈ B²` sums;
- a fiber of size `s ≤ h` adjoined (closed table, `π·s` entries).

On a curve of order `N = h·r`:
- A probe hits with probability `≈ 2π·s/N`.
- About `B` relations are needed, so the relation phase is `≈ c·B·N/(2π s)`
  probes, plus a set-up of `π·s` additions.
- The total `≈ N/(B s) + B² s` is minimised at `B ∝ (N/s²)^{1/3}`.
- That gives `IC ∝ N^{2/3} s^{−1/3}`.

Rho (negation) costs `∝ (N/h)^{1/2}`.  Hence, with the whole fiber `s = h`:

    IC / rho  ∝  N^{1/6} · h^{1/6}

The prediction is that the ratio **grows** with the field: by `2^{1/6} ≈ 1.12`
per bit of `N`, so by 4× from 18 to 30 bits.  The fiber improves the
constant against the blind table (`s = 1`), by `h^{−1/3}`, but never the
exponent.  No size above the measured range crosses rho.

**Predictions.**
- **P-S1 (the question).**  For each cofactor, take per curve the best IC
  arm (the envelope over the `B` sweep).  The slope of `log₂(S_IC/S_rho)`
  against `log₂ N` has a 95 percent interval that contains `1/6` and
  excludes 0.
  - A slope interval at or below 0 would mean the gap is closing.  That is
    the outcome that would make this a speedup lead.
- **P-S2.**  At matched size, the cofactor-72 family's envelope ratio
  exceeds the cofactor-4 family's, by about `(72/4)^{1/6} ≈ 1.6`.
- **P-S3.**  Blind (`s = 1`, `B = 128`) over the closed envelope is at least
  `h^{1/3}` (1.6 at `h = 4`, 4.2 at `h = 72`) wherever `B = 128` is the
  blind arm's best size.

## Inputs

- **Curves.**  Random prime-field curves `y² = x³ + ax + b` with
  `#E = h·r`, `r` prime, `h ∈ {4, 72}`.
  - Two curves per size at 18, 22, 26 and 30 bits.
  - Found by `examples/fiber_scaling_curves.rs` (order by baby-step
    giant-step on the Hasse interval, confirmed on a second point).
  - `ecbench`'s `prime_explicit` re-checks each curve on load.
- **Arms.**
  - `rho.negation`.
  - `ic.pipeline` with `prime-abscissa:size=B` and
    `mitm-fiber:fiber=closed,negation_folded=1,fiber_gen=coset` for
    `B ∈ {16, 32, 64, 128, 256}`.
  - The blind table (`fiber=blind`) at `B = 128`.
- **Measurement.**  3 targets per curve, 2 measured rounds, alternate
  order, 900 s timeout.  One session per cofactor, each audited with
  `ecbench verify --replay 12`.
- **Analysis.**  `examples/fiber_scaling_tables.rs` computes:
  - the per-curve envelope;
  - the least-squares slope of `log₂(S_IC/S_rho)` on `log₂ N`, with a 95
    percent interval from the residuals;
  - P-S2 and P-S3.

  `ecbench bound fit` records each fixed-`B` arm's exponent beside it.

## Stop

- A session cut by a reboot is kept under `interrupted/` and rerun once.
- No arm or size is added after a session has been read.
- If P-S1's interval reaches 0 or below, that is reported as the lead.  A
  frozen challenge against rho at larger sizes then follows in a new
  protocol.

## Visual record

Each session gets a run report with a graph and a PDF under
`reports/run-<UTC>/`.  The README carries the envelope-against-size graph
and a PDF.

## Amendment 1 (2026-10-11, after reading `scaling_h72`; before running the new session)

On six of the eight cofactor-72 curves, the cheapest closed arm was the
smallest factor base swept (`B = 16`).  The derivation puts the optimum at
`B ∝ (N/h²)^{1/3}`, about 4 to 10 at 18 to 26 bits.  So the envelope there
is not resolved, and the slope fit is biased low.

Session `scaling_h72_small_b` reruns the same curves and targets with the
closed table at `B ∈ {4, 8}` and with rho.  It is analysed in two ways:
- with the envelope taken over both sessions;
- with the original `scaling_h72` alone, kept as the pre-registered test.

This amendment was added after reading the data and is reported as such.
