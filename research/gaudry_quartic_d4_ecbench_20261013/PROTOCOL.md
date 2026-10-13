# Protocol: Gaudry's `k = 4` index calculus with the `D₄` oracle, end to end in ecbench

Registered 2026-10-13, before the measured session.  Preceded only by a
six-trial plumbing run of an unregistered scratch spec (same curve family,
`max_trials = 6`, exhausted, not cited) and by decomposition-rate counts
with `examples/gaudry_quartic_d4.rs --find 16 --hit-rate`, used to pick `p`.

## Question

`research/gaudry_quartic_d4_20261013` measured one decomposition (`C₄`)
and projected the whole method on §11.16's model.  Here the whole method
runs once, end to end, in ecbench: factor base, relation trials, linear
algebra mod `r`, verification of `[d]G = Q` by the runner, every phase in
the harness's unit `S` (group-addition equivalents / √r), against rho with
the negation map on the same target.

## What was added to ecbench

- Curve construction `fp4_torsion { p, seed, max_cofactor }`
  (`gaudry_quartic_d4::find_d4_curve`): `E/F_{p⁴}`, `F_{p⁴} = F_p[t]/(t⁴ − c)`,
  rational 2-torsion `T`, `#E = h·r`.  Points travel as packed words.
  Identity: `curve_id::extension` (the Rust form of `scripts/curve_id.py`'s
  `extension_id`; the slug agrees with Python's on this curve).  Registered.
- `CountedGroup` for `E(F_{p⁴})`; rho and the generic methods run on it.
- `ic.pipeline` on this curve with base `fghr-y-line` and oracle
  `gaudry-d4:macaulay=D`.  Pricing, all in counted `F_p` multiplications
  converted at the curve's own multiplications per addition (97 here):
  the oracle's solver work (`priced_by = counted`), its sign-checking
  additions (charged as additions), base lifts (square roots at their
  measured mean count), linear-algebra row operations (one multiplication
  each).  Nothing is wall-priced.

## The run (`specs/SPEC-k4-p29.json`)

`fp4_torsion(29, 1, 16)`: `#E = 708,208 = 16 · 44,263`, `log₂ r = 15.43`,
7 columns.  One target, one round, no warm-up, `isolation_required = L2`,
`--cpus auto` (one reserved core, one thread).  Arms: `rho.negation`
(reference), `ic.pipeline` + `fghr-y-line` + `gaudry-d4:macaulay=14`,
`max_trials = 20,000` (candidate), `rho.negation` A/A (control).

## Predictions

- E1. Every arm verifies `[d]G = Q`; the audit replays every run identically.
- E2. The IC arm needs 600–4,000 relation trials (decomposition rate
  0.003–0.008 on this curve, about `columns + 2` relations).
- E3. `S_IC` lies between `3·10⁸` and `10¹⁰`, and `S_IC / S_rho` between
  `10⁸` and `10¹⁰`; §11.16's model with the measured `C₄`, the measured
  rate and `h = 16` puts the ratio near `1.3·10⁹`.
- E4. The relation phase is at least 99 % of `S_IC`.

## Decision

This is the method run once at a size where it is expected to lose to rho
by nine orders of magnitude; the point is the whole pipeline working and
priced, and the model's inputs checked against an end-to-end count, not a
speed-up.  Class: **engineering**.  A miss on E3 by more than 10× either
way means §11.16's projection of this method needs correcting before it
is quoted.
