# K_0 n=83 paired arithmetic stage gate

Registered before execution on 2026-10-04. This is a component timing
diagnostic, not a complete IC candidate or one-target DLP run.

Build the pinned K_0 curve, its standard dimension-8 polynomial-subspace
base, and the public cofactor-projected subgroup-usable base. Take exactly
the first 64 points in the base's deterministic order. The target is the
public `T001` point in `gate-m83-T001.json`. Verify it lies on the curve and
in the declared subgroup before any timings. Construct a complete 2,080-pair
index, retaining all exceptional sums.

Compare all unordered pair sums made by scalar `point_add` and
`batch_add_fixed` byte for byte. After one warmup of each method, time each
three times, alternating method order by round. Compare exact four-summand
query answers on T001 using scalar residual additions and batched residual
additions through the same pair index; time each three times with the same
warmup and alternating order. Retain all raw durations, correctness status,
and witness or no-witness outcome. Use a release build. The host is
unisolated, so all CPU ratios are exploratory. A 64-point subset is too small
to estimate natural relation yield or establish F6/F4/F5 or IC/rho speed.
