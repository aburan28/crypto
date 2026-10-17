#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set heading(numbering: "1.")

= Strong-rho canonicalization, round 2 (2026-10-16)

*Study:* `research/rho_canon2_20261016/` · *Protocol:* `PROTOCOL.md` (declared before the changes; two addenda, one recorded deviation) · *Baseline:* `72843c83a` (`main`) · *Host:* one x86-64 cloud VM (Intel Xeon @ 2.10 GHz, 4 vCPU), no hardware PMU.

The strong rho is the reference every Koblitz `vs_rho` comparison prices index calculus against. This round makes its canonicalization cheaper without changing any walk. Every output line at $n$ = 41, 53 and 61 is identical to the baseline's, timing aside, and all 13 strong-rho unit tests pass with debug assertions on.

== Changes

+ `MulMod::mul` casts through `i64`: one instruction per conversion instead of a multi-instruction `u64`↔`f64` sequence.
+ `rotate` drops its `k = 0` branch. A `debug_assert` pins $v < 2^n$, $k <= n$.
+ `canonicalize` reuses the rotations the least-rotation search and the sign test already computed.
+ The orbit multipliers $plus.minus lambda^k$ are precomputed with their signs.
+ Shoup's fixed-operand product for $a dot w mod m$: exact for $a, w < m < 2^63$.
+ `canonicalize` is inlined into `solve`.

== Instructions (measured, Callgrind, $n$ = 53, 4 fixtures)

#table(columns: 3, align: (left, right, right),
  [*arm*], [*instructions*], [*per walk step*],
  [baseline], [2,847,169,163], [938],
  [changes 1–4], [2,580,125,642], [850],
  [changes 1–5], [2,557,912,731], [843],
  [*final, 1–6*], [*2,369,298,800*], [*781*],
)
The final build runs 16.8% fewer instructions per step (1.20×), with an identical output digest in every arm.

#figure(image("instructions-per-step.svg", width: 90%), caption: [Instructions per walk step by arm (`runs/callgrind/`).])
#figure(image("instructions-by-function.svg", width: 90%), caption: [Where a step's instructions go, baseline vs final.])

== Native time (measured)

A/A (baseline × 6), then 5 interleaved A/B rounds. 32 processes under `isolated_bench`, none contended.

#table(columns: 5, align: (left, right, right, right, left),
  [*curve*], [*baseline ns/step*], [*final ns/step*], [*ratio*], [*A/A spread, verdict*],
  [`icv1-f2m53-…`], [130.8], [111.5], [1.173], [20.6%: not established],
  [`icv1-f2m61-…`], [139.7], [123.1], [1.135], [12.1%: outside the spread],
)
#figure(image("native-time.svg", width: 90%), caption: [Native ns per walk step, median and range (`runs/time/`).])

*IPC:* the host has no hardware PMU. `perf_ipc.sh` recorded `unsupported`, so no IPC is claimed.

== Scope

- *What remains.* The step is about 37 field multiplications' worth against a floor near 6–7. Most of it is four basis conversions, the rotation search, the `solve` body and the walk-start inversions. Going further would change the walk definition (pinned by committed replays) or need cache-cost measurement this host cannot make.
- *Counted units.* GAE records do not change. Future wall-time and instruction comparisons against this rho get harder for index calculus.
