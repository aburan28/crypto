#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set heading(numbering: "1.")

= Strong rho, round 3: the walk step (2026-10-16)

*Study:* `research/rho_step3_20261016/` · *Protocol:* `PROTOCOL.md` (declared before the changes; one addendum, recorded before any timing) · *Baseline:* `ed282c5b9` (round 2) · *Host:* one x86-64 cloud VM (Intel Xeon @ 2.10 GHz, 4 vCPU), no hardware PMU.

The strong rho is the reference every Koblitz `vs_rho` comparison prices index calculus against. This round makes its walk step cheaper without changing any walk. Every output line at $n$ = 41, 53 and 61 is identical to round 2's, timing aside, and all 14 strong-rho unit tests pass with debug assertions on.

== Changes

+ *Target-feature walk.* `solve` runs its body inside a `pclmulqdq` target-feature method when the field is a folding `Gf2` on x86-64, so field products inline into the step instead of being calls to `mul_fold`. The declared closure form did not inline and was replaced (addendum).
+ *Batched walk starts.* The starts $c + j s$ are formed `lanes` at a time with one batched inversion, instead of a full Itoh–Tsujii inversion (about 640 instructions) per walk. Same points, same order, same charged additions.

== Instructions (measured, Callgrind, $n$ = 53, 4 fixtures)

#table(columns: 3, align: (left, right, right),
  [*arm*], [*instructions*], [*per walk step*],
  [round 2 (`ed282c5b9`)], [2,369,297,848], [781],
  [change 1], [2,323,045,147], [766],
  [*changes 1+2*], [*2,178,454,252*], [*718*],
)
The final build runs 8.1% fewer instructions per step (1.09×) than round 2, and 23.5% fewer than round 1's baseline of 938. Output digests are identical in every arm.

#figure(image("instructions-per-step.svg", width: 90%), caption: [Instructions per walk step by arm (`runs/callgrind/`).])
#figure(image("instructions-by-function.svg", width: 90%), caption: [Where a step's instructions go, by arm.])

== Native time (measured)

A/A (round 2 × 6), then 5 interleaved A/B rounds, 32 processes under `isolated_bench`. Four were contended and are excluded: two A/A at n61 and two A/B round-2 runs at n53.

#table(columns: 5, align: (left, right, right, right, left),
  [*curve*], [*round 2 ns/step*], [*final ns/step*], [*ratio*], [*A/A spread, verdict*],
  [`icv1-f2m53-…`], [106.8], [105.9], [1.008], [5.6%: not established],
  [`icv1-f2m61-…`], [115.0], [112.1], [1.025], [2.4% (4 runs): marginal, not claimed],
)
#figure(image("native-time.svg", width: 90%), caption: [Native ns per walk step, median and range (`runs/time/`).])

The instruction cut does not show as a clear time gain on this host. The removed calls and inversions were cheaper in time than in instructions.

*IPC:* the host has no hardware PMU, so no IPC is claimed.

== Scope

- *What remains.* 718 instructions per step: about 600 in the `solve` body with canonicalization inlined (four basis conversions and the least-rotation search dominate), 71 in the walk's batched inversion, 25 in its single inverse. Going further would change the walk definition, which committed replays pin.
- *Counted units.* GAE records do not change. Future wall-time and instruction comparisons against this rho get harder for index calculus.
- *Not measured.* Hosts other than one x86-64 VM; the m = 83 wide instantiation, which gets change 2 only.
