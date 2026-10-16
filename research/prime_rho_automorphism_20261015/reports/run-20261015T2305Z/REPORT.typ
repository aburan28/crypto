#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set heading(numbering: none)

= Run report: automorphism rho on the round-2 crater curves

*Run start (UTC).* 2026-10-15T23:05Z. *Report type.* Descriptive run report
(observations only). Interpretation is in the study README.

== What ran

#table(
  columns: (auto, 1fr),
  [Protocol], [`research/prime_rho_automorphism_20261015/PROTOCOL.md`, frozen before the run],
  [Specs], [`specs/cubic.json` (`ECS1haae1b01f9883`), `specs/quartic.json` (`ECS1h26c7bf2d77bd`)],
  [Sessions], [`ECBS1h5e1b9912a46d` (cubic), `ECBS1hf8d1171784f1` (quartic)],
  [Comparisons], [`ECC1hb42c42de3b9f` (cubic), `ECC1hcc414017408e` (quartic)],
  [Source], [commit `c9db37d01bc776a1e8594ce5495e58ba2ccf5a1d`, `git_dirty: false`],
  [Binary], [`ecbench` SHA-256 `0a8a48d7…1473b6`],
  [Host], [AMD EPYC 9R05 VM, 2 logical cores, 1 NUMA node, Linux 6.18 x86_64; env class `ECBENV2h0e1f15f7c61e`],
  [Arms], [`rho.negation` (reference, A = 2) vs `rho.cubic` (A = 6) on six j = 0 curves; `rho.negation` vs `rho.quartic` (A = 4) on six j = 1728 curves],
  [Workloads], [the twelve crater curves of round-2 E2, 4 targets each, target seed 20261011],
  [Measurement], [3 alternating rounds, 1 warm-up, seed 1, 900 s timeout, L2 isolation requested],
)

== Completion

Each session planned and wrote 192 records, and all 192 are `verified`. There
were no failures and no timeouts. `ecbench verify --replay 12` reproduced 12 of
12 replays with 0 problems for each session.

*Isolation earned: L0 on all 144 paired attempts* (L2 was requested). The
harness recorded a p95 foreign busy share of 82.5% (cubic) and 83.0%
(quartic). Wall-time comparisons are therefore *descriptive, not results*,
because 0 pairs were admitted and at least 5 are needed. The harness's
measurement quality score is 0 / 100 for both.

== Measured values

The unit is `ecbench.gae` (S, group-addition equivalents). Values are ratios of
mean S, candidate over reference. Intervals are 95% cluster-bootstrap intervals
with 4,000 resamples, 72 same-seed pairs, same session.

#table(
  columns: (auto, auto, auto, auto),
  [*comparison*], [*ratio (measured)*], [*95% interval*], [*predicted (derived)*],
  [`rho.cubic` / `rho.negation`], [0.6112], [[0.505, 0.739]], [1/√3 = 0.577],
  [`rho.quartic` / `rho.negation`], [0.7910], [[0.673, 0.923]], [1/√2 = 0.707],
)

The predicted values are derived from the expected-steps model √(πr / 2A) with
A = 6 or 4 against A = 2, as frozen in the protocol. The harness labels both
ops ratios *"diagnostic only (unpriced work; ratio does not bound full
cost)"*. The cost of mapping each point to its class representative is not
priced in S.

Wall median ratios (descriptive only, not admitted): 0.980 (cubic) and 0.990
(quartic).

#figure(
  image("ratios.svg", width: 100%),
  caption: [Per-curve and pooled ratios of mean S with 95% cluster intervals.
    Each curve has 12 pairs. Comparisons ECC1hb42c42de3b9f and
    ECC1hcc414017408e. Plotted from the saved comparison JSON. No values were
    recomputed.],
)

== Raw records

- `research/prime_rho_automorphism_20261015/sessions/{cubic,quartic}/records.jsonl`
- `…/sessions/{cubic,quartic}/comparisons/*.json`
- `…/sessions/{cubic,quartic}.audit.json` and `.verify.log`
- `…/sessions/{cubic,quartic}.table.md`

== Missing data

- No L2-admitted wall-time interval exists. The host was busy with foreign load.
- There is no full-cost (priced) comparison. `ecbench price-units` covers
  Koblitz curves only.
