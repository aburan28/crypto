#set page(width: 21cm, height: 29.7cm, margin: 1.8cm)
#set text(size: 10pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[Koblitz compact-orbit single-target panels: n=71 and n=83 hosted isolated reruns]
  #v(0.3em)
  #text(size: 10pt)[2026-10-14 • run ids `gha-38585913645-1` (n=71), `gha-38597929540-1` (n=83) • status `PENDING_INDEPENDENT_VALIDATION`]
]

= What ran

Independent-host reruns of the `koblitz.compact_orbit.n71_single_target` (n=71, a=0) and
`koblitz.compact_orbit.n83_single_target` (n=83, a=1) autolab beats on hosted GitHub runners
(AMD EPYC 7763, `ubuntu-24.04`), under `tools/isolated_bench.py reserve` with producers pinned
to the reserved CPU: 16 one-target eval workloads each, fixed K=600, strong-rho rung 3
comparison arm, one process per arm per workload, empty DP table per process. Targets are
public hash-to-curve points with no known scalar; verification is cross-arm scalar agreement
plus out-of-process Python certificate replay.

Both runs executed with the point-encoding normalization of PR #1827; without it the local
runs of the same beats finished `VERIFICATION_FAILURE` with 16/16 identical scalars (rho's
`published_q` decimal strings compared against Python-int targets).

= Measured values (hosted, isolated, 0 contended samples)

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  inset: 5pt,
  [Beat], [K], [W], [online speedup rho/IC (median, 95% CI)], [cold IC/rho], [IC / rho online median ms],
  [n=71 a=0], [600], [16], [*7.29x* (4.64--13.38)], [154.4], [518.6 / 3764.5],
  [n=83 a=1], [600], [16], [*12.29x* (3.26--32.21)], [133.7], [557.5 / 5442.8],
)

`verify` PASS (256 files, 0 mismatches each). Claim-checks 16/16 PASS per run. Local contended
Apple M4 Pro comparison runs read 5.6x (n=71) and 9.46x (n=83) — agreement within the
bootstrap intervals, but wall numbers from the contended host are secondary.

#figure(
  image("figures/n71-n83-hosted-online-ms.svg", width: 100%),
  caption: [Per-workload paired online walls (log scale) for both hosted runs, recorded values only.]
)

= Scope and non-claims

- Public synthetic hash-to-curve targets only; no key recovery, no asymptotic sub-rho claim.
- The IC's reusable setup (~500 s in process, mostly rank PDP) is outside the online interval;
  the cold IC/rho ratio is ~154x (n=71) and ~134x (n=83) with setup paid per process.
- One host class per arm of the comparison; no claim for other hardware classes.
- Out-of-process replay on the same runner; independent-author replay and a Bernstein--Lange
  reference remain outstanding. No ledger promotion.

Raw records: `autolab_runs/gha-38585913645-1/`, `autolab_runs/gha-38597929540-1/`,
`verify.json`, `isolation.jsonl`, `lscpu.txt` beside this report.
