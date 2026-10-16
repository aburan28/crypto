# Item 9 protocol: online IC wins at equal advice, and the free-relation census

Frozen 2026-10-14 before any measured data. The only runs made before this
file were instrument smoke tests on 8 throwaway targets at n=41 (labels
`smoke`, not reused), to confirm that both arms recover `[d]G = Q` and to read
memory. Those runs are not part of the result.

## Check first (what already exists)

- `docs/ic/BOUNDARY_TARGETS.md` (n=61 `vs_rho` row) lists "a Bernstein–Lange
  rho reference arm at n=61+" as the next target. It has not been run.
- `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:127-149` gives only a textbook
  calculation (BL online ≈ 490 steps at n=61 against ≈ 4·10⁴ IC probes) and says
  the control "has not yet been run".
- `research/prime_fast_bench_20261008/probeA.md` ran preprocessing rho on
  *prime* curves only, with no IC arm.
- `examples/koblitz_rho_batch_ks.rs` has a `KIC_RHO_PRECOMPUTE_WALKS` BL mode
  (32 jumps, old packed arithmetic). It was never run against an IC table at
  equal bytes. No result in `research/` or `benchmarks/` uses it that way.
- E-3 (`research/pdp_preregistrations_20261007/PROTOCOL-E3-index-self-collisions.md`)
  is PENDING with no instrument.
- `git log` shows no commit for preprocessing rho at equal advice, or for a
  self-collision census.

## Hypotheses

- **H1 (Lead A).** At equal advice bytes, Bernstein–Lange preprocessing rho
  solves one fresh target in less online time than compact-orbit IC, on every
  rung measured. The textbook prediction is that BL online costs
  `≈ c·√(N/T)` steps with `N = r/(2n)`, while IC online costs
  `≈ r/(c'·K²n²)` probes. Both ratios are expected to be far below 1.
- **H2 (Lead B).** The number of verified, non-degenerate free relations
  (self-collisions of canonical S3 roots) in the K-column pair-state table
  follows the uniform-collision law `M²/(2·r/(2n))`, where `M` is the number of
  root keys. That is about `K⁴n³/(4r)` for unordered states with both roots
  stored. The leads note predicts `n⁴K⁴/r` (≈ 11,000 at n=61, K=600); E-3
  predicts `K⁴n³/r` (≈ 180).

## Arms (same targets, same host, single thread, `nice`, `prlimit --as=1 GiB`)

- **A. IC (incumbent).** `examples/koblitz_orbit_dlp_fast_online.rs`
  `construct:<n>:0:<K>`, unchanged and with default env, as in the admitted
  n=61 row (pair symmetry off, no Bloom filter). All targets run in one
  process after setup. The per-target `online_ms` covers the query through the
  `[d]G = Q` check, and `probes` and `root_calls` are recorded.
- **B. BL preprocessing rho at ≤ equal bytes.** The new
  `examples/ic_next10_09_equal_advice.rs bl` arm uses signed-Frobenius classes
  (normal-basis least rotation, sign by y), 2,048 `G`-multiple jumps, and
  4-step look-ahead against fruitless 2-cycles. The online walk is from
  `Q + [c]G`. The precompute and online walks use the **identical map**. The
  table holds `T = bytes/16` entries (8-byte canonical x plus 8-byte log),
  sorted, with binary-search lookup. The DP density is `2^-t` with
  `t = round(log2 √(N/T))`. The online time covers the per-walk `[c]G`
  scalar multiplications, the walk, the lookups and the `[d]G = Q` check. The
  per-walk cap is 20·2^t steps, and the per-target budget is 400·2^t steps (a
  target that exhausts it counts as a failure).
  Prefix tables (the first `T/4` and `T/16` distinct DPs) come from the same
  precomputation. Because the online walk sequence is a deterministic function
  of the target and the seed, a superset table can never take more steps on
  the same target. So BL at less than equal bytes is a conservative
  (rho-unfavourable) stand-in for BL at equal bytes.
- **C. Plain strong rho (reference of the published wins).**
  `examples/koblitz_rho_batch_ks_strong_online.rs` rung 3, 32 lanes, DP bits
  4, one process per target, `KIC_RHO_TARGET_POINT` as in the n=61 panel. This
  arm reproduces, on this host, the direction of the published online win
  (IC < plain rho), which shows that the IC arm runs at its admitted level.
- **Positive control (BL instrument).** BL online steps against the table size
  `T` over the prefixes should fall with a log-log slope in [-0.75, -0.25]
  (law −1/2).
- **Negative/null control (map identity).** BL at n=41 with
  `BL_ONLINE_MAP_SEED=mismatch`: the table is built with one jump map and the
  online walk uses another. The prediction is that the success rate within the
  400·2^t budget is < 25% (the table is useless without the identical map,
  which is Probe A's warning).

### Advice bytes of the IC arm (what B is given)

`bytes_IC = 16·next_pow2(4K²n)` (root-table slots as allocated)
`+ 24·K²n` (states)
`+ 32·2Kn` (base points and labels)
`+ 8·Kn` (shifted x).

The base `by_x` hash map and the allocator overhead are left out, which
under-counts IC and favours it.

## Sizes (a = 0 Koblitz curves; the largest that fit 1 GiB per arm)

| rung | IC K | IC advice (formula) | BL bytes given | precompute |
|---|---:|---:|---|---|
| n=41 | 100 | ≈ 43.7 MB | = equal (prefixes /4, /16) | ≈ 1.3·10⁸ steps |
| n=53 | 200 | ≈ 320 MB | = equal (prefixes /4, /16) | ≈ 2·10⁹ steps |
| n=61 | 256 and 360 | ≈ 365 MB and ≈ 730 MB | **¼ of K=256** (≈ 91 MB; prefixes /4, /16) = 1/8 of K=360 | ≈ 2.8·10⁹ steps |

At n=61 the full equal-bytes BL table (≈ 5.5·10⁹ precompute steps) is
optional and secondary. It runs only if the 45-minute job cap and the shared
host allow it. The published n=61 point (K=500) and N83 K3240 (14 GiB) do not
fit in 1 GiB: K=500 at n=61 allocates a 1 GiB root table. They are not run.
The IC online cost at K=500 is bracketed by the trend over K ∈ {256, 360} and
by the published 64-target panel. This is labelled as projection, never as a
measured win.

Targets: 64 fresh public points per rung from
`ic_next10_09_equal_advice targets <n> 0 64 eval-v1`, the same file for every
arm. Scalars are kept apart for the verifier only.

## Metric

The primary metric is the per-target online wall time (ms), with both arms on
the same targets on the same host, in separate processes run back to back. It
is reported as the ratio `IC_online / BL_online`: the paired median of
per-target ratios and the ratio of means, with a 10,000-resample target
bootstrap 95% CI. The secondary metric is native counts (IC probes and root
calls, BL steps and walks) and ns per native operation.

## Kill criteria (from the report)

- **Lead A.** If preprocessing-rho online cost is below IC's online cost at
  equal bytes on every rung, the verdict is that the online-win rows
  (13.9x n=61, 54.1x N83) should be retired from `docs/ic/LEADERBOARD.md` as
  wins, and no new rungs should be landed. The coordinator edits the
  leaderboard; this item only recommends.
  Decision rule: POSITIVE for H1 (the audit kills the wins) if the upper end of
  the 95% CI of `BL_online/IC_online` is below 1 on every measured rung.
  NEGATIVE if any rung's lower CI end is above 1.
- **Lead B.** If fewer than 100 verified free relations appear at n=61
  (K=600), close the free-relation lever.
  The census is run at n=61, K ∈ {256, 360, 600, 800}, at n=53 K=200, and at
  n=41 K=100. It reports pairs, verified non-trivial, trivial, failed and
  rank mod r.

## Deviations allowed

Only the scale may change, and only downward under memory or time, recorded
in the README with the exact wall.
