#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 7.5pt)
= Run report: GHS cover census n = 3 (larger q) and n = 5, run 20261014T0050Z

Descriptive record of one run. Interpretation is in the study `README.md` / `REPORT.pdf`. Conductor T-2467.

== What ran
- Spec: `PROTOCOL.md`, frozen at commit `2d95bde3f` before any panel; amendment `338cd011a` (hypothesis H) after panel D and before panels A, B, C, E produced output. Protocol-versioned at commit, not independently preregistered.
- Tool: `cover_census_n.rs` (dependency-free Rust, `rustc --edition 2021 -O`), sha256 in the environment block. Orchestration `run.sh` (pilot, reproduction arms, panel D) and `run_panels.sh` (panels A, B, C, E); comparisons `compare.sh` (`jq` + `cmp`).
- Isolation: none. Shared Apple M4 Pro host, load average 60--260 throughout; 6 threads. Wall times are bookkeeping, not results. All recorded values are exact, deterministic counts.
- Seeds: fixed-seed stream of the merged tool (seed 20261012; tag q for n = 3, tag 100q + n for n = 5).
- Interruptions: the session that wrote the protocol ended during a timing pilot (`n = 5, q = 5, sample:1`, output unread); its scratch outputs, including the first `n = 3, q = 3` pilot records, were lost, and the pilot was rerun unchanged in this run. A later session interruption did not stop the running panel process. Panel E was killed once by a host restart about 2h15m into its first attempt, before it wrote any curve record; it was rerun from scratch with the same command and seed (recorded in the environment block). No other arm was restarted.
- Raw records: `results/run-20261014T0050Z/*.jsonl` (`rec` = field, curve, group, mr, derived_checks, totals), stderr files (all empty), `reproduction-checks.txt`.

#raw(read("../../results/run-20261014T0050Z.env.txt"), block: true)

== Arms requested and completed
#table(columns: 5, stroke: 0.4pt,
  [*Arm*], [*Requested*], [*Completed*], [*Structures per model*], [*Notes*],
  [pilot], [n=3, q=3, both methods], [56 classes each], [819], [pilot, char 3],
  [R1], [n=3, q=5, 7, all classes], [252 + 690], [16275 / 120099], [compared with merged run 20261013T0403Z],
  [R2], [n=3, q=11, sample 120], [120 + 127 neighbours], [1786323], [compared with merged sample],
  [D], [n=5, q=3, all classes], [488], [597861], [panel],
  [R3], [n=5, q=3, cocycle method], [488], [597861], [compared with D],
  [A], [n=3, q=11, all classes], [2662], [1786323], [panel],
  [B], [n=3, q=13, all classes], [4400], [4855539], [panel],
  [C], [n=3, q=17, sample 200], [200 + 187 neighbours], [24221379], [panel],
  [E], [n=5, q=5, sample 60], [60 + 64 neighbours (after one restart: host restarted 2h15m into the first attempt, no records written; rerun from scratch, same seed)], [254313125], [panel],
)

== Reproduction checks
#raw(read("../../results/run-20261014T0050Z/reproduction-checks.txt"), block: true)

== Recorded values
Unit: number of isomorphism classes of curves (exhaustive arms) or sampled curves. `mismatches` (genus formula against the sum over quadratic subcovers) is 0 in every arm; `derived_checks` reports 0 violations in every arm. The plotted heights are recorded `group` counts divided by the recorded group size.

#figure(image("n3_min_genus_by_mod4.svg", width: 100%), caption: [Run 20261014T0050Z, n = 3. Minimal genus over the isogeny class (q = 5, 7, 11, 13: every class) or within one 2-isogeny step (q = 17: sample of 200), by \#E(K) mod 4; n under each bar is the recorded curve count.])

#figure(image("n5_min_genus_by_2torsion.svg", width: 100%), caption: [Run 20261014T0050Z, n = 5. Minimal genus by (\#E(K) mod 4, \#E(K)[2]); q = 3 every class, q = 5 sample of 60.])

#figure(image("n5_m_r_q3.svg", width: 100%), caption: [Run 20261014T0050Z, n = 5, q = 3: every structure of every curve by (r, m) and \#E(K)[2], from the `mr` records.])

#figure(image("n5_m_r_q5.svg", width: 100%), caption: [Run 20261014T0050Z, n = 5, q = 5: every structure of the 124 sampled and neighbour models by (r, m) and \#E(K)[2].])

== Missing data
- Isogeny-class minima for the sampled arms C and E.
- Hyperellipticity of the n = 3, m = 2 genus-3 covers (not decided by branch data).
- No timing result of any kind.
