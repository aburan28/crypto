#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 10pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[Run report: pkm_supersingular_aux, 2026-10-15T23:15Z]
  #v(2pt)
  research/pkm_supersingular_aux_20261015 · observations only; interpretation is in `../../README.md`
]

*What ran.* `run_q3.sh` (ecbench sessions s20, s22, s24), `run_q2.sh` (F4 cells) and `run_q1.sh`
(structure), in that order, from commit `8d719fddf` (`PROTOCOL.md`), binaries in `binaries.sha256`
(`ecbench` `92e39467…`, `pkm_ss_tower` `f277e003…`), on a shared 4-vCPU, 15 GiB container: one thread,
`nice -n 10`, 4 GB address-space cap, ecbench at L0 (`--cpus none --allow-busy`). Start
2026-10-15T23:15:13Z; every start and end is in `progress.txt`.

*Failures.* A container restart killed `run_q2.sh` after 23:26:36Z, during `T31-ordinary-m34` (8
partial rows kept in `q2/interrupted/`). The script was made resumable (completed runs skipped) and
restarted; that run was repeated whole. `run_q1.sh` ran in resumable chunks of 300 in `p`. No system
timed out or stopped for size; no recount mismatched.

*Q3 (measured).* 3 sessions × 168 records (56 warm-up): 504/504 verified; `verify --replay-all`
reproduced 112 replays per session, audit OK ×3. Mean $S$ and paired ratios are in
`q3/tables/summary-arms.tsv` and `summary-ratios.tsv` (ecbench's ratios and 95 % cluster-bootstrap
intervals, 16 pairs on 8 workloads).

#figure(image("../../figures/q3_s_and_ratios.svg", width: 100%),
  caption: [Q3: mean $S$ per arm and paired ratios, sessions `q3/sessions/s20`, `s22`, `s24`.])

*Q2 (measured).* 336 systems in 24 run files (`q2/runs/*.jsonl`, `q2/cells.tsv`), two random targets
per cell, all inconsistent (no decomposition), 296 with an exhaustive recount that agreed.

#figure(image("../../figures/q2_d_vs_n.svg", width: 100%),
  caption: [Q2: F4 solving degree against $N$ at $m = 2$, per prime and tower family.])

*Q1 (exact).* 10 chunks, 103 primes $p equiv 3 mod 4$ in $[103, 2999]$, 195 $(p, t)$ cells, 12,581 rows (`q1/structure.jsonl`); brute-force agreement 906/906 at $p <= 131$; positive control 100 % in every cell. Chunk times 14 s to 23 min (the first carries the brute force).

#figure(image("../../figures/q1_equivalence.svg", width: 90%),
  caption: [Q1: fraction of sets Möbius-equivalent to a torus set (or, for the positive control, to a
  $y^2 = x^3 + x$ set), pooled over primes.])
