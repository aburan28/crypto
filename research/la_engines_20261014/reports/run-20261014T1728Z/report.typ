#set document(title: "Run report: IC linear-algebra engine grid")
#set page(paper: "a4", margin: 2cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8.5pt)
= Run report: `research/la_engines_20261014`, 2026-10-14T17:28Z – 2026-10-15T01:14Z

Descriptive record of every run; interpretation is in `../../REPORT.pdf`.

*Commands.* `scripts/grid.sh` (frozen protocol, `PROTOCOL.md`, commit `404251aee` frozen, instrument fixes `8a82cf286` before data), `scripts/grid_amend1.sh` (`abdf48826`), `scripts/grid_amend2.sh` (`5731279cf`), `scripts/rerun_contended.sh` (amendment 3), `scripts/grid_amend4.sh` (`dca78f544`), `scripts/grid_amend5.sh` (`76f86d667`), `scripts/e2e.sh` (`ic logs`, configurations `n23-legacy:23:2:0:legacy`, `n31-search:31:3:0:search`, `n37-search:37:3:0:search`). CADO-NFS `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b`, Release, GF(p) layers `p1;p3`.

*Host.* 4 vCPU Intel Xeon @ 2.80 GHz (AVX-512), 15 GiB, Linux 6.18.44, rustc 1.97.0; shared cloud VM, no core reservation: isolation L0. Rayon pool 4 threads (1 in the one-thread pass); CADO `thr=2x2` (`1x1`).

*Workload.* Synthetic relation systems, m = 3 (m = 2, 4 at D = 16,384), D = 2#super[10]…2#super[18], seeds 1–3, modulus 2#super[62] − 57; coefficients ±1, or +1 for amendments 1, 2 and the `-pos` arms. Real: `ic logs` on Koblitz curves, seeds 1–3, three linear-algebra modes.

#table(columns: 6, align: (left, right, right, right, right, left),
  table.header([records], [runs], [failed], [verified], [timeouts], [file]),
  [frozen grid], [210], [4], [192], [3], [`raw/grid.jsonl`],
  [amendment 1], [72], [1], [70], [0], [`raw/grid_amend1.jsonl`],
  [amendment 2], [57], [0], [57], [0], [`raw/grid_amend2.jsonl`],
  [amendment 3 re-runs], [38], [1], [35], [0], [`raw/rerun_contended.jsonl`],
  [amendment 4], [243], [0], [242], [0], [`raw/grid_amend4.jsonl`],
  [amendment 5], [135], [0], [134], [0], [`raw/grid_amend5.jsonl`],
  [`ic logs` end to end], [27], [0], [27], [0], [`raw/e2e.jsonl`],
)

*Failures and unverified runs, with their recorded reasons.* Timeouts (1,800 s): `core-bw`, `core-bw` A/A and `cado-bwc` at D = 262,144 seed 1 (the stop rule then skipped seeds 2–3). `cado-bwc` m = 2 seed 1: exit 1 on a 25-column core. Amendment 1 `cado-full` D = 32,768 seed 3 and its re-run: `bwc` balancing abort ("expects a column permutation replicated on rows"); the first amendment-1 start used a stale binary (`raw/grid_amend1.falsestart.*`). Unverified: frozen-grid `core-bw`/`cado-bwc` on three singular unmerged cores (D = 8,192 s3, 16,384 s3, 32,768 s1; CADO `rows_ok` true); CADO-NFS pipeline D = 1,024 seed 2 in every session.

*Measured values.* Medians and min–max per (arm, D) for every pass are in `../../summary.json` (computed by `scripts/plots.py` from these records, no other arithmetic); per-run values, CADO per-step wall times (`steps_wct`) and filter reports are in the records.

#figure(image("../../figures/amend5_pipeline.svg", width: 100%), caption: [Amendment-5 session (records `raw/grid_amend5.jsonl`), median of 3 seeds, bars min–max, `solve_s`.])

*Controls.* Every arm is checked against the planted logarithm vector. A/A: `core-bw` repeated on seed 1 at every D (frozen grid). Contention: 38 runs overlapping builds or diagnostics re-run on the idle host (amendment 3).
