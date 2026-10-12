#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)

= Run report: W8 interleaved session, run-20261011T2131Z

*Started* 2026-10-11 21:31 UTC. *Host* shared Apple M4 Pro (macOS), load average 100–341 (per-process
`uptime` in `raw/driver.log`). *Isolation* L0. This report describes what ran; interpretation is in the
study README.

== Arms and binaries

- *original*: `orig-bin` (sha256 `1a7b78ed…`): the W8 harness source compiled as a scratch crate
  against the program-baseline library, path dependency on the shared worktree
  `wt/lit2026-baseline` at `4b6e09b78` (nothing there edited); run with `--arms sym`.
- *previous*: `sym` in the candidate binary `w8bin5` (sha256 `b42e8494…`, library identical to `4b6e09b78`).
- *candidate*: `eg` in `w8bin5`, interleaved in-process with *previous*.

Driver `speed.sh` (thin shell): 5 rounds, ABAB with the process order flipped each round
(`orig cand`, `cand orig`, ...), each process one round over 8 targets per cell, `--no-degrees`,
`KIC_F4_MAX_DEGREE_ONLY=0`, seed `0x5EED`. Cells: K1/F2^15, K0/F2^15, K1/F2^17, K0/F2^23, K1/F2^23,
all `m = 2` (curve ids and subgroups as in the study README). Ratios computed by the harness:
`lit2026_w8_summation_models compare --inputs orig=raw/r1-orig.json,…,cand=raw/r5-cand.json --out compare.json`.

== Recorded values

#figure(image("eg_speed_session.svg", width: 100%), caption: [`compare.json`: per-round ratio of total ms
over 8 targets (dots) and ratio of medians (bar), and the deterministic word-op ratio, for original/candidate
and previous/candidate.])

#table(columns: 6, inset: 4pt,
  [*cell*], [*word ops orig = prev*], [*word ops cand*], [*ratio*], [*wall orig/cand median [min–max]*], [*wall prev/cand*],
  [K1/F2^15], [6 750 384], [4 218 414], [1.600], [3.73 [0.31–14.8]], [1.60 [0.59–8.62]],
  [K0/F2^15], [6 235 447], [3 797 654], [1.642], [0.96 [0.26–6.75]], [1.43 [0.80–4.45]],
  [K1/F2^17], [169 687 396], [117 819 615], [1.440], [3.11 [2.06–4.28]], [2.27 [1.46–2.86]],
  [K0/F2^23], [11 054 138 493], [8 205 050 470], [1.347], [1.62 [0.95–2.11]], [1.25 [1.21–1.26]],
  [K1/F2^23], [10 866 191 751], [8 107 487 123], [1.340], [1.84 [1.18–5.67]], [1.22 [1.12–1.40]],
)

Word-op totals were identical in all 5 rounds for every arm, and original = previous exactly.
Root-set digests (FNV-1a over sorted `w`-tuples, per target) identical across all three arms and all rounds.
Runs requested 10 processes, completed 10, exit 0.

== Observation (not part of the session)

`raw/observation-wide-path-K0-23.json`: `sym` and `eg` at K0/F2^23 with
`KIC_F4_ACTIVE_MULTIPLIERS=1 KIC_F4_SOLVER_LINEAR_TAIL=1 KIC_F4_FLAT_MATRIX=1` (3 rounds, `w8bin5`):
median per solve `sym` 2 976 ms / 1.02e9 word ops, `eg` 2 510 ms / 0.72e9; root digests equal to the session's.
