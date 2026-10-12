#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.5pt)

= Run report: lit2026 W11, run-20261011T2049Z

*What ran.* Measurements for workstream W11 of the lit2026 speedup program
(branch `lit2026/w11-generic-engineering`, stacked on aburan28/crypto\#1681).
Observations only; interpretation is in the study README / `results.pdf`.

- *Host.* Apple M4 Pro, 14 cores, macOS 26.6, shared with other sessions; 1-minute
  load 72–155 (recorded with every BSGS and class-walk sample). *Isolation L0.*
- *Arms.* `original` = program baseline `4b6e09b78` (baseline worktree, read-only);
  `previous`/`candidate` = one binary of this branch with the switch off/on.
- *Raw records.* `raw/*.jsonl` (one JSON object per sample), `raw/classwalk_19_1/`
  (samples.tsv, round-1 stdout and snapshots of each arm), `raw/summary.json`
  (medians, min–max, per-round ratios, identity checks; computed by
  `examples/lit2026_w11_summarize.rs`). Preliminary runs: `raw/prelim/`.
- *Requested vs completed.* All requested runs completed. One preliminary biext run
  stopped at `r_bits = 56` (no prime `p = c·r − 1 < 2^62` with `c ≤ 64`); the final run
  covers 16–52 bits.

== (a) BSGS step kernels

Counted field operations per visited point (`raw/bsgs_opcount.jsonl`, W = 8, giant
phase shown in the graph's right panel); every kernel visited the same `(index, x)` pairs
as the chain (`same_points_as_chain = true` in every row).

#image("bsgs_opcount.svg", width: 100%)

Wall seconds per solve batch, 5 rounds × 4 arms × 5 curves (`raw/bsgs_time.jsonl`);
answers digest and table entries identical across arms at every size.

#image("bsgs_wall.svg", width: 85%)

#table(columns: 6,
  [curve], [original s], [previous s], [paired s], [coz s], [original/paired (min–max rounds)],
  [secp40k1], [4.67], [4.28], [3.20], [3.49], [1.459 (1.054–1.917)],
  [secp44k1], [5.63], [6.12], [4.27], [4.22], [1.316 (1.015–2.002)],
  [secp48k1], [8.63], [9.55], [6.64], [8.16], [1.301 (1.299–1.603)],
  [secp52k1], [15.55], [19.40], [13.58], [14.77], [1.145 (0.996–1.976)],
  [secp56k1], [45.84], [48.99], [27.44], [35.32], [1.670 (1.349–1.670)],
)

== (b) x-only binary isogeny evaluation

`raw/binary_sect163.jsonl` (ℓ = 3, 5, 7; 8 steps each; 2000 reps × 7 rounds per
step), `raw/binary_torsion19.jsonl` (458 kernels, 20 reps × 9 rounds),
`raw/classwalk_19_1/` (12 rounds × 3 arms; snapshots identical without `seconds`).

#image("binary_xonly.svg", width: 100%)

Torsion-walk step: `F_(q^8)` mults 627 460 → 419 070, inversions 916 → 458, µs per
kernel 1564 → 1211 (medians; per-round ratios 0.94–1.93). Class walk (walk seconds,
medians): original 1.280, previous 1.164, candidate 1.198.

== (c) k = 2 pairing loops

`raw/biext_k2.jsonl`: values of the three loops identical at every size and equal to
`eds_mov::modified_tate` at (1283, 107); raw cubical value = square of the pairing on
every instance.

#image("biext_k2.svg", width: 75%)

== (d) DP statistics

`raw/dp_stats.jsonl` (K0/F2^53: 64 solves per dp; K0/F2^61: 32 per dp).

#image("dp_stats.svg", width: 100%)

== (e) Challenge ladder

`raw/ladder_derive.jsonl` (levels 30–60), `raw/ladder_solve.jsonl` (720 verified
solves). Pooled iterations/√(πn/4): mean 1.077, s.d. 0.546, s.e. 0.020, median 0.998.

#image("challenge_ladder.svg", width: 80%)
