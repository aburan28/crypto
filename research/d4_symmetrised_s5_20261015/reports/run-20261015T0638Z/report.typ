#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
= Run report: d4_weighted_s5 cells, 2026-10-15 06:38Z–07:43Z

*Command and build.* `raw/run_all.sh`, the cells of `PROTOCOL.md` §3, in order. The binary is `examples/d4_weighted_s5` built at commit `87fe78b63` (`raw/binary_commit.txt`).

*Host and isolation.* Claude Code cloud container, 4 vCPU Xeon, shared with four other lanes. Isolation level L0. Each run used one thread (`RAYON_NUM_THREADS=1`), `nice 10` and `prlimit --as` 2.5 GB.

*Runs requested and completed.*

#table(columns: 2,
  [requested], [completed],
  [controls: 2 invariance checks and 9 generic systems], [11 / 11],
  [tiny cells: $p in {13, 17, 29, 37}$, 6 targets each, arms w and c8, brute force and MITM], [24 / 24],
  [measured cells: $p in {269, 1033, 4093}$, 8 targets each, arms w and c8, MITM], [24 / 24],
  [w-only cells: same $p$, seed 2, 12 targets each], [36 / 36],
)

*Failures.* At $p = 17$, target 4, the quotient has dimension 514. Both arms report it unsolved ("border unreachable at this degree"), and brute force finds 4 decompositions there. No crash, timeout or out-of-memory event occurred.

*Values (measured).* Unit: $F_p$ multiplications per solve.
- Arm w: median $1.41$–$1.44 dot 10^9$ at $p = 269$, 1033 and 4093.
- C8: median $2.54$–$2.56 dot 10^10$.
- Paired ratio: medians 17.8–18.2, minimum 7.7.
- Lift-back: at most $2.8 dot 10^5$.

Full tables are in `summary.txt` and `README.md`. Raw records are `raw/{controls,tiny,measured,w_only}.jsonl` and their `.log` files.

#figure(image("../../figures/c_per_solve.svg", width: 100%), caption: [Every target's cost, both arms. Run window 06:38Z–07:43Z.])
