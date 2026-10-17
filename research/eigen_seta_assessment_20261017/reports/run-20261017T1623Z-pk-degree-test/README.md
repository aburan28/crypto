# Run report: run-20261017T1623Z-pk-degree-test

Observations only; interpretation is in [`../../REPORT.md`](../../REPORT.md).
PDF: [`report.pdf`](report.pdf).

* **Command.** `target/release/examples/pk_degree_test <SET> <wrong> <reps> <inspect>`
  run as `L1 100 7 1`, `L3 30 3 0` and `L5 15 3 0`, sequentially.
* **Code.** `eigenx` commit `94fe2d9` (branch
  `claude/eigenx-scheme-assessment-0z5d0u`, clean worktree), source
  `examples/pk_degree_test.rs` and `src/cryptanalysis/pk_test.rs`.
* **Host.** Intel Xeon @ 2.10 GHz, 4 vCPU shared cloud VM, kernel
  6.18.44, rustc 1.97.0, release profile. Isolation: none (L0, no pinning).
  During some L1 timing repetitions a Lean toolchain install and a Lean
  compile ran on the same host. Those reps are not identifiable after the
  fact, and the L1 ranges include them.
* **Arms.** Per set: one key from `keygen_with_trace` (OsRng), the true degree
  (both eigenvalue assignments), wrong candidates (uniform odd, prime to 21,
  and near misses f + j·2^k, k ∈ {16, 64, ⌊f_bits/2⌋}, j ≤ 4, with multiples
  of 3 or 7 skipped), and interleaved timing of one public-key test (theta
  chain) against one `EigenSeta::decaps` of an honest ciphertext.
* **Requested / completed.** L1: 100 uniform + near misses, 105 tested, 7
  timing reps, 1 true and 1 wrong structural inspection. L3: 35 tested, 3
  reps. L5: 18 tested, 3 reps. No failures, panics or timeouts.
  Every decapsulation returned the encapsulated key.
* **Values (measured, seconds, median [range]).**

| set | E | e (true f) | wrong tested / accepted | test | decaps | reps |
|---|---|---|---|---|---|---|
| L1 | 1076 | 569 | 105 / 0 | 0.79 [0.53–0.91] | 1.28 [1.11–1.57] | 7 |
| L3 | 1594 | 832 | 35 / 0 | 2.72 [2.64–2.78] | 5.54 [5.49–5.75] | 3 |
| L5 | 2106 | 1087 | 18 / 0 | 6.84 [6.66–7.03] | 13.34 [11.54–14.26] | 3 |

* **Structural inspection at L1 (measured).** True f: glue ok, no early split,
  563 Richelot steps, final split. Wrong f (uniform): glue ok, no early split,
  566 Richelot steps, no final split.
* **Records.** [`../../records/run-20261017T1623Z-pk-degree-test/`](../../records/run-20261017T1623Z-pk-degree-test/)
  (`L1.tsv`, `L3.tsv`, `L5.tsv`, `host.txt`). An earlier smoke run that was
  killed during its second inspection is kept as
  `../../records/smoke_l1_inspect_true_partial.txt` (true f, e = 566,
  final split).
* **Graph.** [`../../figures/pk_test.svg`](../../figures/pk_test.svg).
