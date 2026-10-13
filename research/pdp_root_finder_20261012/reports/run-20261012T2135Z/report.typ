#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)
#show raw: set text(size: 8.5pt)
#show figure.caption: it => it.body

= Run report: the isolated timing stage of `pdp_rootfind`, run-20261012T2135Z

*Observations only.* Interpretation: `../../README.md` and `../../REPORT.pdf`.

== What ran

- *Command:* `research/pdp_root_finder_20261012/run_timing.sh 5` from the repository root, started 2026-10-12T21:35:53Z, finished 2026-10-12T22:25:23Z. Every timed command ran through `tools/isolated_bench.py run --cpus 3 --label <label>` (CPU 3 reserved, every other thread moved off it, a 2 s settle check before each run, PSI and other-process CPU recorded; a run whose record says `contended` is listed and not pooled), with `RAYON_NUM_THREADS=1`.
- *Binary:* `target/release/examples/pdp_rootfind` (release build of `examples/pdp_rootfind/` as committed in this PR, base commit `cd0286878`); the `analyze` subcommand of the same source tree is not part of any timed arm.
- *Host:* Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU (no SMT siblings reported), L2 2 MiB a core, L3 260 MiB, 15 GiB RAM, kernel 6.18.44, Python 3.13 (the isolation tool). Isolation level: a reserved core on a shared virtual machine (the tool's records in `../../raw/isolated.jsonl`).
- *Rounds:* five, arms interleaved within a round (A B C D per rung, rung after rung), an A/A pair of the baseline at $ell = 10$ opening every round.
- *Binary workload* (`bin-time`): arms `baseline` (PS-0's `decompose`), `group-scalar` (PS-2's row kernel), `rf-pairs`, `rf-mitm`; rungs $ell = 8 dots 14$ ($n = 3 ell$) on PS-0's corpora ($ell <= 12$) and this study's ($ell = 13, 14$), and the (23, 8) cell; the first 16 refuted uniform targets of a rung, 3 passes ($ell <= 11$), 1 pass ($ell = 12$), 4 targets and 1 pass ($ell >= 13$). The unit is ns per target (total wall of all calls ÷ calls); the table arm's table build is recorded apart as `setup_secs` and is not in the per-target figure.
- *Prime workload* (`prime-ladder`): Kummer towers at $p_1 = 2013265921$, seed 20261012, 16 uniform + 8 planted targets a rung; arms `pairs3`, `mitm3` ($m = 3$, $t = 4 dots 13$), `subtract2` ($m = 2$, $t = 4 dots 12$), `mitm4` ($m = 4$, $t = 4 dots 11$); the first 8 uniform targets timed in full-scan mode (every candidate probed, the cost of a refutation), passes repeated until 0.3 s. Unit: ns per full scan.
- *Requested and completed:* 185 commands requested (170 binary `bin-time`, 15 `prime-ladder`); 177 ran, all with exit status 0; 8 were refused before starting by the isolation tool's settle check (other processes used 0.25–0.38 CPU s in the 2 s window, limit 0.20; the top process in every refusal was the session's own agent process `claude`): `baseline-l8` (round 5); `group-scalar-l8` (rounds 2, 4, 5); `rf-pairs-l8` (round 5); `baseline-l9` (round 3); `group-scalar-l9` (round 3); `group-scalar-l12` (round 1). Of the 177 runs the tool flagged 7 as contended (other-process CPU above the limit during the run; the same agent process), listed and not pooled: `round2-aa2-baseline-l10` (0.31 s of other CPU in 2.64 s); `round2-baseline-l8` (0.14 s of other CPU in 0.45 s); `round2-group-scalar-l9` (0.01 s of other CPU in 0.04 s); `round2-rf-mitm-l12` (0.11 s of other CPU in 0.33 s); `round4-baseline-l8` (0.36 s of other CPU in 0.45 s); `round4-rf-pairs-l9` (0.01 s of other CPU in 0.05 s); `round4-rf-mitm-l10` (0.06 s of other CPU in 0.02 s). Every prime-ladder command ran in every round.

== Measured values

#figure(image("rounds.svg", width: 100%), caption: [Every timed row of run-20261012T2135Z: µs per target (binary) and per full scan (prime) on the reserved core, one dot per round, medians joined; hollow rings mark runs the isolation tool flagged as contended.])

Binary, ns per refuted target (median and minimum over the uncontended rounds; the counted operations of the same targets beside them):

#include "table_binary.typ"

Prime, ns per full scan (median and minimum over the uncontended rounds):

#include "table_prime.typ"

A/A control (the baseline against itself at $ell = 10$, every round): 171.0 / 159.1, 156.7 / 155.1, 172.9 / 176.6, 170.1 / 161.1, 174.7 / 169.7 ms per target (rounds 1–5; the round-2 second run is one of the contended runs above); within a pair at most 7.5 % apart, 155–177 ms across the ten runs.

== Records

`../../raw/bin_time.jsonl` (one row per timed binary command), `../../raw/prime_time.jsonl` (one row per rung and arm), `../../raw/isolated.jsonl` (the isolation tool's record of every run: label, command, reserved CPUs, settle check, PSI before and after, other processes' CPU seconds, the `contended` flag, exit status), `../../raw/timing_progress.txt`. Graph source: `plot.py` (reads the raw records only). Medians, minima and ratios: `../../ANALYSIS.json` and `../../ANALYSIS.md`, written by `pdp_rootfind analyze`.
