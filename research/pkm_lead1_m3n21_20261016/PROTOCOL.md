# PKM lead 1 at `m = 3`, `N = 21`: protocol

**Status: frozen** by the commit that adds this line, after the identity run
of §7 and before any `N = 21` system has run. Card:
[`backlog/XP-20261015-pkm-m3-n21-m4-n20`](../../backlog/XP-20261015-pkm-m3-n21-m4-n20/README.md),
cell `m3-n21`. Encyclopedia §8 lead 1. Background:
[`RESEARCH_PKM_TOWER_ORACLE.md`](../notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md)
§11.6 (rule A1) and §16 (round 7, `N = 18`). This protocol does not cover the
card's `m = 4, N = 20` cell.

Every statement carries one label: *recorded* (a committed row), *derived*
(a calculation shown here) or *predicted* (stated before the runs).

## 1. Question

Does the Kummer tower's F4 solving degree `D` stay at 7 at `m = 3`,
`N = 21` (`t = 7`, `|V| = 128`), or does it rise?

The recorded line is `D = 6, 6, 7, 7` at `N = 9, 12, 15, 18` at `p₁`, with
two random targets per size and identical rows for both targets at every size
(*recorded*: round 6's `replay/M3-kummer-m3-p1.jsonl` and round 7's
`runs/M3b-t*.jsonl`).

## 2. The cell

| | value |
|:--|:--|
| engine | `f4_fp_tower` (`--engine tower`), full-rank exit on (default) |
| prime | `p₁ = 2013265921 = 15·2^27 + 1` |
| tower | Kummer, `m = 3`, `t = 7`, `N = 21` |
| targets | the cell's two random targets, `--random 2 --planted 0`, each run as its own process (`--targets 0`, `--targets 1`) |
| budget | `--budget 345600` (96 h) per process |
| memory caps | `--max-dense C` and `ulimit -v A`, from the host (below) |
| threads | all of the host's cores unless `RAYON_NUM_THREADS` is set; recorded |
| binary | `pkm_tower_pilot` built in release mode from the commit that adds this file |
| commands | [`run.sh`](run.sh): `host`, `identity`, `main 0`, `main 1`, `confirm T D`, `verify` |

**Memory caps.** These are fixed now, as round 7 fixed its rule for `C`
(note §16.3).
- `A = 0.92 · MemTotal / SHARE`, in bytes. `SHARE` is the number of
  processes running at once on the host. `ulimit -v` takes `A` in KiB.
- `C` is the largest multiple of `10⁸` with `4 · 1.45 · C ≤ A − 4·10⁹`.
  - *Derived*: `4` is bytes per `B'` entry. `1.45` is round 7's ratio of peak
    resident memory to the largest `B'` in bytes:
    `12,961 MiB / (4 · 2.25·10⁹ B)`. `4·10⁹` bytes is headroom.
- A run that the machine ends with an allocation failure, rather than the
  engine, is recorded as such, with the lower bound its trace gives. It is not
  re-run (round 7 §16.3).

**Host.** The host must have at least 192 GB of RAM for one process
(`SHARE = 1`), or 384 GB for both targets at once (`SHARE = 2`).
- *Derived*: the card's projection is 133–146 GB peak for one system.
- `run.sh host` records the host, its CPU, `MemTotal`, `A`, `C`, the commit,
  the binary's sha256 and `rustc` before any measured run.

## 3. Stages

1. **Identity.** `run.sh identity` runs the recorded `N = 15` system, target 0,
   on the measuring host and binary.
   - It must match round 6's committed row (`replay/M3-kummer-m3-p1.jsonl`,
     `N = 15`, target 0) on `solving_degree_max`, `max_cols`, `steps` and
     `inconsistent`.
   - Its `muladds` may only fall (note §15.1). The recorded value is
     3,185,681,037,100.
   - The wall clock, thread count and fields added since round 6 are exempt.
   - A failure stops the round until it is explained.
2. **Measurement.** `run.sh main 0` and `run.sh main 1`, each once. They run
   one after the other with `SHARE = 1`, or at once with `SHARE = 2`.
3. **Confirmation.** A target that finishes with `D` is re-run with the
   degree bound at `D − 1` (`run.sh confirm T D`). It must end with pairs
   above the bound, without refuting and without a staircase stop. That is
   round 7's confirmation (note §16.3).
4. **Verification.** `run.sh verify` runs the native checker
   `pkm_tower_check verify` over every row.
   - It recounts each finished tower row by exhaustive search over `V^m`.
   - *Derived*: at `t = 7` that is `128³ = 2,097,152` triples per row, which
     is cheap. The card's note that exhaustive search is infeasible at this
     size is wrong for `m = 3`.

## 4. Predictions (before the runs)

1. **Identity holds** (stage 1).
2. **`D = 7` on both targets, and both refute.** Confidence is low. The
   `m = 3` line rose once (`N = 12 → 15`) and then held (`15 → 18`). `m = 4`
   held at 7 for `t = 3, 4`.
3. **Both targets give the same `D`, width and step count**, as every
   recorded cell has.
4. **Size** (*derived* in PR #1893, card table; not run):
   - width ≈ 360,000–475,000 columns;
   - 2.5·10¹⁵–7.4·10¹⁵ multiply-adds;
   - 133–146 GB peak.

   These are too low if `D` rises to 8.
5. **Time** (*derived*): round 7 did 8.99·10¹³ multiply-adds in 176.4 min on
   four threads, about 8.5·10⁹ per second. At that rate one system takes
   82–242 h on four threads. With 32 threads at an assumed four- to six-fold
   speed-up, it takes 14–60 h. The speed-up is not measured.

## 5. Decision rule, fixed now

**A1** (note §11.6), per cell, over `N = 9, 12, 15, 18, 21`. The upper half
is `N = 18–21`.

| outcome at `N = 21` | final plateau `L` | reading |
|:--|--:|:--|
| `D = 7` on both targets | 6 (`N = 15 → 21`) | inconclusive: one more level without a rise. H1a needs `L > 10`, that is `D = 7` through `N ≥ 27` |
| `D ≥ 8` on either target | 0 | H0 for this cell (a rise within 4 in `N` over the upper half) |
| stopped with `D ≥ 7` only | — | inconclusive; nothing new |
| stopped after a productive step at degree ≥ 8 | 0 | H0, from a lower bound |
| allocation failure or timeout without such a step | — | recorded; not evidence about `D` |

**Instrument first.** These void the reading until they are explained:
- a failed identity check;
- a verdict that `pkm_tower_check verify` contradicts;
- a confirmation that does not behave as stage 3 requires.

**Inadmissible:**
- changing the engine, the cell, the budget, the seeds or the rule for `A`
  and `C` after stage 2 starts;
- dropping a target;
- re-running a system in the hope of a different answer.

A container or host loss is not an allocation failure. The run is repeated
whole with the same flags and binary, as rounds 3 and 7 did, and the lost
attempt's files are kept.

## 6. What was run before this protocol (disclosure)

- **Read, not run:** the committed rows and traces of rounds 6 and 7, and the
  card.
- **Run with the freeze commit's binary:** the identity stage on the
  authoring container (4 vCPU, 15 GB). §7 records its result. Nothing at
  `N = 21` has run on any build.

## 7. Pre-freeze identity run on the authoring container

Run on the authoring container (4 vCPU, 15 GB; `runs/host-vm.txt`), commit
`7513735f1`, binary sha256 `0aa2184a…8570f`, 2026-10-16 20:38–20:44 UTC,
`A = 15,161,070 KiB`, `C = 1.9·10⁹`. Row: `runs/identity-N15-t0-vm.jsonl`; trace:
`runs/identity-N15-t0-vm.log`.

**Result: the identity check passes.** Against round 6's committed row
(`replay/M3-kummer-m3-p1.jsonl`, `N = 15`, target 0), every compared field is
equal:

| field | new | recorded |
|:--|--:|--:|
| `solving_degree_max` | 7 | 7 |
| `max_cols` | 31,397 | 31,397 |
| `steps` | 24 | 24 |
| `inconsistent` (refuted) | true | true |
| `muladds` | 3,185,681,037,100 | 3,185,681,037,100 |
| `max_rows` | 52,129 | 52,129 |
| `max_residual_rows` | 30,474 | 30,474 |

The wall clock differs, as §3 allows: 371.8 s against 363.0 s. The native
checker agrees with the row (1 tower row checked against exhaustive search, 1
agrees). The protocol is frozen here. The measuring host repeats the identity
run (`run.sh identity`) before `main`.
