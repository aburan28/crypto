# R01: baseline v0

**Declared 2026-10-01, before any of its runs.** The plan is
`research/notes/index-calculus/IC_TOOL_PROGRAM.md`; this is its §12's
first round. Nothing below changes after the first run, except by a
dated amendment appended at the end.

## Question

What does v0 cost on suite v1, phase by phase? How noisy is this host's
timing? R01 has no candidate and claims no speedup. Its class is
**accounting**: it is a baseline.

## v0

- **Code.** v0 is `main` at `46ae2014`. Its `src/` tree is `003badc2`,
  the same tree the binary was built from (at `4afd2990`).
  - Its Koblitz code is identical to §23's binary (`0bf67f16`, `src/`
    tree `574b9a87`).
  - The only `src/` differences between the two are in
    `hyperelliptic_index_calculus.rs` and `hyperelliptic_ic_bench.rs`,
    which the Koblitz pipeline does not call.
- **Binary.** `ic` built with `cargo build --release`, rustc 1.94.1,
  sha256 `c776cc04a4e7f2657505f2ab192030ca7e3fb03fff5853709a4778106b726797`.
  It is kept outside the tree.
- **Provenance.** A report's own `commit` field is the working
  directory's `HEAD` (IC_TOOL_PROGRAM.md §9). `runs/host.json` records
  the build commit and the binary's hash beside it.

## Frozen inputs

Suite v1 (`research/ic_tool_program/suite/v1/`, frozen in this PR).
`make_suite.py --check` re-derives every file before every step:
- 88 S rows: 11 sizes × seeds 201–204 × 2 targets;
- 2 smoke rows (`E_0/GF(2^31)`).

`M1`'s 12 rows at §23's six sizes are byte-for-byte §23's `T01`/`T02`
parameter files.

## Steps, in order

1. **Host manifest** (`runs/host.json`). It includes THP's settings.
   This host's are `madvise` for both `enabled` and `defrag`.
2. **The pin, untimed.**
   - v0 runs on `M1`'s 12 rows at §23's six sizes, under
     `taskset -c 2`.
   - Their outputs must equal §23's `T01-R1` and `T02-R1` reports:
     - counts;
     - both arms' recovered scalars;
     - rho's counts;
     - the verification and agreement flags.
   - If the pin holds, §23 is v0's comparison under the single-target
     rule, at those six sizes. If it fails, R01 stops before step 4 and
     reports the rows that differ.
3. **The smoke tier, untimed.** Its two rows must be complete and
   verified.
4. **The profile pass.** v0 runs once over the whole S suite: 88
   isolated processes.
5. **The A/A.** v0 runs against a byte-identical copy (`A2`) on `M1`'s
   22 rows: five rounds, with the order alternating. That is 220
   isolated processes.
6. **A huge-page probe.** This is a diagnostic, not a candidate.
   - v0 runs as itself (`4k`) and with
     `GLIBC_TUNABLES=glibc.malloc.hugetlb=1` (`thp`). Under that
     setting, glibc's malloc asks for transparent huge pages on its
     large mappings. Nothing in `src/` sets an allocator or calls
     `madvise`.
   - Rows: `M1`'s rows at the three largest sizes (`k0n53`, `k1n59`,
     `k0n61`), for five rounds, with the order alternating. That is 60
     isolated processes.
7. **Memory calibration.**
   - The tool is the earlier study's `mem_calib.c`
     (`research/notes/index-calculus/cachegrind_n41_20260929_run/`),
     built with `gcc -O2`, with its hash recorded.
   - It measures random pointer-chase latency at 16 KiB to 1 GiB,
     plus memory-level parallelism at 256 MiB. It does so with and
     without the tunable, through `tools/isolated_bench.py`.
8. **Callgrind, untimed, after every timed step.**
   - v0 runs with `--repeats 1 --repeats-fast 1` on `M1-T01` at `k0n41`
     and `k0n61`.
   - The tool is `valgrind --tool=callgrind --cache-sim=yes
     --I1=32768,8,64 --D1=32768,8,64 --LL=2097152,16,64`, the earlier
     study's lower last-level size.
   - It records instructions, and D1 and LL misses, per function.

## Figures (`analyse.py`)

- **Per size:**
  - the set-up, online and rho times, and the set-up's and online
    interval's phase shares;
  - `unit_ns`, which is v0's pinned unit (IC_TOOL_PROGRAM.md §5);
  - nanoseconds and units per scanned summand, and per stored pair;
  - `S` cold, online and set-up, in v0's unit.
- **The A/A, per size:**
  - the geometric mean and 95% `t` interval of the ten paired
    `A`/`A2` cold-time ratios (two rows × five rounds);
  - the same for the online interval.
- **The huge-page probe, per size:** the paired `4k`/`thp` ratio of
  cold time and of the collection phase, with the user and system
  seconds from the isolation records.
- **The calibration:** latency by working set, and the huge-page
  effect.
- **Callgrind:**
  - the functions with the most instructions;
  - the scan's and the build's D1 and LL misses per scanned summand and
    per stored pair.
- **v0's row** in the programme's ledger (IC_TOOL_PROGRAM.md §7).

## Complete, stop, and what follows

- **Complete** when:
  - the pin holds;
  - every timed row has a clean, complete and verified figure;
  - the A/A is measured at every size.
- **Stop.**
  - A verification failure stops that size's remaining rows, as in §23.
  - A failed pin stops R01 before step 4.
- **What follows.**
  - R02's protocol takes its phase and lever from this profile.
  - The A/A band sets how many rounds R02 needs: a size whose A/A
    interval is wider than ±5% gets ten rounds instead of five.

## Accounting

- One thread, isolated.
- Cold time is the set-up plus the online interval, at the median
  in-process repetition (IC_TOOL_PROGRAM.md §5).
- A process's figure is its first clean, complete attempt. Contended
  and failed attempts are kept and counted, never pooled.

## Inadmissible

- Changing the suite or a row after a run.
- Re-running a row whose clean figure exists.
- Pooling contended runs with clean ones.
- Quoting the huge-page probe as a speedup. Its arm is a tunable, not a
  candidate binary; a candidate that asks for huge pages is a later
  round's.
- Quoting callgrind's cache simulation as measured hardware behaviour.
- Any claim about the method against rho beyond §23's.

## Cost

About 1.5 hours of isolated runs and 30 minutes of callgrind on this
4-core host. Nothing else runs on the machine meanwhile.
