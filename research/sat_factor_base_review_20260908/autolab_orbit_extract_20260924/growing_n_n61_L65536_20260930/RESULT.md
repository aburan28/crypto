# n=61, L=65,536: compact orbit vs historical KS v2 batched rho (2026-09-30)

Restart of the n=61 L=65,536 panel that PR #830 left aborted mid-corpus
(`../growing_n_n61_L65536.sh`). Status **`PENDING_INDEPENDENT_VALIDATION`**.
This is a **multi-target batch diagnostic against a superseded comparator**,
measured on a shared, unisolated macOS host. It does not promote any ledger row.

## What this is and is not

- Comparator: frozen Kuhn–Struik batched signed-Frobenius rho v2
  (`examples/koblitz_rho_batch_ks_v2_n61.rs`, `KIC_RHO_DP_BITS=4`). The
  2026-09-29 strong-rho ladder and sweep superseded it as the operative
  reference. Against strong rho R3, IC/rho is **2.297** at n=61 L=1,024 in
  retired instructions
  ([`RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md`](../../../notes/index-calculus/RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md)).
  Nothing below changes that verdict. This panel was not run against strong rho.
- One batch of L=65,536 targets per run. Under the AGENTS.md one-target rule
  this is a historical multi-target diagnostic, not the primary single-target
  comparison. The `vs_rho` claim-check fails closed, as it should (below).
- No key recovery. No asymptotic claim: at optimal K both arms scale as
  √(L·r/n), so any lead is a finite constant factor.

## Host and timing conditions

Apple M4 Pro, 14 logical CPUs, 48 GiB, macOS 26.6 arm64, rustc 1.93.1, single
thread, git `639886506` (`host.json`). **`tools/isolated_bench.py` does not run on
macOS** (it raises `AttributeError: os.sched_getaffinity`; it needs `/proc`,
CPU affinity and PSI). So the walls below come from `/usr/bin/time -l` with the
load average recorded per run. They are **not** AGENTS.md §10 evidence. Another
agent's job (`ecc2k83-m33`, ~21 GB RSS) shared the machine, and swap stayed
14–16.6 GB of 16–17 GB throughout.

Retired instructions (from `/usr/bin/time -l`) are the counted unit. On macOS
they include kernel instructions spent on page faults, so a run that swaps
overstates its count, as block 0's IC arm shows.

## Protocol and deviations

1. Corpora `n61-ks-growing-tune-65536-v1` (tune) and `n61-ks-growing-65536-v1`
   (eval), seed 531310, generated with `KIC_RHO_GENERATE_ONLY=1`. This uses
   the same blake3 fixture derivation as the walking run; every rho block's
   published scalars and Q match the corpus and the IC targets.
2. K tune on the disjoint 65,536-target corpus, lowest wall wins, candidates
   1,400 / 1,800 / 2,200. The frozen script listed 1,400–2,600; 2,600
   (~30 GB estimated) was dropped up front.
   - K=1,400: 2,478 s wall (339 s user, **1,061 s sys**: swapping).
   - K=1,800: **censored**. Killed at 4,908 s elapsed, more than twice
     K=1,400's wall, with 2.5 GB resident of ~18.6 GB estimated. Its wall is
     at least 4,908 s, so K=1,400 has the lowest wall either way.
   - K=2,200: **skipped**. Estimated RSS 27.75 GB exceeded the 23.7 GB
     available.
   - **K=1,400 is therefore a memory-capped choice, not a tuned optimum.**
     At L=16,384 the sweep already peaked at its top value, K=1,400.
3. Three paired eval blocks, alternating arm order (IC→rho, rho→IC, IC→rho).
   Each stage writes a `.done` marker, so an abort keeps the corpus and the
   finished blocks.

## Results (compact / rho, same 65,536 targets per block)

| Block | Order | IC wall s (user / sys) | Rho wall s (user) | Wall ratio | User-CPU ratio | Instructions ratio | Load at IC start |
|---:|---|---:|---:|---:|---:|---:|---|
| 0 | IC→rho | 1,798.2 (325.3 / **520.9**) | 913.7 (859.4) | **1.968** | 0.379 | 0.201 | 34.96 |
| 1 | rho→IC | 295.3 (274.0 / 17.0) | 859.8 (834.7) | 0.343 | 0.328 | 0.135 |  9.14 |
| 2 | IC→rho | 295.0 (273.9 / 16.7) | 839.1 (826.7) | 0.352 | 0.331 | 0.135 |  9.24 |

- **Wall ratio:** median 0.352, range 0.343–1.968. Compact is faster on wall in
  2 of 3 blocks. **Block 0 is a compact wall loss.** Its IC index build took
  1,259 s instead of the ~58 s it took in blocks 1–2, with 521 s of system
  time, because the 11.6 GB IC footprint was swapping under the foreign job.
  The loss is reported as measured and is not dropped.
- **User-CPU ratio:** median 0.331 (0.328–0.379).
- **Retired-instruction ratio:** median 0.135 (0.135–0.201). Rho's count is
  stable to 0.05% across blocks (16.861–16.870 T). IC is 2.2745/2.2756 T in
  the clean blocks and 3.39 T in block 0, which includes kernel paging work.
- Peak memory footprint: IC 11.55 GB, rho 2.47–3.27 GB. The IC arm uses
  3.5–4.7× rho's memory.
- Correctness: 3 × 65,536 targets solved in each arm, 0 failures. IC
  `recovered_matches_published` and `group_verified` hold for every target, as
  does rho's `recovered == published`. IC and rho walk the same Q in every
  block. Records with timing fields removed are byte-identical across the three
  blocks in each arm.
- Counts (deterministic): rho takes 394,269,676 walk steps, 420,718,431 group
  additions and 26,216,908 distinguished points, with 65,527 cross-target
  solves. IC has 119,560,000 regular index states (K²·n), a mean of 11,238
  rank probes per relation, and 11,179.7 target probes per target. These are
  different operations; no common operation unit (S) is defined, so no ratio
  of them is claimed.

### Trend against the L=16,384 panel (same K=1,400, same comparator)

| L | Wall ratio median | Instructions ratio (clean blocks) |
|---:|---:|---:|
| 16,384 | 0.205 | 0.109 (1.14 T / 10.52 T) |
| 65,536 | 0.352 | 0.135 (2.27 T / 16.86 T) |

At fixed K the compact lead over v2 **narrows** as L grows. IC instructions
double for 4× the targets, while rho's rise 1.6×, because batched rho's cost
per target falls with L and compact's does not at fixed K. Recovering the
L=16,384 gap would need a larger K, which this 48 GiB shared host cannot hold.

## Independent replay

`../independent_replay.py` (pure-Python GF(2^61) and Koblitz arithmetic, no
code shared with the producer) checks every IC record of every block against
`base_n61_K1400.jsonl`: `replay_b{0,1,2}/independent_replay.json`. **Result:
196,608 / 196,608 pass, 0 fail** (65,536 per block). Rho's logs are not
separately replayed.
Rho's scalars equal the corpus scalars, which the IC replay verifies as
[d]G = Q on the identical points.

## Claim-check

- `end_to_end_dlp`: `claim_end_to_end_dlp.json` →
  `claim_check_end_to_end_dlp.json`. **PASS** (all required fields present;
  this checks schema completeness, not whether the beat was won).
- `vs_rho`: `claim_vs_rho.json` → `claim_check_vs_rho.json`. **FAIL by design**
  (18 missing stage fields, `target_count must equal 1`).
  The schema requires one target, a single-target online wall interval, the IC
  phase split, a rho replay certificate and matched resource envelopes. This
  batch producer emits none of these, and no values were filled in.

## Files

Committed: scripts, `host.json`, `SHA256SUMS` (binaries), `SHA256SUMS.corpora`,
both scalar corpora, `panel.log`, `memory_checks.log`, per-run `.time`, `.load`,
`.avail`, IC `.summary.json`, rho summary lines (`ks_n61_b*.summary.json`),
block-0 raw records (`ic_n61_b0.jsonl.xz`, `ks_n61_b0.jsonl.xz`), the base
(`base_n61_K1400.jsonl.xz`), replay receipts, `growing_n_summary.json`, claim
drafts and receipts. Blocks 1–2 raw records (about 54 MB IC and 45 MB rho
each) are not committed. Their SHA-256 and byte counts are in
`raw_manifest.json`, and the summary's `untimed_record_sha256` shows they equal
block 0 apart from timing fields. Regenerate them with
`../growing_n_n61_L65536_20260930.sh` (`BLOCKS=3`) and summarise with
`../growing_n_n61_L65536_20260930_summary.py`.

## Left undone

- An isolated (§10) rerun on a Linux host through `tools/isolated_bench.py`, with
  A/A spread.
- K above 1,400 at L=65,536, which needs about 19–30 GB free for the IC arm.
- This L against strong rho R3, which decides the operative verdict.
- Single-target online instrumentation for the `vs_rho` schema.
