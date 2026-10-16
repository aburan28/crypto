# The scan keyed from slopes, pipelined: explorations on the Cascade Lake host

**What this is.** Timed explorations of the collection scan and the
pair-table build (plan §8's A2 and A1) on the host this container moved to
on 2026-10-07: a Cascade Lake, outside the programme's reference class.
They ran on 2026-10-08 from 02:46 to 04:44 UTC, before any round was
declared on them. They are records, not a round's measurement, and none of
their runs is pooled with any round's. They led to R09's declaration
([protocol](../../rounds/R09-slope-keys/PROTOCOL.md)), which discloses them.

**The outcome.** The final candidate runs **1.119 [1.093, 1.146], 1.112
[1.060, 1.166] and 1.103 [1.071, 1.136]** times faster in cold time than
v3 at `2^44.3`, `2^44.5` and `2^47.2` (eight clean pairs each), with every
logarithm v3's. Collection is 1.10–1.15× faster and the build 1.08–1.14×.
**R09 is declared on it**, for this class.

## The host and its class

`host.txt`: `Intel(R) Xeon(R) Processor @ 2.80GHz` (Cascade Lake), four
cores, 15 GB, kernel `6.18.44-fc-v80`, with AVX-512 F, DQ, CD, BW, VL and
VNNI and PCLMULQDQ, but no VPCLMULQDQ, GFNI, VBMI or VBMI2.
- **R06 and R08 cannot run here** (R06's amendment 1).
- **The scan's subtraction is the scalar product here** at every suite
  size: the eight-lane kernel needs VPCLMULQDQ.
- **Every process** ran under `isolated_bench` on CPU 2, `RAYON_NUM_THREADS=1`,
  between Track B's campaign's processes, under the same lock. A process
  its isolation record marks contended is left out of every figure.

## The scan's stages on this host

R04's probes on v3 (`ic-r04probes-on-995ea207-65cc4f6e`, the build the
[v3 stage record](../scan-stages-v3-20261006/README.md) used), `M1`'s two
rows at the three largest sizes, three rounds ([`probe.sh`](probe.sh),
[`stages.sh`](stages.sh), [`stages.tsv`](stages.tsv)). Nanoseconds a
scanned summand, the mean over clean processes:

| curve | `log₂ r` | processes | subtract | key | filter | admitted | scan |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 5 | 15.3 | 12.0 | 6.4 | 4.6 | 38.3 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 6 | 15.5 | 13.4 | 6.7 | 5.2 | 40.8 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 5 | 15.5 | 13.7 | 8.6 | 3.3 | 41.0 |

- **The subtraction leads here,** at 15.3–15.5 ns against 11.3–12.1 on
  the reference class (the scalar path there too at the two wide-tail
  sizes).
- **The scan is 57–59% of cold time** at the first two sizes and **85%**
  at `2^47.2`; the build is most of the rest.

## What was tried, in order

### 1. The scan pipelined across trials, alone

Candidate `7a1d132e`: a run's trial `t` prepared (rests and keys) before
trial `t − 1`'s filter is read, with `t − 1`'s filter words prefetched into
the second-level cache in step with `t`'s keys. One round of `M1`'s rows
([`explore-pipe.sh`](explore-pipe.sh)), stopped after it: cold time read
1.039, 0.963 and 1.017, and collection was slower at the two smaller
sizes.

The micro-benchmark below found why: the prefetch takes about 5 ns a
summand out of the probe pass and puts 2–4 ns into the preparation. The
words are fetched earlier, not more cheaply.

### 2. A2 re-checked here

A2's patch (`ic-a2-thp-on-995ea207-b2e9d484`, huge pages advised on the
pair tables), set aside on the reference class, re-checked because this
host's second-level TLB is smaller. Two rounds
([`explore-a2.sh`](explore-a2.sh); its candidate arm is named `pipe` in
the run tree): 0.950 at `2^44.3` (the first touch dearer) and 1.022
[1.002, 1.042] at `2^47.2`. Not the translation cost the probe pass was
suspected of; set aside again.

### 3. The micro-benchmark

[`examples/pipe_scan_bench.rs`](../../rounds/R09-slope-keys/candidate.patch)
times the scan over the same trials four ways in one process, alternating,
on a real folded table at `icv1-f2m53-tm56619371-dac20a85` (3,000 trials)
and `icv1-f2m61-t158598901-ab42b6c5` (1,500): the one-call scan (`call`);
the two halves with and without the prefetch (`pipe`, `bare`); and the
three-stage driver (`stages`). Each run was isolated; the minimum over its
rounds is the figure, because this host's noise only adds
([`microbench.tsv`](microbench.tsv)). Nanoseconds a scanned summand:

| binary | `2^44.3`: call | pipe | stages | `2^47.2`: call | pipe | stages |
|:--|--:|--:|--:|--:|--:|--:|
| the three stages (`stages-wip`) | 32.1 | 31.1 | 31.4 | 36.6 | 34.1 | 34.4 |
| plus the rotation step (`srlv-wip`) | 31.9 | 31.1 | 31.0 | 34.5 | 33.6 | 34.4 |
| plus keys from slopes (`slopes-wip`) | 33.0 | **27.8** | 28.4 | 35.2 | **30.0** | 30.6 |

- **Pipelining alone moves 3–7%** of the scan, and the three stages add
  nothing over two here.
- **The rotation step** helps the one-call scan at `2^47.2` (36.6 to 34.5)
  and is within the noise at `2^44.3`.
- **Keys from slopes** take the scan 15–16% below the one-call scan in the
  same binary.

These are stage diagnostics: they price the scan, not cold time.

### 4. The disassembly, which led to the slopes

`objdump` of v3's subtraction (`add_many_lazy_clmul`, `batch_inv_clmul`):
five folded products a summand, fifteen carry-less multiplies, each
seven cycles here and three of them chained in each product; the four
inversion chains cannot hide a 23-cycle link; and the slope loop squares
`λ` and stores a 24-byte rest for every summand, though the filter rejects
97–98% of them. In a normal basis `coords(λ²) = rotl(coords(λ))`, so a
rest's key can be formed from its slope with the basis change it already
pays, and the squaring, the sum and its store skipped.

### 5. The slope-keyed scan, without the build

Candidate `1e44a603`: one round, stopped when the build was keyed the same
way ([`explore-slopes.sh`](explore-slopes.sh)): 1.160 and 1.060 on two
pairs each.

### 6. The final candidate

`0b1afe38`: the three-stage scan keyed from slopes, the rotation step,
eight inversion chains, and the folded build's rows keyed from slopes. Four
rounds of `M1`'s rows at the three sizes, the order alternating, 48
processes ([`explore-final.sh`](explore-final.sh), [`pairs.tsv`](pairs.tsv),
[`intervals.txt`](intervals.txt), [`phases.tsv`](phases.tsv)). None was
contended or failed, and every pair recovered the same logarithm:

| curve | `log₂ r` | cold time, base over candidate | 95% interval | cold, ms | build, ms | collect, ms |
|:--|--:|--:|:--|:--|:--|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | **1.119** | [1.093, 1.146] | 931.7 → 837.2 | 273.9 → 254.0 | 584.6 → 509.2 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | **1.112** | [1.060, 1.166] | 956.8 → 875.6 | 295.7 → 260.4 | 574.7 → 520.7 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | **1.103** | [1.071, 1.136] | 3960.9 → 3614.8 | 578.2 → 505.6 | 3294.8 → 2980.5 |

(Medians over each arm's eight processes.) R09's candidate, `1d66b7dc` and
its formatting commit, adds only a gate that keeps the eight-lane kernel
where VPCLMULQDQ runs it; on this host the gate is always open, so the
code this exploration ran is the code R09 measures.

## Its limits

- **Eight pairs a size,** `M1`'s two rows only: an exploration, not a
  round. R09 runs 40 suite rows and 24 fresh holdouts a size, five
  rounds.
- **One host, one class.** No claim is made for the reference class, for
  CPUs without AVX-512, or for any other.
- **Shared with Track B's campaign** under one lock: every process
  isolated, none contended, but the explorations' processes interleave
  with the campaign's in time.
- **The micro-benchmark's spread is wide** on this host (the same binary's
  one-call scan read 36.6 ns in one process and 42.6 in another); only
  its minima are quoted, and only as stage diagnostics.

## Files

| file | what it holds |
|:--|:--|
| `probe.sh`, `stages.sh`, `stages.tsv` | the stage diagnostic and its table |
| `explore-pipe.sh`, `explore-a2.sh`, `explore-slopes.sh`, `explore-final.sh` | the four timed explorations, in order |
| `pairs2.sh`, `pairs-tsv.sh`, `phase-medians.sh` | A13's analysis scripts, unchanged |
| `pairs.tsv`, `phases.tsv`, `intervals.txt` | the final exploration's pairs, phase medians and intervals (the others' intervals in `intervals.txt` too) |
| `microbench.tsv` | every micro-benchmark round |
| `host.txt`, `binaries.sha256` | the host and every binary's SHA-256 |
| `runs.tar.xz` | every run tree (`scanp-cl`, `a2cl`, `pipex`, `pipex2`, `pipex3`, `pipeb`), SHA-256 in `runs.tar.xz.sha256` |
