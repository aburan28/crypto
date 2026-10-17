# Matched point-count run 20261016T0205Z

Observed run, 2026-10-16 UTC. Original cf2745214d83e7bbf737f07d758e486d6af8fa7f
versus candidate 44f8d17f1c5b8988a57a16e9dddc8e5ce099d325; native measurement
source frozen at protocol commit bdcf3fca. Worker hashes and source identity
are retained in [frozen-hashes.txt](../../evidence/frozen-hashes.txt), with
Rust 1.99.1 (58cf8f9c3), release opt-level 3. The later CM generation fix is
outside these timed revisions. It has separate untimed regression evidence.

## Execution

Seven workloads, each five A/A pairs followed by five alternating A/B pairs:
140 requested records, 140 completed, exit 0, zero timeouts, zero mismatches,
zero contention exclusions. A is always the original binary. An individual
child has a 120-second limit; the suite has a 60-minute limit. Cache setup
was prepared before the timed sequence. The original preparation recorded
82.296 seconds setup, separately from the measured cached-stage calls.

The host is an AMD EPYC 9V74 VM, Linux x86_64, five visible CPUs. The native
runner holds /tmp/crypto-bench.lock, pins the child to CPU 4, evicts movable
visible threads, observes process CPU use, and restores affinities. PSI is
unavailable and retained as null. Unmovable system threads are recorded in
[host.json](../../evidence/run-20261016T0205Z/host.json). Cloud neighbors and
frequency are uncontrolled; no exclusive-host isolation is claimed.

## Inputs and included work

| Workload | Registered model | Calls per record |
| --- | --- | ---: |
| bsgs24 | icv1-fp24-t6316-2e757e7a | 1000 |
| bsgs32 | icv1-fp32-t6738-7a2b626e | 300 |
| bsgs40 | icv1-fp40-tm1043-a68449f5 | 100 |
| sea40_cached | icv1-fp40-tm1043-a68449f5 | 100 |
| sea192_cached | icv1-fp192-t31607402316713927207482677199-52e4af59 | 1 |
| count192_cold | icv1-fp192-t31607402316713927207482677199-52e4af59 | 1 |
| cm256_control | icv1-fp256-t432420386565659656852420866390673177327-76dadd18 | 10 |

Seeds are 1+i%4 within a record; one-call records use seed 1. Every result
asserts the frozen registered order; complete API calls require both order
certificates. The analyzer checks paired work, seeds, exact model and order,
then independently certifies every reference after timing. This independent
check shares the library arithmetic and probable-prime assumptions; it is
not an independently verified arithmetic implementation.

BSGS rows time complete curve::order calls. Cached rows time sea::sea after
canonical Phi-cache import and a warm-up; setup is separately recorded. Cold
API counting includes fresh Phi construction, counting and certification.
The unchanged CM count control uses the complete API on secp256k1. Validated
input model construction, process launch and output serialization are excluded.
Full accounting is in [MEASUREMENT.md](../../MEASUREMENT.md).

## Native summary

Times are measured microseconds per call. Ratios and envelopes are derived
by the native analyzer: median paired A/B ratio; observed paired minimum and
maximum; maximum absolute A/A ratio deviation from one. There are five A/B
pairs per row. An effect within its A/A envelope is inconclusive under the
frozen policy. The envelope is descriptive, not a confidence interval.

| Workload | A median us/call | B median us/call | A/B median ratio | Paired range | A/A max deviation | A/B pairs |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| bsgs24 | 61.747 | 31.495 | 1.98652 | 1.88259..2.14328 | 26.575% | 5 / 5 |
| bsgs32 | 327.251 | 117.994 | 2.79019 | 2.73356..3.14247 | 17.594% | 5 / 5 |
| bsgs40 | 1936.899 | 591.519 | 3.24310 | 3.07922..3.42764 | 8.710% | 5 / 5 |
| sea40_cached | 264.650 | 314.525 | 0.81228 | 0.80232..0.93396 | 49.308% | 5 / 5 |
| sea192_cached | 1904284.547 | 1825125.454 | 1.02393 | 0.90254..1.10098 | 6.766% | 5 / 5 |
| count192_cold | 98449066.339 | 97670607.289 | 1.00469 | 0.99858..1.01812 | 2.869% | 5 / 5 |
| cm256_control | 19048.878 | 17948.105 | 1.03405 | 0.98103..1.17123 | 19.957% | 5 / 5 |

The native comparison admits the three BSGS effects beyond their A/A envelope.
The other four effects are within it. Cached SEA40 has ratios below one in
all five pairs. No complete-API or SEA improvement is established by this run.

## Evidence and replay

[runs.jsonl](../../evidence/run-20261016T0205Z/runs.jsonl) retains every
observation, process accounting and output. [summary.json](../../evidence/run-20261016T0205Z/summary.json)
retains the unrounded native statistics, including minimum arm times and
setup medians. [Reference certificates](../../evidence/run-20261016T0205Z/reference-certificates.json)
retain the independent post-timing order checks. No failed samples were omitted.

With Rust 1.99.1 installed, run bash replay-benchmark.sh from the study
directory. It builds both frozen source worktrees, generates the excluded
canonical cache, runs native reserve and analyze, and retains a fresh output
directory. The original command used the baseline worker under
/workspace/sea-original and the candidate worker under this study's native
target/release, with output evidence/run-20261016T0205Z. Exact binary hashes
are retained; rebuilding with the later constructor fix is not the timed arm.

The [results report](../../REPORT.md) includes the formal checks and CM bug
correction. Positive TLC/Lean results do not extend the measured workloads
or prove the concrete Rust implementation. render_reports.py renders this
report and its chart from the existing native summary without recomputing
statistics or rerunning experiments.
