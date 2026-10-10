# Synthetic performance studies

This standalone, standard-library Python suite answers three bounded questions:

1. Does a change in a finite random process reduce its operation-count distribution?
2. Does a numerical evaluation change improve throughput with exact reference agreement?
3. Does that improvement survive initialization, host serialization, verification, process overhead, and retained failures?

The implemented candidates are **known positive controls**, not novel algorithms:

- **Coverage:** nearest-neighbor random walking versus independent uniform refreshes on rings of 32 and 128 states. The task is visiting every state, not finding collisions or recovering a logarithm. The transition law changes; this is not an interchangeable cryptographic walk implementation. One transition is the operation-count unit, not a claim of equal machine work.
- **Arithmetic:** direct modular power-sum evaluation versus Horner evaluation for degree-24 polynomials over the fixed field F_257. An independent Python integer expression checks every field element for each polynomial. No cryptographic-sized arithmetic is present.
- **Full cost:** identical generated inputs, randomized within-pair execution order, separate fresh worker processes, host encode/decode phases, complete output checks, separate memory passes, and bounded child watchdogs.

There are no curve/public-key inputs, external targets, imported problem instances, discrete-log solvers, or GPU kernels. These measurements do not establish improvements to PSCKangaroo, JeanLucPons/Kangaroo, or oritwoen/kangaroo. Those repositories' published speed claims motivated the measurement questions only.

## Run and reproduce

From the repository root, using Python 3.10 or newer:

```sh
python3 -m unittest discover -s research/synthetic_performance -v
python3 research/synthetic_performance/study.py --output /tmp/synthetic-study-new
python3 research/synthetic_performance/report.py /tmp/synthetic-study-new
```

The output directory must not already exist. Runs write a manifest first, flush every trial to `trials.jsonl`, and write the final summary only after all scheduled trials have been attempted. A process interruption can leave an incomplete ledger without a summary; such a directory is not a completed study. Failed arithmetic comparisons return a nonzero exit code after retaining the observations. A statistically inconclusive performance result is a valid completed study and does not fail CI.

Defaults are 128 coverage trials per split/size/variant (1,024 total), 24 arithmetic pairs per split (96 timing workers), and one separate memory worker per split/variant (4 workers). Each arithmetic timing worker evaluates 4,096 inputs. Configurations are bounded, with a 15-second watchdog per child by default. Seeds are derived from study/split/instance identifiers. The discovery and holdout fixtures are disjoint; no candidate or parameter is selected using holdout results.

Optional plotting uses matplotlib, tested locally with version 3.10.3:

```sh
MPLCONFIGDIR=/tmp/synthetic-study-mpl python3 research/synthetic_performance/report.py /tmp/synthetic-study-new --plot
```

The report command verifies the ledger hash, study-source hash, and regenerated summary before writing the report. It prefers a run-local `study_source.py` snapshot when present. It uses the current summarizer; incompatible future summarizer changes should be accompanied by an explicit schema migration. Dependency installation is not necessary for measurements, tests, or Markdown reports.

## Statistical and measurement contract

- Retain all capped coverage trials. Report empirical completion quantiles only when the full empirical completion CDF reaches the requested probability. A missing quantile is unknown, not the cap and not a successful-trials-only quantile.
- Compare restricted means E[min(T, 8192)], not inferred uncensored means. A ratio under this cap does not estimate unrestricted completion-time speedup.
- Report baseline/candidate ratios of paired sums, with descriptive 95% percentile bootstrap intervals from 1,000 deterministic resamples. Intervals do not correct for multiple exploratory comparisons and are not a universal significance claim.
- Treat arithmetic correctness and pairing as gates. Any failed, timed-out, or incorrect timing or memory trial blocks a speedup claim for that split. Paired fixture and output hashes must match.
- A finite end-to-end improvement requires a lower confidence bound above one for the parent-measured full wall envelope. Compute-phase improvement alone is insufficient.
- Full arithmetic timing includes input generation, host encode/decode, warmup, computation, output encoding, verification, startup/teardown, subprocess IPC, and parent receipt parsing. Final report generation and artifact writes lie outside per-trial timing. Suite duration is recorded before writing the final summary.
- Memory instrumentation runs separately; its timings never enter speed ratios. Python traced peak and process high-water RSS measure different things. RSS includes the interpreter. One diagnostic memory sample per split/variant does not establish a memory speedup.
- **Transfers are host serialization and process IPC only. GPU/PCIe/network transfer costs are unmeasured.** CPU frequency, affinity, thermals, and other host workloads are uncontrolled. Fresh subprocesses do not imply cold OS caches. Paired random ordering reduces, but cannot eliminate, timing confounding.

## Archived local evidence

- [Final study report](results/final-20260911/REPORT.md): final timing boundary includes parent receipt parsing; complete raw trials, manifest, summary, and holdout chart are retained.
- [Initial study report](results/local-20260911/REPORT.md): retained preliminary run; the parent timer originally stopped before receipt parsing. Its exact source snapshot is archived and hash-checked. This run is not pooled with the corrected run and has not been discarded because a confidence interval was inconclusive.

The second run was triggered by the timing-accounting correction, not by choosing a favorable result. Neither run establishes a cryptanalytic result or an improvement to a production recovery tool.

## Validation and CI

The focused workflow runs correctness, censoring, identity-control, wrong-output, fixture-mismatch, real child timeout, and real child failure tests, followed by a small measurement smoke run. It uploads the smoke receipts even on failure. CI enforces correctness and complete evidence, not a speed threshold on shared runners. No repository runtime code is changed.
