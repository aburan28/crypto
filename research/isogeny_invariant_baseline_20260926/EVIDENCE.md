# Evidence required for this benchmark

An assistant's answer is not a receipt. A hard-coded expected number is a regression fixture, not independent verification. A SHA-256 hash binds bytes; it does not establish that a calculation happened or that the mathematics is correct.

## What is now externally inspectable

The `isogeny-invariant-baseline` GitHub Actions workflow checks out the exact PR head, runs the benchmark from scratch, checks the archived and new summaries against the frozen result, then runs `verify_independent.py`. It uploads an evidence artifact even on failure. Follow the PR's **exact-replay** check to its public run and download `isogeny-evidence-<run-id>-<attempt>`. Artifacts have a 90-day retention period; the code and original raw baseline remain in Git. The receipt is not signed by a theorem prover, and artifact retention is not permanent archival.

Each evidence bundle contains:

- `receipt.json`: source commit, worktree cleanliness, source and input SHA-256 hashes, CI run URL and attempt, Python/platform metadata, commands with exit codes, start/end times, per-command diagnostic wall time, artifact hashes, negative controls and explicit scope limits.
- `fresh.json`: the independently hosted run's complete outputs, including zero-coverage targets.
- `summary.json`, `archive-summary.json`: fresh and archived audits.
- `independent.json`: results from the second implementation.
- `command-*.stdout` and `command-*.stderr`: actual subprocess logs, plus `failure.txt` if a check fails.

A pass requires equality of all three audited summaries and rejection of three deliberately corrupted inputs: a pair count, a curve point, and a Frobenius orbit count. If the run fails, its receipt is marked `failed`; the same bundle must not be presented as successful evidence. Python optimization is rejected because the original benchmark uses assertions.

## Independent finite computation

`verify_independent.py` imports no benchmark arithmetic. It constructs F81 arithmetic using explicit polynomial reduction and exhaustive inverse search, enumerates every curve point, and computes sums by removing roots from the cubic obtained by intersecting the curve with a chord or tangent. It reconstructs the fraction set using all nonzero denominators rather than the benchmark's normalized projective denominator representatives.

It recomputes the full per-target pair-count vectors and distinct-x vectors for **246 factor-base rows** (41 bases on each of six curves), including field/point orbits for invariant bases. This is stronger than checking aggregate totals or expected headline numbers. Its three corruption controls must fail. Both implementations were written by the same assistant, so correlated specification errors remain possible; external review or a third mathematical system would add assurance.

## Claim labels

| Label | Required evidence | Present here? |
|---|---|---|
| Proposed | Explicit hypothesis and scope | Yes |
| Executed | Public CI run, source/input hashes, raw outputs and logs | Workflow produces this on each run |
| Independently recomputed | Separate algorithm reproduces complete finite outputs | F81 incidence only |
| Mathematically proved | Reviewed derivation or proof-assistant certificate with stated assumptions | Not supplied by these computational receipts |
| Performance improvement | Matched complete workload, accounting and held-out results | Not established |

The F13 degree-7 construction remains checked by the original implementation, not this independent verifier. Neither verifier establishes an efficient quotient solver or a full-DLP speedup.

## Reproduce locally

From a clean checkout:

```sh
python3 research/isogeny_invariant_baseline_20260926/replay_with_receipt.py \
  --require-clean --output /tmp/isogeny-evidence-new
```

Use a new output directory each time. Inspect `receipt.json`, verify its source commit against the PR head, recompute its file hashes, and require successful commands, passing negative controls and the independent verification scope before citing the result. For future calculations, reports should link the run and receipt alongside every numerical claim rather than saying only “tests passed.”
