# Result: fixed-addition rho screen toward 12 B/s on G7

## Decision

Reject the non-equivariant fixed-add walk and retain the selected self-Frobenius
walk. Fixed add removes the runtime coordinate transforms and raises the raw
kernel median from **6.321422** to **7.601382 billion scalar updates/s**, but it
requires **1.399081x** as many complete iterations in the matched GF(2^23)
collision study. The collision-adjusted candidate rate is **5.433126 billion
iterations/s**, below selected.

## Matched G7 screen

Each row executed 34,359,738,368 complete scalar updates on the same local
AWS g7.2xlarge, NVIDIA RTX PRO 4500 Blackwell Server Edition, at a 165 W power
limit. Process observation found no external GPU process.

| repetition | selected B/s | fixed-add B/s | raw ratio |
|---:|---:|---:|---:|
| 0 | 6.381188 | 7.619793 | 1.194103 |
| 1 | 6.321422 | 7.601382 | 1.202480 |
| 2 | 6.300687 | 7.543261 | 1.197212 |
| **median / paired result** | **6.321422** | **7.601382** | **1.197927** |

The paired 95% Student-t interval on log speed ratios is **[1.187458,
1.208488]**. The raw improvement is reproducible but reaches only 63.34% of the
12 B/s target.

## Collision work and adjusted result

Two thousand planted GF(2^23) solves per mode used identical seeds, DP threshold,
eight-way partition and orbit key. Every scalar was recovered.

| mode | solves | bad | mean iterations | median | p90 |
|---|---:|---:|---:|---:|---:|
| selected self-Frobenius | 2,000 | 0 | 172.331 | 165 | 291 |
| fixed add | 2,000 | 0 | 241.105 | 229 | 422 |

The candidate/selected mean-work ratio is **1.399081**. Dividing the raw paired
ratio by this measured work ratio gives **0.856224**, with 95% interval
**[0.848741, 0.863773]**. The fixed points break covariance under Frobenius, so
walks in the same orbit no longer follow the same quotient path. This explains
the collision regression observed in the matched study.

Generic work remains `sqrt(n/262) * 1.399081` for this rejected candidate. A
full-DLP S value remains null because no end-to-end recovery was measured.

## Correctness and audit evidence

- The final built-in suite passes all supported curve and collision-algebra rows.
  The two earlier fixture failures and their corrected binaries/logs are retained.
- An independent full-state oracle matched 2,048 GPU points across seven steps,
  for zero coordinate mismatches.
- The DP client emitted 8,406 records; 16 records replayed on CPU, and the
  checkpoint and DP multiset hashes are frozen in `validation.json`.
- Compute Sanitizer memcheck, initcheck and synccheck pass.
- The GF(2^23) collision study recovered all 4,000 planted scalars.

The experiment was closed by the collision-adjusted regression, so promotion-only
partial-population, restart-boundary and fresh DP34 confirmation gates were not
run. The selected production binary and production sources were not changed.

## Reproducibility

`comparison.json` contains the machine-readable paired calculation. Raw timing,
telemetry, process receipts, collision output, validation logs, final isolated
source manifest, candidate patch and both corrected build attempts are retained
in this directory.


## Covariant fixed-add follow-up

A preregistered follow-up made the fixed addend rotate with a phase extracted
from the normal-basis x coordinate. Exhaustive GF(2^23) verification passed all
8,388,606 nonexceptional coordinates and proved the transition commutes with
Frobenius. It did not pass the CPU continuation gate.

With the existing restart mechanism set to 1,024 steps, the matched 100-trial
selected control averaged 166.640 complete iterations; the covariant candidate
averaged 197.260, a **1.183749x** work ratio. All scalars and checked full states
were correct, but one walk required restart and three exact additive collisions
had zero solve denominator. The offending deterministic seed has an exact
four-state DP-free cycle with minimum x weight 12 versus DP threshold 10.

Because 1.183749 exceeds the preregistered 1.10 limit, the follow-up was closed
before a 2,000-trial study or CUDA implementation. See `FOLLOWUP.md`,
`covariant-comparison.json`, and `covariant-cycle-witness.txt`.
