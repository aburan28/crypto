# Selector admission from the verified W64 cold phase ledger

**Decision: prioritize a single-index, low-overhead factor-base or probe policy at n41/L1,024; defer a two-full-index per-instance selector.** Class of change: **accounting** under AGENTS.md §3. This is a cost-admission calculation from the already verified [six-cell v2 outcome](../disjoint_cold_v2_outcome_20261001/RESULT.md), not a measured selector and not an ECC2K-130 speed claim. The one-shot n131 m10 representation gate remains separately pending in [PR #1112](https://github.com/aburan28/crypto/pull/1112); these W64 source-curve timings do not predict native-leaf PDP yield.

The archived n41 and n53 batch cells used the same public Q within each `ic_a`/rho/`ic_b` block. Both compact arms in every block solved 1,024 Q, reached full rank in exactly K successful attempts with zero failed or nongaining rank rows, and agreed on the base hash. [derive.py](derive.py) verifies the SHA-256 of both raw tarballs and every consumed member against the archived manifest, replays ten compact phase records and fifteen `wait4` CPU receipts per cell, and confirms the five-block paired CPU medians in the published analysis.

| Verified cell | Median complete IC/rho CPU | Median index-build share of compact elapsed time | Median rank + target share | Rank + target reduction required with **zero** selector overhead | Required reduction if a second full index is built |
| --- | ---: | ---: | ---: | ---: | ---: |
| n41 / L1,024 / K255 | 1.3137× | 35.27% | 63.23% | 37.86% | 93.97% |
| n53 / L1,024 / K440 | 1.6054× | 28.40% | 70.61% | 53.85% | 93.54% |

For each paired block, let `R` be complete compact/rho child CPU, `I` the geometric-mean index-build elapsed-time share of compact `process_total`, and `V` the corresponding rank-plus-target share. A counterfactual selector that leaves the other phases unchanged needs variable-stage saving `s ≥ (1−1/R)/V` if it has no overhead. If it builds **one additional full index** before choosing a base, the threshold becomes `s ≥ (1+I−1/R)/V`, even when its scan, score, selection, memory and output are free. The table reports medians of the five blockwise thresholds, not thresholds calculated from independently pooled medians. [RESULT.json](RESULT.json) retains every block and its ranges.

This calculation assumes internal elapsed-phase shares proxy charged CPU shares. That assumption is plausible for the pinned single-thread children, but is **not** a measured per-phase CPU attribution or a lower bound on all selector implementations. A different base may build a faster index; it may also change the rank and target costs. The second-index column is therefore a deliberately favorable *scenario requirement*, not a mathematical impossibility or a measured negative experiment. The zero-overhead column is likewise an optimistic planning threshold. Any actual selector must charge its scan, score, chosen index, rejected candidates, rank, all Q descents, and scalar checks in one cold process and compare against same-Q rho. A target-blind score learned from these v2 Q would be invalid on them; new orbit-disjoint Q and a separately frozen score are required.

The next bounded experiment should use n41/K255/L1,024 first and build **one** complete index for the selected 255 useful signed-orbit columns. Candidate scoring must precede index construction, use only a separately frozen training stream or public base features, and preserve a same-Q default-base control and matched rho. Require all arms to recover every Q and full rank; compare complete CPU with the previously measured 1.3137× gap as the admission reference. If n41 does not improve its complete ratio under this charged policy, do not spend an n53 batch on the same selector. If it does, preregister n53/K440 and then the n83 wider-backend confidence gate; no n131 transfer follows from n41 alone. Preserve negative, timeout and censored runs and independently replay all recovered scalars.

Reproduce the deterministic phase calculation from the repository root:

```sh
python3 research/notes/ecc2k130/selector_phase_admission_20261001/derive.py \
  --out /tmp/ecc2k130-selector-phase-admission.json
cmp /tmp/ecc2k130-selector-phase-admission.json \
  research/notes/ecc2k130/selector_phase_admission_20261001/RESULT.json
```
