# Review guide

The generated catalogs and raw evidence make this PR large. Review the maintained
implementation and claim boundaries first; the generated files are independently
rebuildable and the raw sealed records remain available for audit.

| Area | Source | Main invariant |
|---|---|---|
| CM arithmetic | [`cm_orbit.rs`](../../src/cryptanalysis/ic_framework/cm_orbit.rs) | Derived degree-5/11 maps act as the recorded subgroup scalar; direct scalar control uses the same public seeds. |
| Native adapter | [`cm_folded.rs`](../../src/cryptanalysis/ecbench/cm_folded.rs) | Required binary hash; independent public-base, relation, rank and target replay; five exclusive online phases; unpriced work stays explicit. |
| Public fixtures | [`workload.rs`](../../src/cryptanalysis/ecbench/workload.rs) | Measured-child transport preserves G and Q; prime nonsingular cofactor-one curves, subgroup membership and Hasse-unique order checks reject malformed inputs. |
| Generic control | [`folded_recurrence.cpp`](native/folded_recurrence.cpp) | P[k+2] = P[k+1] - [5]P[k] generates the same orbit without a floor CM map. |
| Statistics | [`campaign_analysis.rs`](native/campaign_analysis.rs) | Only sealed sessions; all outcomes retained; exact same-crater comparisons; whole-gap bootstrap with recorded seeds. |
| Independent IC replay | [`replay_campaign.rs`](native/replay_campaign.rs) | Reconstruct the base by direct scalar arithmetic, check SHA-256, rows, ranks, query indices and [d]G=Q; reject a corrupted sign. |
| Fast codomain | [`sqrt_velu.cpp`](native/sqrt_velu.cpp), [`batched variant`](native/sqrt_velu_batched.cpp) | Kernel-moment identities, complete known-reference comparisons, explicit caps; no point-map or dual-map claim. |
| Website admission | [`build_lab_browser.py`](../../scripts/build_lab_browser.py), [`regressions`](../../scripts/test_lab_browser_cost_admission.py) | An unpriced zero charged bound is unknown, not an eligible zero-cost speed score. This is UI metadata tooling. |

## Results to cross-check

- [`paired-primary.json`](results/paired-primary.json): 500 verified attempts;
  native CPU ratio 0.980 [0.866, 1.111]. Endpoint waves and shared-host activity
  confound CPU timing. This interval is not a causal hardness estimate.
- [`replication-metrics.json`](results/replication-metrics.json): 265 verified
  attempts. The ten-pair D=-43 signal is exploratory, unadjusted for multiple
  comparisons and accompanied by extra queries.
- [`posthoc-replay.json`](results/posthoc-replay.json): 473 independently
  replayed IC transcripts. [`target-replay.jsonl`](results/target-replay.jsonl):
  797 independently checked recoveries, including repeated controls.
- [`isolated-summary.json`](results/isolated-summary.json) and
  [`local reservation receipt`](results/isolated-receipt.jsonl): CM 5.464 ms,
  floor recurrence 6.857 ms, direct floor scalar 42.141 ms. Local CPU reservation
  does not establish physical-host or frequency isolation.
- [`Lean sources`](lean/README.md): three conditional lemmas with no axioms;
  CM classification, rational-map homomorphism and order-lattice bridges remain
  explicitly unformalized.

## Reproduction boundaries

The protocol and frozen binaries define the historical experiment. Current
sources add stronger fixture checks after the completed full library run;
targeted regressions and pinned Rust 1.98 strict Clippy pass on those changes.
The new CI workflow independently replays evidence and runs fresh small controls;
fresh CI runs do not replace or silently amend the historical measurements.

Use [`native/README.md`](native/README.md) for backend builds and
[`evidence-manifest.sha256`](evidence-manifest.sha256) for checkpoint hashes.
Run `sha256sum -c` with that manifest from the repository root. The PDF is
rendered by `python native/render_report.py` from the study directory; research
aggregation, replay and statistics are native Rust. Resource-limited upper
construction attempts are retained, and a floor is counted only after an
independent PARI certificate. All completed large-gap ECDLP timings are 53-bit.

Full common-cost IC/rho comparisons, a prime-field algebraic IC backend,
the original non-rational-kernel 96-bit floor, and the listed formalization
bridges are open research obligations. See the [report ledger](REPORT.md).
