# AES: learned native proposal ranker

**Blocked proposal, not an executable trial.** The original advisory card is `learned-native-proposal-ranker` in [hypotheses.json](../../research/aes_distinguisher_20261017/hypotheses.json). Its literature/source identifiers resolve in [sources.json](../../research/aes_distinguisher_20261017/sources.json). Novelty is unverified; no scientific status is promoted.

## Question and decision

Object: learned statistic of differences/quartets compared with classical feature likelihoods.

Prediction: incremental discrimination persists on untouched keys at a fixed false-positive rate.

If this prediction survives the specified independent controls, extend the exact constraint model for this mechanism. If its falsification condition holds, revise the observable before attempting larger AES variants. Retain per-key replication and complete costs separately from any claim of scientific novelty. Acceptance: single untouched final evaluation with charged generation/training/inference and reproducible model snapshot.

## Proposed discriminating experiment

<=2^20 training examples, 2^18 validation examples, <=2^20 final examples over >=64 final keys, <=3 frozen models.

Controls:

- classical histograms
- linear baseline
- shuffled labels
- consistent PRP
- wrong input difference
- raw versus difference-only ablations

Falsification: gain vanishes under key-disjoint testing, comes from generator artifacts or is matched by a cheaper baseline.

Lossy-projection obligation: a learned score is an empirical observable without a deterministic propagation guarantee.

## Prerequisites and resources

Current implementation boundary: native key-disjoint dataset exporter implemented; training not performed. Implement the extension and an independent native checker, then freeze the exact protocol, source revision, inputs, sampling and alpha allocation before marking this card ready. The argv/build names in card.json are proposed interfaces; their binaries do not yet exist and must not be submitted to a runner.

no GPU training launched in this task. The card proposes one CPU and 1 GiB RAM, with a one-hour process watchdog as machine protection. Time estimates are unknown, not a claim that the experiment fits that watchdog; derive them from the implemented operation counts and a measured pilot before admission. Retain censored cells. Funding units rank these nine proposals and do not transfer money or impose a cap on the research program. No cairn environment, task or lease is created.

## Done means

single untouched final evaluation with charged generation/training/inference and reproducible model snapshot. Preserve negative, failed and inconclusive outcomes. A software pass or component observation does not establish a full-AES attack, universal key property or official campaign admission. Existing GOAL-AES-001/002/003 controls remain separate.
