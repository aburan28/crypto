# AES: cross ratio superbox

**Blocked proposal, not an executable trial.** The original advisory card is `cross-ratio-superbox` in [hypotheses.json](../../research/aes_distinguisher_20261017/hypotheses.json). Its literature/source identifiers resolve in [sources.json](../../research/aes_distinguisher_20261017/sources.json). Novelty is unverified; no scientific status is promoted.

## Question and decision

Object: quartet cross-ratio plus the information needed across real AES affine/mixing layers.

Prediction: local inversion identity survives its stated hypotheses but requires additional coupled information across a complete round.

If this prediction survives the specified independent controls, extend the exact constraint model for this mechanism. If its falsification condition holds, revise the observable before attempting larger AES variants. Retain per-key replication and complete costs separately from any claim of scientific novelty. Acceptance: complete published predicate retrieved and implemented, value-level counterexamples or independent exact witnesses.

## Proposed discriminating experiment

exhaustive GF(16) quartets and <=2^20 selected GF(256) quartets, then <=2^18 multi-byte structures/key on <=16 keys.

Controls:

- inversion alone
- actual AES output affine map
- affine undone
- random bijection
- denominator/pole and ordering checks

Falsification: identity fails inside its declared domain, apparent survivor matches nulls, or retained object reconstructs the whole input tuple.

Lossy-projection obligation: test equal-cross-ratio distinct quartets through each real layer; do not substitute a lossless pair representation.

## Prerequisites and resources

Current implementation boundary: local component ablation implemented; super-box extension not implemented. Implement the extension and an independent native checker, then freeze the exact protocol, source revision, inputs, sampling and alpha allocation before marking this card ready. The argv/build names in card.json are proposed interfaces; their binaries do not yet exist and must not be submitted to a runner.

component observations do not authorize whole-round or full-AES claims. The card proposes one CPU and 1 GiB RAM, with a one-hour process watchdog as machine protection. Time estimates are unknown, not a claim that the experiment fits that watchdog; derive them from the implemented operation counts and a measured pilot before admission. Retain censored cells. Funding units rank these nine proposals and do not transfer money or impose a cap on the research program. No cairn environment, task or lease is created.

## Done means

complete published predicate retrieved and implemented, value-level counterexamples or independent exact witnesses. Preserve negative, failed and inconclusive outcomes. A software pass or component observation does not establish a full-AES attack, universal key property or official campaign admission. Existing GOAL-AES-001/002/003 controls remain separate.
