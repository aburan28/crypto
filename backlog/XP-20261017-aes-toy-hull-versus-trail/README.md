# AES: toy hull versus trail

**Blocked proposal, not an executable trial.** The original advisory card is `toy-hull-versus-trail` in [hypotheses.json](../../research/aes_distinguisher_20261017/hypotheses.json). Its literature/source identifiers resolve in [sources.json](../../research/aes_distinguisher_20261017/sources.json). Novelty is unverified; no scientific status is promoted.

## Question and decision

Object: complete endpoint differential distribution and compatible intermediate paths.

Prediction: beam model mass and exact fixed-key endpoint counts differ in explainable, key-dependent ways.

If this prediction survives the specified independent controls, extend the exact constraint model for this mechanism. If its falsification condition holds, revise the observable before attempting larger AES variants. Retain per-key replication and complete costs separately from any claim of scientific novelty. Acceptance: exact count agreement between independent toy evaluators and explicit distinction between model and fixed-key probabilities.

## Proposed discriminating experiment

complete 16-bit toy domain, <=16 keys, <=8 frozen input differences, rounds 1..5.

Controls:

- exact brute force
- random 16-bit permutations
- independent-round-key model
- random bijective S-box

Falsification: claimed excess disappears under full enumeration or depends on incompatible paths.

Lossy-projection obligation: difference/activity forgets the base state and cannot generally propagate deterministically through SubBytes.

## Prerequisites and resources

Current implementation boundary: not implemented. Implement the extension and an independent native checker, then freeze the exact protocol, source revision, inputs, sampling and alpha allocation before marking this card ready. The argv/build names in card.json are proposed interfaces; their binaries do not yet exist and must not be submitted to a runner.

stored states <=2^20; preserve censored cells. The card proposes one CPU and 1 GiB RAM, with a one-hour process watchdog as machine protection. Time estimates are unknown, not a claim that the experiment fits that watchdog; derive them from the implemented operation counts and a measured pilot before admission. Retain censored cells. Funding units rank these nine proposals and do not transfer money or impose a cap on the research program. No cairn environment, task or lease is created.

## Done means

exact count agreement between independent toy evaluators and explicit distinction between model and fixed-key probabilities. Preserve negative, failed and inconclusive outcomes. A software pass or component observation does not establish a full-AES attack, universal key property or official campaign admission. Existing GOAL-AES-001/002/003 controls remain separate.
