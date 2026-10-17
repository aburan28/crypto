# AES: sat witness and unsat

**Blocked proposal, not an executable trial.** The original advisory card is `sat-witness-and-unsat` in [hypotheses.json](../../research/aes_distinguisher_20261017/hypotheses.json). Its literature/source identifiers resolve in [sources.json](../../research/aes_distinguisher_20261017/sources.json). Novelty is unverified; no scientific status is promoted.

## Question and decision

Object: exact value-level round constraints and checkable witness/UNSAT certificates.

Prediction: the backend accepts every real S-box pair, rejects every incorrect pair and agrees with full toy enumeration.

If this prediction survives the specified independent controls, extend the exact constraint model for this mechanism. If its falsification condition holds, revise the observable before attempting larger AES variants. Retain per-key replication and complete costs separately from any claim of scientific novelty. Acceptance: external certificate replay and complete toy agreement; timeout kept separate from UNSAT.

## Proposed discriminating experiment

all 8-bit S-box relation pairs; <=32 toy cells with <=5 rounds, 300-second watchdog/cell.

Controls:

- independent direct evaluator
- planted witnesses
- inconsistent constraints
- random S-box
- round/state-order fixtures

Falsification: a witness fails direct evaluation, a valid witness is excluded or preprocessing changes the predicate.

Lossy-projection obligation: full-state SAT is lossless infrastructure, not a new cryptanalytic object.

## Prerequisites and resources

Current implementation boundary: not implemented. Implement the extension and an independent native checker, then freeze the exact protocol, source revision, inputs, sampling and alpha allocation before marking this card ready. The argv/build names in card.json are proposed interfaces; their binaries do not yet exist and must not be submitted to a runner.

<=2 GiB/process; no full-AES conclusion from bounded failure. The card proposes one CPU and 2 GiB RAM, with a one-hour process watchdog as machine protection. Time estimates are unknown, not a claim that the experiment fits that watchdog; derive them from the implemented operation counts and a measured pilot before admission. Retain censored cells. Funding units rank these nine proposals and do not transfer money or impose a cap on the research program. No cairn environment, task or lease is created.

## Done means

external certificate replay and complete toy agreement; timeout kept separate from UNSAT. Preserve negative, failed and inconclusive outcomes. A software pass or component observation does not establish a full-AES attack, universal key property or official campaign admission. Existing GOAL-AES-001/002/003 controls remain separate.
