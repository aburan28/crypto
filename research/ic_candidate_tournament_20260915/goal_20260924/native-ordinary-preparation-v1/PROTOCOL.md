# Native ordinary preparation: independent mathematical gate v1

Status: implementation and retained-data controls only. No new solver has run,
no new ordinary-query panel is admitted, and no fresh comparison is registered.
This protocol does not reopen any consumed SAT/F5 registration or confirmation set.

The objective remains complete, source-bound one-target F4/F5 and CryptoMiniSat
IC pipelines, natural relation accounting, then a fresh paired comparison with
the incumbent and a strong rho reference. This gate addresses preparation
mathematics; it does not replace the source, runtime, timing or comparison gates.

## Hypothesis and boundary

A preparation transcript can independently establish which ordinary queries
have verified geometric witnesses, which have exact geometric negatives, which
remain inconclusive, and whether its projected relations determine all factor
logs. Solver status, relation count and a claimed rank alone cannot establish
this. A duplicate relation or dependent row contributes no novel rank.

The mathematical acceptance boundary is full rank over the exact subgroup
modulus with independently replayed column logs, a complete preregistered query
panel, and an explicit producer log claim matching the independent solution.
Reject a false witness, false negative, inconsistent row, wrong input law,
missing chronology, incomplete base or unsupported log claim. Preserve an
interrupted panel and every declared failure; full rank alone does not complete
that panel. There is no performance hypothesis or speedup acceptance threshold
for a data audit. Online wall time and speedup remain null.

## Exact disclosed synthetic mathematics

Only the educational n17 Koblitz fixture is admitted: binary polynomial basis
with modulus 131081, curve coefficient a=1 and b=1, group order 131174,
prime subgroup order 65587, cofactor 2, generator (43693,23339), and Frobenius
eigenvalue 17184. Targets and target seeds are empty; no target scalar is
constructed. The independent bit-arithmetic oracle verifies irreducibility,
curve/group order, generator and Frobenius action. No challenge or third-party
target is an input to this protocol.

Reconstruct all affine points with x in the six-dimensional polynomial
subspace, x<64, rather than trusting a nominal dimension or supplied count.
For x=0 the sole point is (0,1). For nonzero x, solve z²+z=x+a+1/x² with
the odd-degree half trace, test it, and include both y=xz and y+x when they
exist. The supplied ordered list must contain exactly the resulting 63 points
without duplicates. Cofactor projection has 62 distinct nonidentity images;
sign/Frobenius folding produces 29 independently reconstructed columns.
Retain the ordered geometric-set digest and actual column points/projections.
These counts are separate; the 63-point geometric list is not an IC1 fb count.

## Input and independent audit

`icprog ordinary-preparation-audit --input INPUT.json --out NEW.json` reads at
most 8 MiB, starts no solver or subprocess, rereads the exact input bytes before
publishing, and creates its output without overwriting evidence. Use the
repository busy wrapper for builds and audits. All implementation and controls
are Rust. The immutable old F5 certificate is read only in tests; its historical
Python provenance and original seal remain unchanged.

The strict JSON input contains exactly:

* `schema_version: 1`, `question: "ordinary-preparation-mathematics-n17-v1"`;
* `plan`: `family` (`matrix_f5` or `cryptominisat`), integer `algorithm_seed`,
  `planned_queries` between 1 and 512, and
  `query_law: "independent-probe-scalar-rand08-v1"`;
* the exact target-free `fixture`, ordered integer-coordinate `geometric_base`;
* chronological `attempts`, each with integer `trial`, integer `scalar`,
  `outcome`, and `indices` (three geometric indices for a witness, null otherwise);
* `stop: "panel_complete"` exactly when the planned count is retained, otherwise
  `"interrupted"`, and `claimed_column_logs` (integer list or null).

Unknown fields and JSON floats are rejected. The family is a declaration, not
source admission. Per trial, reconstruct the production collector's keyed
rand-0.8 StdRng law with an independent PCG seed expansion/ChaCha12/rejection
sampler. The audit never calls the producer sampler or its curve arithmetic.
The trial key is `seed XOR 0x50524f4245534551 XOR
rotate_left_17(trial * 0x9e3779b97f4a7c15 mod 2^64)`; sample uniformly in
1..65587 with the registered rejection rule.

Re-add each three-point witness to `[scalar]G` before cofactor projection.
Exact negatives use an independent finite pair-sum membership oracle,
including repeated points, torsion and identity pair sums. Retain `incomplete`,
`unsupported`, `invalid_model`, `unresolved`, `timeout`, `transport_failure`
and `oom` without turning them into negative results or rows. Independently
feasible inconclusive inputs are diagnostic missed witnesses, never wins.

Reconstruct each augmented relation modulo 65587 with RHS `2*scalar`.
Match the existing raw duplicate policy `(scalar, sorted geometric indices)`;
also count repeated projected augmented rows separately. Independent modular
elimination records rank after every attempt and rejects inconsistent rows.
Back substitution is permitted only at full rank. Verify every recovered
column log by scalar multiplication and compare any producer claim exactly.

Passing output explicitly reports `source_bound_execution_admitted: false`,
zero native solvers and new queries, false fresh/headline/promotion/full-goal
flags, and null online cost/speedup. A passing partial panel is mathematical
consistency only. This input schema carries no solver models, costs or runtime
custody and therefore cannot independently admit them.

## Meaningful controls before integration

Required Cargo controls compare the independent sampler with the existing
sampler across seeds/trials; replay the retained 216-attempt F5 certificate
(61 witnesses, 155 exact negatives, 32 dependent rows, rank 29); reject changed
laws, chronology, geometry, targets, witnesses, negatives and logs; preserve a
partial full-rank panel and timeout; distinguish duplicate and dependent rows;
reject inconsistent equations, unknown fields and floats. This is replay of
existing disclosed data, not new natural yield. Keep validation failures and
their fixes in the PR record.

## Subsequent execution gates, still pending

1. Build a genuinely target-free native producer for each exact pipeline.
   The current full-DLP worker requires a target and is not this producer.
   MatrixF5 must identify its actual bounded engine and limits; do not describe
   it as a complete incremental F5 algorithm. CryptoMiniSat must use the
   accepted external source/exporter/assets, not the worker's Rust CDCL backend.
2. Freeze a new protocol, full source/dependencies/native assets, binaries,
   limits and one-use registration before dispatch. A proposed common natural
   panel is 512 independent queries with seed 2026100311; that proposal is not
   an executable registration. Complete the fixed panel even after rank 29;
   rank-based early stopping changes the sampling question. Preserve interruption
   rather than refill or silently repeat missing trials.
3. Retain original source models for SAT witnesses and original solver receipts,
   independently check encoding/model/group agreement, and verify exact negative
   claims. Retain every failure and per-attempt integer phase costs. Native
   custody and frozen postexecution audit must bind all of these to the producer.
   Report verified yield and novel-rank rates with uncertainty on the frozen
   law; no planted decomposition estimates natural yield. Do not claim
   full-source runtime admission from the mathematical gate above.
4. Complete native preparation with no imported factor logs or target answer,
   source-bound rank/log replay, setup/base/memory/matrix/LA costs and failed
   attempt accounting. Only then freeze a new fresh public-point protocol with
   an exposure census, complete IC1/workload identities and source-bound
   prepared state. Existing exposed control points remain controls.
5. Use native ecbench and a calibrated resource envelope for paired independent
   one-target solves by F5, CryptoMiniSat, the incumbent and strong rho. Measure
   target-dependent online work through independently verified scalar recovery,
   with five exclusive online phases summing to the charged interval. Keep
   reusable preparation separate. Preserve all failures; report paired target
   costs, ratios and uncertainty. No throughput average replaces this question.
6. Maintain diverse pipeline candidates and interaction tests rather than
   promoting isolated stage winners. A stage winner cannot establish a global
   optimum. Fresh confirmation can support only the fastest verified pipeline
   within its disclosed tested portfolio, instance distribution and resources.

These six steps remain open. This implementation neither admits a new ordinary
panel nor completes the broader goal.
