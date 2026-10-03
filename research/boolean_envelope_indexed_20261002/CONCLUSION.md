# Indexed envelope: correctness verified, no dramatic discovery gain

The degree-indexed and occupancy-based envelopes reconstruct every declared
Boolean matrix exactly. The native Linux ARM64 discovery is complete,
resource-qualified and hash-bound in `qualified_discovery_01`: 128 cells,
8,960 A/B observations, 2,560 A/A observations and 244,800 verified output
matrices. One earlier build-host PSI refusal remains sealed in
`failed_discovery_01` and supplies no admitted timing.

The degree-only arm passed **0/32** incremental and **0/32** dramatic discovery
groups. The indexed arm passed **16/32** incremental and **0/32** dramatic
groups against the frozen five-arm reference, after the A/A noise gate.
Discovery cannot pass the preregistered primary holdout gate; no holdout was
executed. `BOUNDARY_CORRECTION.md` explains why even a positive comparison
against these five arms would not establish a gain against packed direct.

| Arm | n12 coefficients, batch 64, cold ms | Fastest old / arm | Correctness |
|---|---:|---:|---|
| Direct | 5.890810 | 0.983 | PASS |
| Verified layout reuse | 5.834268 | 0.992 | PASS |
| Exact product schedule | 5.795539 | 0.999 | PASS |
| Exact packed-matrix cache | 5.789168 | 1.000 | PASS |
| Original support envelope | 6.859047 | 0.844 | PASS |
| Degree-indexed envelope | 6.787266 | 0.853 | PASS |
| Indexed envelope | 5.054149 | 1.145 | PASS |

Values are the means of the two fixed seed medians, in milliseconds per
complete batch including validation. The ratios are descriptive and use the
pooled fastest of the five old arms; the decision uses each of 32 separately
paired lower confidence bounds and its A/A floor in `qualified_discovery_01/results.json`.
No old-host absolute time is used as a speed ratio here. The counted n12
quadratic plan visits fall from 299 to 13, but that mechanical saving did not
produce a 2x cold construction result across the discovery groups.

This is a bounded mathematical construction study on generated public Boolean
systems. Its full-IC cost, calibrated-operation ratio, relation-yield effect
and rho ratio are null. There is no key-related result or cryptanalytic
breakthrough.
