# Follow-up validation and accounting

Measured source commit:0cbc86b2c828f1aa7fae62282913cce9ea6e099d.
Measured binary SHA256:50a5739fadcb08ab8f170d4fe03ba3427fb0fd387c1a5918843fcda4d68ff1fb.
The protocol, specs and algorithm sources were committed before measurements.
The dispatcher command manifest records every run, audit, table, resource and
comparison invocation, including the initial wide-reader failures.

Six immutable sessions contain373 attempts:336 verified,36 exhausted and one
legacy-width error. All336 verified outcomes pass standard exact replay; all37
unsuccessful outcomes pass the native supplemental replay with the exact original
measured binary. The wide session has no successful target and no matched rho
cost. Its tiny S values are censored work, not attack costs.

The N83 session's initial reader failure, original table error and resource error
are retained. The corrected reader uses the historical power table and scaling
algorithm. The versioned receipt is `wide-n83.audit-v2.json`; `wide-n83.audit.json`
is an exact copy so the common analysis can locate it. The old commands produced
no audit/resource JSON. `wide-n83.table-v2.md` preserves the corrected table beside
the failed original. Auditor SHA256 is
5375956bdbe597a4d710455ab81ea3c3bc2394fc5e055665a0e2334a9647fe43;
the auditor was built with Rust1.98 and the reader correction, while measurement
used Rust1.99. The14 reference vectors come from independently pinned serde_json
1.0.145 with only std, outside root feature unification. Existing v1 bit patterns,
signed zero, tiny/subnormal and large finite values remain checked.

The frozen method implementation, seed law, point set and counts are unchanged.
No tolerance, integer comparison, record seal or replay gate was relaxed. Source
and commands for the reader correction are in AUDIT-READER-AMENDMENT.md and the
decimal-reference crate. New analyses are Rust; Python draws/layouts the recorded
results and renders the existing catalog prose. It does not run experiments or
derive new measurements, ratios, statuses or research verdicts.

Host class:ECBENV2h96c931cc33fa; AMD EPYC9V45, x86_64, five exposed logical CPUs,
managed virtual machine, L0. The environment-class digest differs from the parent
campaign, but this is not an independent-host replay or an isolated hardware
performance result. Native work is unpriced. No exponent, full-cost ratio,
wall speedup, new theorem or N131 transfer is claimed.

Eight constructed algebraic systems and16 paired old/v2 solver checks validate
shape metadata. The rank panels comprise56 IC recoveries and40 matched rho/A-A
controls. All224 prime-panel attempts verify. The4 N83 attempts exhaust query
caps with zero checked relations and aggregate rank0; the previous pipeline's
refusal is retained as its own control.

Library/integration tests, formatting, lint, generated-metadata checks, catalog
checks and visual inspection receipts are recorded beside this document. CI and
review state are reported separately; local correctness is not merge approval.

The corrected pre-integration full library invocation completed with4095 passes,
zero failures and117 preexisting ignored tests. The earlier4094/one-failure run
is retained as `test-library-v1-failed.log`. The new actual-width point witness
and all14 independent decimal-token reference cases pass.
