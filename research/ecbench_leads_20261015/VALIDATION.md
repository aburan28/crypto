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

## Integration checks

Integration commit7b0654c768f5b6775918d759446f400a6aebc03d merges exact main
a80b15bbc5c1ea4f2db2b294a9bcddc31f06fdd7 without rebasing or rewriting a sealed
session. Both conflicts concerned generated/index views. The native publisher
preserves main's110 sessions and adds the two campaigns'15, totaling125.

- `cargo +1.98.0 test --release --lib`:4115 passed, zero failures,117 existing
  ignored; `test-library-integrated.log` and exit0. All five orbit-IC tests,
  including the constructed actual-N83 point witness, pass. The witness does
  not assert an unknown-target N83 recovery.
- The standalone publisher's test and warnings-as-errors clippy pass. Native
  metadata freshness, the exact historical-input guard, catalog display,
  curve-name registration and three cover-catalog tests pass. The path filter
  check covers all214 recorded presentation inputs; its nine tests pass.
- The unchanged historical leaderboard and browser checks pass in a detached
  exact-a80 worktree. No new campaign is analyzed by their Python builders.
- Legacy autolab control-plane suite:47 pass. Its first import failure is
  preserved; the retry materializes the tracked `oracle_u128.py` sparse input.
- Every changed Rust file passes `rustfmt +1.98.0 --edition 2021 --check`.
  The first scoped check found only the new decimal-reference layout, corrected
  without changing token values. Whole-repository cargo fmt still reports the
  same two incoming examples as exact a80; normalized outputs match byte for
  byte. The unchanged leaf-m10 verifier also fails at the same assertion on a80.
  Both are linked from the existing baseline-CI problem, not bypassed.
- The current four-page result PDF and six one-page session PDFs have been
  visually inspected. Pages1/4 were rechecked after the held-out point-count and
  prose corrections. Benchmark copies and PDF hashes were regenerated together.
- Chromium154 desktop,390px/320px mobile, dark-theme, fragment navigation,
  no-JavaScript and no-external-request assertions pass. The runtime's managed
  policy blocks `file:` navigation, preserved in the first browser logs. The
  retry uses CDP's in-memory document fixture with the same local HTML and every
  original assertion; it changes no browser policy. `presentation/browser-check.json`
  records `document_mode: inline` and the exact page hash. Current desktop and
  mobile screenshots were visually inspected, including the new evidence table.

All frozen source/spec/Cargo hashes revalidate against measured commit0cbc;
`frozen-source-check.log` also verifies the original measured binary digest.
Integration audits use a separately built current verifier and add new receipts;
the original audit/replay files and failed outputs remain unchanged.

The final integrated verifier SHA256 is
244b1ba1fc1a4cd2eb996a672ddf2e45551c5876aa144153dc8b70532c575e12.
All15 sessions pass its audits:336 follow-up and867 parent verified outcomes
replay exactly. integration-audits/receipt-index.jsonl indexes the new receipts.
The supplemental37 follow-up failures and70 parent failures were already
replayed with their exact original binaries; no unsuccessful record is omitted.

The integrated CLI/examples invocation passes27 tests:16 ecbench integration,
five large-prime and six example tests. The ecbench binary and two examples
have no unit tests and are compiled, not described as tested cases. The final
leads-example retry passes its one test after the needless-borrow lint fix.
The native catalog binaries pass13 tests, then verify2101 covers and import258
standard records, leaving seven unresolved source records explicit.
Warnings-as-errors clippy passes for the ecbench binary and all six touched
examples with Rust1.98. The first lint failure is retained; no lint was suppressed.
Exact commands are in INTEGRATION-COMMANDS.md.
