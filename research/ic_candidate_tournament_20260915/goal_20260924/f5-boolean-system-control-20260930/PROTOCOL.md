# Disclosed F5 Boolean-system correctness control

Status at registration: pending one execution. This is a small n17 kb1
correctness diagnosis of two previously exposed ordinary queries, not a new
complete-DLP candidate, fresh workload, yield estimate or speed comparison.
The three confirmation rounds and every historical registration remain closed.

The preceding retrospective audit in PR #1046 found finite mathematical S3
chains for ordinary queries 164 and 173 but explicitly left the implemented
Boolean encoding, parameter specialization, variable permutation and F5 row
consequences unproved. This control tests those properties on exactly these
two disclosed inputs. It does not select a budget, order or search policy.

Hypothesis: their full-point chains satisfy the exact implemented ANF; both
the unfiltered F4 and F5-criterion root reductions span the independently
constructed degree-three Macaulay rows on each registered layout.

`inputs.json` fixes the original n17/a1 polynomial-basis curve, subgroup,
dimension-six standard base and two point witnesses from the closed archive
`039e535fa73e5c794aeff666bc02989076f04abbf751b5f609884a2785760678`.
The original source is commit `765c3c5f19032bd852163805f257c56babef2040`,
whose complete root/dependency manifest is
`c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7`.
The runner rejects a different source or dependency manifest. It copies those
inputs to a new build tree and adds only `examples/f5_boolean_control.rs` from
`export.rs`. It never modifies or executes the historical worker. This
diagnostic overlay has its own compiled-source manifest and native binary
registration, retained before native execution; it is not an old candidate ID.

For each query, export all coefficients of the direct constructor, explicitly
instantiated template, reuse path, repeated reuse path and interleaved system.
Independently regenerate both field-valued Boolean S3 polynomials using Python
polynomial-basis arithmetic and square-free monomials. Compare every one of
the 34 coordinate polynomials exactly. This establishes coefficient identity
on the entire 35-variable Boolean cube at the registered parameters, beyond
the twelve disclosed witness assignments. Independently derive and evaluate
all six full-point chain assignments per query, including their renamed
assignments; separately replay their public ordinary-query scalars.

Build native full-readback root matrices at degree three for original and
interleaved layouts, each with `RowCriterion::None` and `RowCriterion::F5`.
These are eight bounded reductions, not eight PDP searches. Independently
enumerate all products with multiplier degree at most `3 - deg(generator)`;
XOR equal Boolean monomials, construct the full row space with Python integers,
and compare both native spans and every disclosed assignment. Exact row-space
equality is the success criterion. A mismatch, oversize, error or timeout is
a terminal negative/incomplete control, retained in a new output directory.

Run the native control once, one Rayon thread, no ambient solver/cache/kernel
overrides, 180 seconds for native execution and 1800 seconds for compilation.
All argv, source/dependency hashes and source archives, compiler/native-tool
receipts, environment, input hashes, raw stdout/stderr, native binary and
pre/post source gates remain available. The output directory refuses reuse.
No verified online time, common operation total or speedup is reported.
Native criterion/elimination word XORs remain bounded matrix diagnostics;
they omit construction/readback and are not a full cost measure.

Passing establishes the registered ANF and full-readback root-row properties
only. The production linear-tail path, specialized branch matrices, cache
observations, traversal, budget policy, lifting completeness, natural yield
and complete one-target DLP admission remain separate gates. The interleaved
layout is an extra correctness control: it was not the historical default
MatrixF5/LowestFree traversal. No claim extends to another field/base/parameter.
Fresh paired qualification and a calibrated rho comparison require those
remaining gates; this control consumes no confirmation attempt.

Before execution, commit this protocol, `protocol.json`, exact inputs, exporter,
runner, independent auditor and fault controls. Then run:

```sh
env PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/goal_20260924/f5-boolean-system-control-20260930/run_control.py \
  --historical-source /absolute/immutable-765c3c5f19032bd852163805f257c56babef2040-checkout \
  --out /absolute/new-f5-boolean-control-output
```

The result PR retains the unmodified execution artifacts and terminal verdict.
Post-execution analysis corrections preserve the first receipt and raw output.
