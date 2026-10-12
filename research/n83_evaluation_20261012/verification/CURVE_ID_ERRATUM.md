# Curve-identity metadata erratum, 2026-10-12

This addendum corrects one identifier in the new native archived replay
`receipt.json` and its `REPORT.md`/`REPORT.pdf`. Those original artifacts are
retained unchanged. Their numerical field and curve inputs were correct;
the manually inserted `curve_model_slug` was copied from a registry entry
for a different field modulus.

The checked fixture has characteristic 2, degree 83, polynomial low terms
`[0,2,4,7]`, and equation `y²+xy=x³+x²+1`.
Its subgroup order is `8569786107849059`, cofactor `1128547018`,
generator `(6614244651398184860292015,8868260410029249391539529)` and
target `(2655898345537166180532960,3881526307991673265292615)`.

The correct explicit-field model slug, computed by the fresh strong
producer's identity path, is:

```
icv1-f2m83-t6151469093347-9eba6cd5
```

Its full model identity is:

```
ICV1:f2m-83-a5fd0d39:6151469093347:
9671406556910881928556062:0x1:-7:unk:r:9eba6cd5d359
```

The line break above is for rendering only; the full identity has no newline.
The incorrect historical alias `icv1-f2m83-t6151469093347-cdcc5432`
belongs to a model represented using low terms `[0,1,2,45]`. It must not
join these observations merely through the shared field degree and trace.

The archived checker used explicit polynomial `[0,2,4,7]`, exact retained
G and Q and the independently verified subgroup. Its scalar/base/witness
arithmetic result is unaffected. Receipt SHA-256:
`7d4e290524ea9b83ff75d48d20b17449fb60a6991bae787fdfee53ef28d461d5`.

The prepared `check_fresh.cpp` verifier checks explicit polynomial, subgroup,
G and Q. It does not depend on either slug and imports no registry entry.
No timed workload, frozen protocol, producer source or original receipt
was changed by this addendum; no checker execution ran during the panel.

Source context: fresh explicit-field identity in
`source_audit/n83_strong_point_eval.rs`; retained inputs in
`crypto:experiments/koblitz-single-target-n83-20261006/frozen/fixture.json`;
the conflicting historical alias in `crypto:docs/curves/registry.json`.
