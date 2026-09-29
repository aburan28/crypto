# Frozen results: no consistent descending-neighbor advantage

All 1,280 systems passed exact F4/F5 row-space, Boolean ideal, root and signed
point-lifting checks; three repetitions give 7,680 verified solver executions.
All four isogeny/dual certificates passed. No algebraic roots were rejected by
lifting in this corpus. The frozen success condition is **false**.

Across 1,024 neighbor/source pairs, ideal-completion degree fell in 133, stayed
equal in 775 and rose in 116. Both matrix variants have exactly the same degree
on each system. Ordered budgets range 4–8; canonical budgets range 4–7.
The canonical constraint lowers the maximum budget here, but that does not
establish a total solver or full-DLP improvement.

Neighbor F5 matrix cost/source ratios range from 0.839 to 1.133 across split and
encoding aggregates. Gains on one support need not persist on the holdout or
other encoding. The four neighbors are Frobenius-related, not four independent
statistical observations. Class: **accounting**, with a negative result for the
predeclared uniform-neighbor criterion. No full-DLP speedup or exponent is known.

## Matrix-stage ledger

One unit: <=64-bit row XOR, elimination plus F5 criterion over every degree
through completion. Counts are per repetition, summed over 64 systems in each
row (32 targets, two decomposition sizes). Completion-certificate costs and
non-matrix work are separate in raw evidence, not hidden inside this unit.
Ratios use the same solver, split and encoding on the source. Full-DLP S,
rho/floor ratios and speedup are null for every row, including the reference.
Correctness: every row 64/64 verified; classification: accounting.

| Split | Encoding | Model | Matrix variant | Degree range | Matrix XORs | Criterion XORs (included) | Cost/source |
|---|---|---|---|---:|---:|---:|---:|
| frozen | canonical | neighbor_224 | f4 | 4–7 | 148915 | 0 | 0.996000 |
| frozen | canonical | neighbor_224 | f5 | 4–7 | 148718 | 112 | 0.995102 |
| frozen | canonical | neighbor_225 | f4 | 4–7 | 127815 | 0 | 0.854875 |
| frozen | canonical | neighbor_225 | f5 | 4–7 | 127359 | 567 | 0.852185 |
| frozen | canonical | neighbor_92 | f4 | 4–7 | 150620 | 0 | 1.007404 |
| frozen | canonical | neighbor_92 | f5 | 4–7 | 150154 | 286 | 1.004711 |
| frozen | canonical | neighbor_93 | f4 | 4–7 | 138967 | 0 | 0.929464 |
| frozen | canonical | neighbor_93 | f5 | 4–7 | 138842 | 382 | 0.929020 |
| frozen | canonical | source | f4 | 4–7 | 149513 | 0 | 1.000000 |
| frozen | canonical | source | f5 | 4–7 | 149450 | 439 | 1.000000 |
| frozen | ordered | neighbor_224 | f4 | 4–8 | 168543 | 0 | 1.078654 |
| frozen | ordered | neighbor_224 | f5 | 4–8 | 168302 | 205 | 1.079343 |
| frozen | ordered | neighbor_225 | f4 | 4–8 | 153385 | 0 | 0.981645 |
| frozen | ordered | neighbor_225 | f5 | 4–8 | 152576 | 789 | 0.978490 |
| frozen | ordered | neighbor_92 | f4 | 4–8 | 177279 | 0 | 1.134564 |
| frozen | ordered | neighbor_92 | f5 | 4–8 | 176675 | 317 | 1.133040 |
| frozen | ordered | neighbor_93 | f4 | 4–8 | 160471 | 0 | 1.026995 |
| frozen | ordered | neighbor_93 | f5 | 4–8 | 160170 | 726 | 1.027192 |
| frozen | ordered | source | f4 | 4–8 | 156253 | 0 | 1.000000 |
| frozen | ordered | source | f5 | 4–8 | 155930 | 594 | 1.000000 |
| holdout | canonical | neighbor_224 | f4 | 4–7 | 159486 | 0 | 1.075653 |
| holdout | canonical | neighbor_224 | f5 | 4–7 | 159072 | 436 | 1.080520 |
| holdout | canonical | neighbor_225 | f4 | 4–7 | 143040 | 0 | 0.964733 |
| holdout | canonical | neighbor_225 | f5 | 4–7 | 142426 | 1215 | 0.967450 |
| holdout | canonical | neighbor_92 | f4 | 4–7 | 141676 | 0 | 0.955534 |
| holdout | canonical | neighbor_92 | f5 | 4–7 | 141161 | 679 | 0.958857 |
| holdout | canonical | neighbor_93 | f4 | 4–7 | 146447 | 0 | 0.987712 |
| holdout | canonical | neighbor_93 | f5 | 4–7 | 145840 | 26 | 0.990640 |
| holdout | canonical | source | f4 | 4–7 | 148269 | 0 | 1.000000 |
| holdout | canonical | source | f5 | 4–7 | 147218 | 482 | 1.000000 |
| holdout | ordered | neighbor_224 | f4 | 4–8 | 140862 | 0 | 0.832139 |
| holdout | ordered | neighbor_224 | f5 | 4–8 | 140443 | 500 | 0.839222 |
| holdout | ordered | neighbor_225 | f4 | 4–8 | 161450 | 0 | 0.953762 |
| holdout | ordered | neighbor_225 | f5 | 4–8 | 160697 | 1139 | 0.960251 |
| holdout | ordered | neighbor_92 | f4 | 4–8 | 178330 | 0 | 1.053480 |
| holdout | ordered | neighbor_92 | f5 | 4–8 | 177693 | 1420 | 1.061811 |
| holdout | ordered | neighbor_93 | f4 | 4–8 | 149907 | 0 | 0.885572 |
| holdout | ordered | neighbor_93 | f5 | 4–8 | 149060 | 144 | 0.890713 |
| holdout | ordered | source | f4 | 4–8 | 169277 | 0 | 1.000000 |
| holdout | ordered | source | f5 | 4–8 | 167349 | 1976 | 1.000000 |

## Paired degree changes

| Neighbor b | Lower | Equal | Higher |
|---|---:|---:|---:|
| 93 | 28 | 202 | 26 |
| 225 | 30 | 196 | 30 |
| 92 | 36 | 188 | 32 |
| 224 | 39 | 189 | 28 |

## Evidence and limitations

- `results/raw.json.gz`: all raw systems, traces, roots, certificate and phase counters.
- `results/summary.json`: frozen aggregates and the falsification result.
- `contract.json`: frozen inputs, metric, bounds and decision rule.
- `README.md`: map formulas, conductor argument, solver provenance and reproduction.

This measures a bounded Python reference of the repository's Boolean matrix-F5
criterion. It does not execute native Rust or a full signature-based F5 solver.
Exhaustive ANF interpolation and root extraction limit the experiment to tiny
systems; the reported degree is ideal-completion budget, not intrinsic d_reg.
Only transported supports are compared; native-neighbor support selection and
invariant-coordinate elimination remain outside this measured result. There is
no inference about production curves or end-to-end cryptanalysis performance.
