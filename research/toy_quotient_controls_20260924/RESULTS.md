# Results: quotient benefit depends on decomposition size

The finite quotient helps the three-summand Python reference, but regresses for
two summands. On fresh supports, invariant/ordered complete-call opcode ratios
are 0.257–0.307 for three summands and
1.060–1.316 for two summands across the five models.
The predeclared requirement to improve every model and both summand counts by
at least 20% is **false**. The three-summand result
is a bounded Python opcode-proxy improvement, not a native runtime or full-DLP
speedup. Classification: accounting; full-DLP S and rho/floor ratios remain null.

All 1920/1920 systems preserve the exact unordered signed-point
decomposition set. Three uninstrumented repetitions agree, and an additional
instrumented replay agrees on each input. Independent Buchberger audits verify
240/240 selected systems. Actual
signature-based SymPy F5B completes 101/240;
139 exhaust the declared work budget. Those are explicitly
incomplete F5B runs, not successful F5B validations.

## Representation effects

Coordinate rescaling alone lowers the matrix ideal-completion budget in
64 cases, leaves it equal in 178, and raises it in
78. Thus degree changes of this kind are not unique to isogeny
descent. Selector relabeling and equation reversal preserve every completion
degree while changing cost, as expected from their unchanged degree filtration.

Within the invariant representation, neighbor/source degrees are lower in
0, equal in 256, and higher in 0 comparisons.
This finite-image encoding does not establish that descent makes the algebra
intrinsically easier. Its rank labeling and exhaustive preprocessing remain
representation choices, fully described in README.md.

## Complete-call proxy ledger

Unit: executed CPython opcode dispatches through setup, encoding, matrix solve
and root extraction, reconstruction, and point verification. Each row aggregates
40 inputs (five curve models, eight targets). Ratios use ordered on the same
split and summand count. This excludes C-internal instructions, arithmetic-table
and fixture generation, and independent audit work; no calibrated machine-operation
total is available. Per-model totals and exclusive phase counts are in the raw
and summary artifacts. Every row is classified as accounting.

| Split | Summands | Variant | Python opcodes | Cost/ordered | Verified |
|---|---:|---|---:|---:|---:|
| fresh_a | 2 | canonical | 5246020 | 0.956401 | 40/40 |
| fresh_a | 2 | invariant | 6260311 | 1.141316 | 40/40 |
| fresh_a | 2 | ordered | 5485169 | 1.000000 | 40/40 |
| fresh_a | 2 | relabel | 5493151 | 1.001455 | 40/40 |
| fresh_a | 2 | reverse_equations | 5479988 | 0.999055 | 40/40 |
| fresh_a | 2 | scale | 5362363 | 0.977611 | 40/40 |
| fresh_a | 3 | canonical | 54923758 | 0.778149 | 40/40 |
| fresh_a | 3 | invariant | 19093504 | 0.270513 | 40/40 |
| fresh_a | 3 | ordered | 70582572 | 1.000000 | 40/40 |
| fresh_a | 3 | relabel | 70222527 | 0.994899 | 40/40 |
| fresh_a | 3 | reverse_equations | 70463470 | 0.998313 | 40/40 |
| fresh_a | 3 | scale | 70475959 | 0.998490 | 40/40 |
| fresh_b | 2 | canonical | 5325716 | 0.983704 | 40/40 |
| fresh_b | 2 | invariant | 6378400 | 1.178144 | 40/40 |
| fresh_b | 2 | ordered | 5413940 | 1.000000 | 40/40 |
| fresh_b | 2 | relabel | 5428569 | 1.002702 | 40/40 |
| fresh_b | 2 | reverse_equations | 5419652 | 1.001055 | 40/40 |
| fresh_b | 2 | scale | 5732985 | 1.058930 | 40/40 |
| fresh_b | 3 | canonical | 62268227 | 0.900789 | 40/40 |
| fresh_b | 3 | invariant | 19174073 | 0.277377 | 40/40 |
| fresh_b | 3 | ordered | 69126293 | 1.000000 | 40/40 |
| fresh_b | 3 | relabel | 69799697 | 1.009742 | 40/40 |
| fresh_b | 3 | reverse_equations | 69067074 | 0.999143 | 40/40 |
| fresh_b | 3 | scale | 68683671 | 0.993597 | 40/40 |
| frozen | 2 | canonical | 5222417 | 0.927072 | 40/40 |
| frozen | 2 | invariant | 6435392 | 1.142396 | 40/40 |
| frozen | 2 | ordered | 5633240 | 1.000000 | 40/40 |
| frozen | 2 | relabel | 5665736 | 1.005769 | 40/40 |
| frozen | 2 | reverse_equations | 5629381 | 0.999315 | 40/40 |
| frozen | 2 | scale | 5721240 | 1.015622 | 40/40 |
| frozen | 3 | canonical | 52933993 | 0.781225 | 40/40 |
| frozen | 3 | invariant | 18804886 | 0.277531 | 40/40 |
| frozen | 3 | ordered | 67757690 | 1.000000 | 40/40 |
| frozen | 3 | relabel | 66631636 | 0.983381 | 40/40 |
| frozen | 3 | reverse_equations | 67732189 | 0.999624 | 40/40 |
| frozen | 3 | scale | 73160461 | 1.079737 | 40/40 |
| previous_holdout | 2 | canonical | 5271525 | 0.950255 | 40/40 |
| previous_holdout | 2 | invariant | 6253184 | 1.127210 | 40/40 |
| previous_holdout | 2 | ordered | 5547487 | 1.000000 | 40/40 |
| previous_holdout | 2 | relabel | 5525486 | 0.996034 | 40/40 |
| previous_holdout | 2 | reverse_equations | 5528986 | 0.996665 | 40/40 |
| previous_holdout | 2 | scale | 5812556 | 1.047782 | 40/40 |
| previous_holdout | 3 | canonical | 56362851 | 0.795766 | 40/40 |
| previous_holdout | 3 | invariant | 19468944 | 0.274875 | 40/40 |
| previous_holdout | 3 | ordered | 70828380 | 1.000000 | 40/40 |
| previous_holdout | 3 | relabel | 69603641 | 0.982708 | 40/40 |
| previous_holdout | 3 | reverse_equations | 70845466 | 1.000241 | 40/40 |
| previous_holdout | 3 | scale | 75358926 | 1.063965 | 40/40 |

## Fresh quotient comparisons by model

| Split | Model | Summands | Opcode cost/ordered |
|---|---|---:|---:|
| fresh_a | neighbor_224 | 2 | 1.071425 |
| fresh_a | neighbor_224 | 3 | 0.256862 |
| fresh_a | neighbor_225 | 2 | 1.315739 |
| fresh_a | neighbor_225 | 3 | 0.303982 |
| fresh_a | neighbor_92 | 2 | 1.136796 |
| fresh_a | neighbor_92 | 3 | 0.273342 |
| fresh_a | neighbor_93 | 2 | 1.156150 |
| fresh_a | neighbor_93 | 3 | 0.258618 |
| fresh_a | source | 2 | 1.060227 |
| fresh_a | source | 3 | 0.264359 |
| fresh_b | neighbor_224 | 2 | 1.155073 |
| fresh_b | neighbor_224 | 3 | 0.267461 |
| fresh_b | neighbor_225 | 2 | 1.168024 |
| fresh_b | neighbor_225 | 3 | 0.259446 |
| fresh_b | neighbor_92 | 2 | 1.231401 |
| fresh_b | neighbor_92 | 3 | 0.274934 |
| fresh_b | neighbor_93 | 2 | 1.177291 |
| fresh_b | neighbor_93 | 3 | 0.281090 |
| fresh_b | source | 2 | 1.162412 |
| fresh_b | source | 3 | 0.306950 |

## Decision and reproducibility

Retain the finite-quotient implementation as a tested toy representation with
an explicit setup/lifting cost. Reject the blanket improvement claim. A future
three-summand-only hypothesis would require a new predeclared protocol and new
holdouts; the present two-summand regressions must remain visible.

- `results/run-002/raw.json.gz`: canonical raw inputs, phase counts, per-degree
  traces, repetitions, exact-set hashes, independent bases and budget receipts.
- `results/run-002/summary.json`: full per-model aggregates and decisions.
- `results/run-001/`: superseded audit instrumentation run, retained with its
  original independent-audit source and contract. It is not accepted F5B evidence.
- `contract.json`, `protocol-v1.json`, `protocol-v2.json`: protocol and revision history.
- `dependency-replay.json`: exact replay receipt for the preceding full corpus.
- `verify_replay.py`: hashes, exclusive accounting, correctness and exact replay gates.

The independent observer reports only degrees seen at signature-reduction
boundaries and final reduced-basis degree. Neither these observations nor the
Boolean matrix ideal-completion budget are intrinsic degree of regularity.
