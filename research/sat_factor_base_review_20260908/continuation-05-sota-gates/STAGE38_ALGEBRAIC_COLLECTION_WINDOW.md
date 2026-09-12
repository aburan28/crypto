# Stage 38: algebraic relation-collection window sweep

Stage 38 integrates the windowed, walked relation collector developed in draft PR #204 with the current Stage-35 public-target and evidence boundary. The upstream draft measured pseudo-random subgroup bases; this stage tests only the algebraically discovered Frobenius-union/two-torsion-saturation base.

A full three-summand pair-table scan finds a triple once with each summand in the scanned position, while retaining only one sorted witness. A cyclic window keeps any witness whose first, second, or third index falls in the window, then verifies the unsorted witness in the group. Additional probes are walked in deterministic runs of 64, so each run pays one scalar multiplication and shared additions. Work-unit boundaries still reproduce the same relation sequence.

The integration retains Stage 35's explicit `descent_summands` report and adds the effective `collection_window`. It also carries `summands_scanned` in every relation-unit document and reports initial and final totals. Older units deserialize the new counter as zero; new measurement admission requires the exact nonzero charge `trials × window`, or `trials × |F|` for the full-scan control.

The hosted sweep runs six cells on one inherited Linux CPU from one clean binary:

- full scan;
- windows 595, 297, 149, 74, and 37;
- 8,192 probes per work unit, one initial unit, and at most 64 units;
- the exact Stage-35 algebraic search, five public targets, two-summand descent, and signed-Frobenius rho in every cell.

The winner minimizes `collection + logs`, the two stages the window changes. Every cell still charges the repeated factor-base search and five-target descent, but their timing noise cannot choose this parameter. Selection occurs only after all six cells finish and does not change their public target stream or factor-base candidates. A later holdout must confirm any chosen window.

A local same-binary feasibility sweep found all six cells complete. Relative to the new full-scan control, windows reduced the final third-summand lookup charge while keeping 29 certified columns. Window 595 gave the lowest local precompute, 3.661 seconds versus 4.373 seconds full scan; window 149 gave 3.632 seconds in a consecutive run. These timings select no claim: the hosted sweep determines the candidate and a separate run must validate it.

This is a setup constant optimization. It does not replace the SAT matrix, change the algebraic factor-base requirement, supply Magma or external review, change the asymptotic exponent, or establish a Koblitz index-calculus SOTA result.
