# Cached fixed spans

Engineering hypothesis, declared before measurements: cache the prefix-independent part of each block relaxation, preserving every accept/reject decision. Use the coefficient-reuse and batch-inversion solver as the unchanged reference and direct S3 as the strongest measured competitor. Chained S3 and symmetrized S4 remain in the cold screen. This is an equivalent family-specific solver-stage suite; the WDSat protocol uses a different curve/interface. No full-DLP gain is implied.

At fixed function coefficient b and k free low bits of z, write the exact circuit as

\[
C'+\sum_{i<k} Z_i z_i+\sum_j W'_j w_j+\sum_{i<k,j}M_{ij}z_iw_j=0.
\]

Let T_k be the binary span of Z_i and M_ij with i<k. It is independent of the fixed high-bit prefix and T_k is contained in T_(k+1). The old relaxation accepts precisely when C' belongs to T_k + span(W'_j). Precomputing a basis of T_k, copying it, then inserting W'_j therefore produces exactly the same subspace. If T_k already has ambient rank, rejection is impossible. Any dual separator of the resulting basis annihilates every original relaxed column and separates C', proving every rejection. The inherited independent checker derives columns by finite differences of the unexpanded identity; it does not trust this basis construction.

All Gaussian additions, basis copies, table construction, failed targets, point extraction and signed-sum verification are charged. Field API and selected binary word counts remain separate. Logical field-slot storage is not peak RSS. Reordering Gaussian pivots may change the separator but must preserve the rejection blocks, leaf visits and complete outputs.

The reference boundary is the strongest direct S3 table on identical targets. The unchanged coefficient-preparation cost is an architecture-specific floor for this candidate; no universal operation floor or calibration is claimed. The signed-triple yield ceiling is unchanged. Full-DLP S, rho, calibrated total-cost and floor ratios remain null.

The screening goal is at least 20% fewer field API calls than reuse on a complete equal-work batch, with all certificates valid. It is not a Semaev win unless direct S3 is crossed in a measured common total-cost unit. Abandon this implementation if setup offsets the saved elimination or correctness fails; retain all regressions. The original three-size 20% calibrated goal remains the promotion target.
