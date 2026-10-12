# Fixed-width reduction of the K_0/K_1 n=83 field polynomial

Registered before baseline on 2026-10-04. The n=83 field uses
`z^83 + z^45 + z^2 + z + 1`. Generic `F2mElement` multiplication and
squaring currently reduce it through the general chunk routine because
`83 - 45 = 38 < 64`. Test an exact specialized reducer that folds the
high polynomial coefficients through this fixed modulus using at most three
`u128` shift-XOR folds. Apply it only when the field degree and complete
low-term list match; all other fields retain the existing path.

Before editing source, run the release `f6_wide_n83_index_probe` on the
actual K_0 dimension-8 and dimension-10 subgroup bases and public T001,
and run `f6_wide_n83_full_index_probe` on the complete dimension-12 base.
Use their existing warmup/repetition policy, exact pair caps, and in-process
peak RSS reading. Then rerun unchanged probes after the reduction change
with the same compiler and host. Save all rows, outcomes, failures, source
hashes, run order, and exact target. Do not mix compilation or base setup
into the reported build/query times.

Correctness gate: compare specialized multiplication and squaring against
the independent schoolbook and bitwise-reduction references on fixed n=83
vectors including high bits, zero, and random-looking operands; also pass
the n=83 pair-closure group tests and preserve all benchmark outcomes.
Retain the optimization only if exploratory median build and query times
at both smaller bases, and the single full-base build and query times, show
no greater than 10% regression. A speedup claim needs an isolation receipt;
these runs are unisolated stage diagnostics. A complete n=83 F6 or IC/rho
speedup remains unknown regardless of this reducer's result.
