# Follow-up decision rule (written before run 2)

Run 1 (`run-20261017T0331Z`, 3 rounds) showed all three switches exact and the
combination faster than the previous defaults on every F5 workload, direct
pack clearly faster, but direct unpack and table reuse mixed: on
`f5_n20_m20_d4` the arm with either switch off was faster than `on` in some
form/thread cells (down to 0.868×, its round range wholly below 1), while the
same cells' `off` arm was slower than `on`. Run 2 (`scripts/run_followup.sh`)
repeats the arms `on`, `nounpack`, `noreuse`, `off` on the F5 cases, 15
rounds, Reduced and SelectiveEchelon output, 1 and 4 threads, and `on` /
`noreuse` on the matrix-F4 cases n ≤ 18, 10 rounds.

A switch stays on by default only if, in every (form, threads) cell,
1. its median arm/on ratio is ≥ 0.97 on `f5_n20_m20_d4` and `f5_n24_m24_d4`
   (turning it off is at most 3 % faster), and
2. its median arm/on ratio is ≥ 1.00 on `f5_n24_m24_d4`;
and, for table reuse, its median arm/on ratio is ≥ 0.97 on every matrix-F4
case with n ≥ 16. Otherwise the switch goes back to opt-in in this PR.
Direct pack is decided by run 1 (arm/on 1.007–2.17, no cell below 1).
