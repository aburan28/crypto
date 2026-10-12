# `(10, 5)` u0, first attempt: stopped by hand, not evidence

Binary `0227cd61` with `KIC_SPARSE_DENSE_BUDGET_MB=11000`, started
2026-09-30T04:13:05Z and stopped (PID kill) after 4 h 26 min.  Peak resident
memory 10.06 GB (`cell-10-5-7.u0.hwm`).  It wrote no outcome.

- The survivors after band 7 did not fit the 11,000 MB budget: band 7 is
  rank-deficient, so more than about 376k of the 836,820 rows survived it.
  The engine therefore went on eliminating band 6 sparsely.
- Two gdb samples of the loop's column counter (callee-saved registers of
  `sparse_then_dense`), about a minute apart, read 495,848 and 495,900:
  inside band 6 (columns 480,700 to 657,800), advancing at under one column
  a second, so about two days to the end of the band.
- It was stopped for the F5-criterion rows (the amendment in
  PREREGISTRATION.md), which drop 171,003 of the 836,820 rows on this draw
  and so the same number of survivors.

`cell-10-5-7.u0.phase-sample.txt` is every tenth minute of the RSS/CPU log.
A stopped run is a resource limit, not a mathematical result.
