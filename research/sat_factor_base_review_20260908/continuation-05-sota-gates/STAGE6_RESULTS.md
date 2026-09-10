# Stage 6: degree-15 undersized-base control

The frozen first extension beyond degree 9 used the single-factor kernel `single:0`, with dimension 4, 31 rational points, and two signed Frobenius orbits. The intended two-summand system had only `m*ell = 8` factor coordinates for a degree-15 target.

The original run was invalid as an index-calculus completion. After 151 SAT refutations and zero factor-base relations, trial 152 reached `aG+bQ=O` and recovered scalar 101 directly. The amendment retained this run, disabled that bypass, and charged every later direct event as a skipped trial.

The corrected successor exhausted 20,000 trials. It made 19,910 SAT calls, recorded 19,910 terminal refutations and 13,040 conflicts, skipped 90 direct events, found no factor-base relation, and did not recover the scalar. It used 2.763971 core-seconds, 4.718951 seconds wall, and 2.97 MiB peak RSS. The matched signed-Frobenius rho run recovered and verified scalar 101 in two iterations and 23 reported group additions, using 0.005979 core-seconds, 0.008396 seconds wall, and 2.03 MiB peak RSS.

This is a scoped negative result for one dimension-4 kernel with `m=2`. It motivated the target-independent dimension rule used in Stage 7; it is not evidence against other algebraic factor bases or summand counts.
