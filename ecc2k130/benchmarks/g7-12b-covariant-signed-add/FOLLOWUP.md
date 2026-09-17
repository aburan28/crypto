# Preregistered follow-up: decorrelated signed covariant addends

Status: closed at the 100-trial CPU gate before CUDA. The candidate averaged
767.390 complete iterations versus 157.500 selected, ratio 4.872317, with 61
overdue restarts.

The single-orbit candidate retains too much algebraic structure. Keep the
exhaustively proven phase and orientation functions, but give each of the eight
Hamming-weight branches its own base point

```
A_j = r_j B + t_j Q
r = [1,3,5,7,11,13,17,19]
t = [2,5,11,17,23,29,31,37].
```

Branch j applies plus or minus `sigma^phase(A_j)` according to the orientation
bit. Branch choice, phase and sign have the same Frobenius/negation behavior as
the first candidate, so the map remains covariant under both quotient actions.
Replay adds or subtracts `r_j*s^phase` and `t_j*s^phase` from its coefficients.
The decorrelated bases should break the witnessed short cycles.

Run the same 100 matched GF(2^23) DLP trials with restart1024 and complete work
accounting. Require all scalars and coefficient states, zero unexplained
infinities, and candidate mean work no more than 1.10x selected before a 2,000
trial study. CUDA is permitted only after that larger study. A CUDA design must
first measure global read-only versus compact shared storage for the 1,048
branch/phase points; no throughput gain is assumed. Promotion and 12 B/s gates
remain those in README.md.


## Result

All 100 candidate scalars and coefficient states recovered, but decorrelation
only reduced the single-orbit failure: 61 walks still crossed restart1024 and
mean complete work was **4.872317x** selected. This exceeds 1.10, so no
2,000-trial or CUDA work follows.
