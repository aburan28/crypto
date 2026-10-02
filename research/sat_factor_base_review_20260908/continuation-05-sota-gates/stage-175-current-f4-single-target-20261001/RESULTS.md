# Stage 175: current repository F4 on one frozen target

The current repository `BlockTables` F4 engine was run three times on the exact
opened `n=59, ell=9, m=3` true-negative used by Stage 174. All repeats
authenticated the source, visited all 512 fixed-X1 masks, completed all 242
rational systems, reproduced equation fingerprint `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`, and
returned exhaustive UNSAT. The algebraic factor-base contract still records no
target-subgroup enumeration and no discrete-log labels.

| backend | wall median (s) | core median (s) | peak RSS median (bytes) | status |
|---|---:|---:|---:|---|
| Stage 174 selected F4 | 26.358218 | 147.841771 | 2756362240 | historical exhaustive UNSAT |
| current repository F4 | 34.680670 | 267.171734 | 4030480384 | exhaustive UNSAT |
| same-binary direct MITM | 0.361851 | 0.351335 | 45072384 | exhaustive UNSAT |

The current F4 is `1.316x` the Stage 174 wall,
`1.807x` its CPU, and
`1.462x` its RSS. It is therefore rejected by the
preregistered rule. Against same-binary direct MITM it is
`95.84x` wall, `760.45x`
CPU, and `89.42x` RSS.

The repeated current-F4 structural counts are exact: 242 calls, 1,204 steps,
`319313687585` row-equivalent word XORs, `147794583858` actually performed
table-assisted word XORs, and no conflicts (not applicable to F4). The clean
build cost 67.121431 wall seconds,
278.361012 core-seconds, and
3460104192 bytes peak RSS. The complete
seven-process Stage 175 charge is 173.937518 wall
seconds, 1081.228415 core-seconds, and
4863803392 bytes maximum RSS.

This is a negative single-target engineering result. It does not change any
SOTA gate and does not measure natural relation yield or a full unknown-scalar
index-calculus run.
