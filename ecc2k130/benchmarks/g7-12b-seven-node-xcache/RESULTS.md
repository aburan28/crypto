# G7 seven-node shared-X cache result

No candidate passed the preregistered speedup gate. The selected executable is
retained.

| Variant | Median B updates/s | Paired ratio | 95% CI | Rate / 12 B/s | Decision |
|---|---:|---:|---:|---:|---|
| selected | 6.262430 | 1.000000 | [1.000000, 1.000000] | 0.521869 | reference |
| xcache-4 | 6.243798 | 0.995699 | [0.987072, 1.004402] | 0.520317 | unconfirmed |
| xcache-5 | 6.295066 | 1.000250 | [0.989535, 1.011081] | 0.524589 | unconfirmed |
| xcache-6 | 6.374050 | 1.009801 | [0.991428, 1.028515] | 0.531171 | unconfirmed |

All twelve timing samples completed on the same g7.2xlarge / RTX PRO 4500 at
165 W, with 34,359,738,368 complete scalar updates per sample. Correctness,
replay, resume, memcheck and shared-memory initcheck passed. Generic work is
`sqrt(n/262)`, ratio 1; full-DLP S is null.
