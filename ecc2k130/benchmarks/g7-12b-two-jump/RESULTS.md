# G7 two-jump rho result

Neither candidate passed the preregistered speedup gate. The selected
 eight-jump executable is retained.

| Variant | Median B updates/s | Paired ratio | 95% CI | Collision work ratio | Adjusted B/s | Decision |
|---|---:|---:|---:|---:|---:|---|
| selected | 6.279988 | 1.000000 | [1.000000, 1.000000] | 1.000000 | 6.279988 | reference |
| two-jump indexed | 6.273318 | 0.992430 | [0.964207, 1.021479] | 1.013080 | 6.192320 | unconfirmed |
| two-jump immediate | 5.418604 | 0.844542 | [0.778245, 0.916487] | 1.013080 | 5.348642 | regression |

All 6,000 planted DLPs in the matched collision cohort were recovered. The
fixed sigma helpers, full states, DP records, CPU replay, checkpoint resume,
memcheck, initcheck and synccheck passed. All nine timing samples performed
34,359,738,368 complete scalar updates on the same g7.2xlarge / RTX PRO 4500
at 165 W. Generic work is sqrt(n/262), the measured two-jump work ratio is
1.013080, and full-DLP S is null.
