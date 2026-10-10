
## research/ic_leads_20261010/runs/20261010T193217Z/report_n13.json

instance `icv1-f2m13-t181-515ee569`  log2 r = 11.0  status complete  rho reference: negation-map r-adding walk: look-ahead, short-cycle escape by doubling, distinguished points, stride starts A=2 mean S=2.67 (plain A=1: S=19.2)
calibration: pinned units ['ns_per_double', 'ns_per_lookup', 'ns_per_row_op', 'ns_per_as_solve', 'ns_per_word_xor', 'ns_per_frobenius', 'ns_per_canon'], measured on this host [] (GAE of a measured unit is host-dependent; the counted columns are not)
skipped: koblitz-orbit[divisor=1] + mitm-frobenius-counted + — + walk: oracle `mitm-frobenius-counted` is not available here
skipped: koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius-counted + — + walk: oracle `mitm-frobenius-counted` is not available here

### medians over repeats

| configuration | repeats | verified | exhausted | signed_pts | columns | pts/col | tried | found | trials/rel | rows | rank | la_work | canon | frob_maps | lookups | tbl_entries | fb_as_solves | setup_adds | decomp_adds | decomp_doubles | decomp_smults | la_row_ops | S | S/rho | wall_s(note) | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 4e+03 | 155 | 25.8 | 16 | 16 | 1 | 16 | 16 | 208 | 0 | 4e+03 | 16 | 0 | 4.1e+03 | 4.014e+06 | 200 | 308 | 34 | 208 | 8.97e+04 | 3.36e+04 | 3.69 | 480 | 4.014e+06 | 510 | 3.52 | 14 | 4.015e+06 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 4e+03 | 2e+03 | 2 | 99.5 | 97 | 1.02 | 97 | 97 | 1.18e+03 | 0 | 4e+03 | 97 | 0 | 4.1e+03 | 4.014e+06 | 392 | 341 | 38 | 1.18e+03 | 8.97e+04 | 3.36e+04 | 3.12 | 480 | 4.014e+06 | 724 | 20 | 14 | 4.015e+06 |

### ratio to baseline `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` (lead / baseline; the repeated baseline row is the noise floor)

| configuration | columns | pts/col | trials/rel | rows | la_row_ops | canon | setup_adds | decomp_adds | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 1 | 1 | 1 | 1 | 1 | — | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk` | 12.9 | 0.0774 | 1.02 | 6.06 | 5.67 | — | 1 | 1.96 | 1 | 1 | 1.42 | 5.67 | 1 | 1 |

### B·D²/N (B = signed points, D = online GAE per target = decomp + la + verify, N = r)

- `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk`: B=4e+03 D=528 → B·D²/N = 5.57e+05
- `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk`: B=4e+03 D=758 → B·D²/N = 1.149e+06

## research/ic_leads_20261010/runs/20261010T193624Z/report_n17.json

instance `icv1-f2m17-tm101-00378d4e`  log2 r = 16.0  status complete  rho reference: negation-map r-adding walk: look-ahead, short-cycle escape by doubling, distinguished points, stride starts A=2 mean S=1.65 (plain A=1: S=5.75)
calibration: pinned units ['ns_per_double', 'ns_per_lookup', 'ns_per_row_op', 'ns_per_as_solve', 'ns_per_word_xor', 'ns_per_frobenius', 'ns_per_canon'], measured on this host [] (GAE of a measured unit is host-dependent; the counted columns are not)

### medians over repeats

| configuration | repeats | verified | exhausted | signed_pts | columns | pts/col | tried | found | trials/rel | rows | rank | la_work | canon | frob_maps | lookups | tbl_entries | fb_as_solves | setup_adds | decomp_adds | decomp_doubles | decomp_smults | la_row_ops | S | S/rho | wall_s(note) | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` | 8 | 8 | 0 | 205 | 7 | 29.3 | 18.5 | 6 | 2.64 | 6 | 6 | 61 | 18.5 | 205 | 18.5 | 0 | 256 | 919 | 283 | 477 | 34 | 61 | 6.75 | 4.1 | 0.00115 | 26.2 | 919 | 761 | 0.842 | 21 | 1.73e+03 |
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 205 | 7 | 29.3 | 19.5 | 6 | 3 | 6 | 6 | 60 | 0 | 205 | 19.5 | 0 | 256 | 1.07e+04 | 286 | 477 | 34 | 60 | 45 | 27.3 | 0.0035 | 26.2 | 1.07e+04 | 760 | 0.828 | 21 | 1.15e+04 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius + — + walk` | 8 | 8 | 0 | 205 | 103 | 1.99 | 188 | 51 | 3.93 | 51 | 51 | 503 | 188 | 205 | 188 | 0 | 256 | 1.07e+04 | 507 | 477 | 34 | 503 | 45.9 | 27.9 | 0.0131 | 26.2 | 1.07e+04 | 1e+03 | 6.94 | 21 | 1.18e+04 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 205 | 103 | 1.99 | 215 | 55 | 3.87 | 55 | 55 | 631 | 0 | 205 | 215 | 0 | 256 | 1.07e+04 | 518 | 480 | 34 | 631 | 45.9 | 27.9 | 0.00501 | 26.2 | 1.07e+04 | 995 | 8.71 | 21 | 1.18e+04 |

### ratio to baseline `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` (lead / baseline; the repeated baseline row is the noise floor)

| configuration | columns | pts/col | trials/rel | rows | la_row_ops | canon | setup_adds | decomp_adds | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 1 | 1 | 1.14 | 1 | 0.984 | 0 | 11.7 | 1.01 | 1 | 11.7 | 0.998 | 0.984 | 1 | 6.66 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius + — + walk` | 14.7 | 0.068 | 1.49 | 8.5 | 8.25 | 10.2 | 11.7 | 1.79 | 1 | 11.7 | 1.32 | 8.25 | 1 | 6.81 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk` | 14.7 | 0.068 | 1.46 | 9.17 | 10.3 | 0 | 11.7 | 1.83 | 1 | 11.7 | 1.31 | 10.3 | 1 | 6.8 |

### B·D²/N (B = signed points, D = online GAE per target = decomp + la + verify, N = r)

- `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk`: B=205 D=783 → B·D²/N = 1.92e+03
- `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk`: B=205 D=781 → B·D²/N = 1.91e+03
- `koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius + — + walk`: B=205 D=1.03e+03 → B·D²/N = 3.31e+03
- `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk`: B=205 D=1.02e+03 → B·D²/N = 3.28e+03

## research/ic_leads_20261010/runs/20261010T193620Z/report_n23.json

instance `icv1-f2m23-t5197-69e76b73`  log2 r = 21.0  status complete  rho reference: signed-Frobenius r-adding walk, distinguished points (koblitz_signed_frobenius_rho_reference) A=46 mean S=1.06 (plain A=1: S=2.71)
calibration: pinned units [], measured on this host ['ns_per_double', 'ns_per_lookup', 'ns_per_row_op', 'ns_per_as_solve', 'ns_per_word_xor', 'ns_per_frobenius'] (GAE of a measured unit is host-dependent; the counted columns are not)

### medians over repeats

| configuration | repeats | verified | exhausted | signed_pts | columns | pts/col | tried | found | trials/rel | rows | rank | la_work | canon | frob_maps | lookups | tbl_entries | fb_as_solves | setup_adds | decomp_adds | decomp_doubles | decomp_smults | la_row_ops | S | S/rho | wall_s(note) | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` | 8 | 8 | 0 | 2.02e+03 | 45 | 45 | 62.5 | 27.5 | 2.61 | 27.5 | 27.5 | 268 | 62.5 | 2.02e+03 | 62.5 | 0 | 2.05e+03 | 4.76e+04 | 430 | 648 | 34 | 268 | 33.9 | 32.1 | 0.222 | 357 | 4.76e+04 | 1.08e+03 | 0.266 | 30 | 4.9e+04 |
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 2.02e+03 | 45 | 45 | 62 | 26.5 | 2.53 | 26.5 | 26.5 | 251 | 0 | 2.02e+03 | 62 | 0 | 2.05e+03 | 1.027e+06 | 432 | 648 | 34 | 251 | 711 | 673 | 8.95 | 357 | 1.027e+06 | 1.07e+03 | 0.248 | 30 | 1.029e+06 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius + — + walk` | 8 | 8 | 0 | 2.02e+03 | 1.01e+03 | 2 | 1.18e+03 | 466 | 2.62 | 466 | 466 | 4.92e+03 | 1.18e+03 | 2.02e+03 | 1.18e+03 | 0 | 2.05e+03 | 1.027e+06 | 2.01e+03 | 648 | 34 | 4.92e+03 | 712 | 674 | 4.17 | 357 | 1.027e+06 | 2.68e+03 | 4.86 | 30 | 1.030e+06 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 2.02e+03 | 1.01e+03 | 2 | 1.23e+03 | 468 | 2.6 | 468 | 468 | 5.14e+03 | 0 | 2.02e+03 | 1.23e+03 | 0 | 2.05e+03 | 1.027e+06 | 2.06e+03 | 652 | 34 | 5.14e+03 | 712 | 674 | 2.58 | 357 | 1.027e+06 | 2.72e+03 | 5.09 | 30 | 1.030e+06 |

### ratio to baseline `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` (lead / baseline; the repeated baseline row is the noise floor)

| configuration | columns | pts/col | trials/rel | rows | la_row_ops | canon | setup_adds | decomp_adds | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 1 | 1 | 0.969 | 0.964 | 0.935 | 0 | 21.6 | 1 | 1 | 21.6 | 0.997 | 0.935 | 1 | 21 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius + — + walk` | 22.5 | 0.0444 | 1 | 16.9 | 18.3 | 18.9 | 21.6 | 4.68 | 1 | 21.6 | 2.49 | 18.3 | 1 | 21 |
| `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk` | 22.5 | 0.0444 | 0.996 | 17 | 19.1 | 0 | 21.6 | 4.8 | 1 | 21.6 | 2.53 | 19.1 | 1 | 21 |

### B·D²/N (B = signed points, D = online GAE per target = decomp + la + verify, N = r)

- `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk`: B=2.02e+03 D=1.11e+03 → B·D²/N = 1.18e+03
- `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk`: B=2.02e+03 D=1.1e+03 → B·D²/N = 1.18e+03
- `koblitz-orbit[divisor=1,no_fold=1] + mitm-frobenius + — + walk`: B=2.02e+03 D=2.71e+03 → B·D²/N = 7.11e+03
- `koblitz-orbit[divisor=1,no_fold=1] + mitm[negation_folded=1] + — + walk`: B=2.02e+03 D=2.76e+03 → B·D²/N = 7.36e+03

## research/ic_leads_20261010/runs/20261010T193711Z/report_n23.json

instance `icv1-f2m23-t5197-69e76b73`  log2 r = 21.0  status complete  rho reference: signed-Frobenius r-adding walk, distinguished points (koblitz_signed_frobenius_rho_reference) A=46 mean S=1.06 (plain A=1: S=2.71)
calibration: pinned units [], measured on this host ['ns_per_double', 'ns_per_lookup', 'ns_per_row_op', 'ns_per_as_solve', 'ns_per_word_xor', 'ns_per_frobenius'] (GAE of a measured unit is host-dependent; the counted columns are not)
skipped: binary-subspace[dimension=11] + mitm-frobenius + — + walk: this instance has no Koblitz structure to fold by

### medians over repeats

| configuration | repeats | verified | exhausted | signed_pts | columns | pts/col | tried | found | trials/rel | rows | rank | la_work | canon | frob_maps | lookups | tbl_entries | fb_as_solves | setup_adds | decomp_adds | decomp_doubles | decomp_smults | la_row_ops | S | S/rho | wall_s(note) | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` | 8 | 8 | 0 | 2.02e+03 | 45 | 45 | 62.5 | 27.5 | 2.61 | 27.5 | 27.5 | 268 | 62.5 | 2.02e+03 | 62.5 | 0 | 2.05e+03 | 4.76e+04 | 430 | 648 | 34 | 268 | 33.8 | 32.1 | 0.168 | 299 | 4.76e+04 | 1.08e+03 | 10.4 | 30 | 4.9e+04 |
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 2.02e+03 | 45 | 45 | 62 | 26.5 | 2.53 | 26.5 | 26.5 | 251 | 0 | 2.02e+03 | 62 | 0 | 2.05e+03 | 1.027e+06 | 432 | 648 | 34 | 251 | 710 | 673 | 3.72 | 299 | 1.027e+06 | 1.07e+03 | 9.75 | 30 | 1.029e+06 |
| `binary-subspace[dimension=11] + mitm[negation_folded=1] + — + walk` | 8 | 8 | 0 | 2.08e+03 | 1.04e+03 | 2 | 1.55e+03 | 362 | 4.28 | 362 | 362 | 3.15e+03 | 0 | 0 | 1.55e+03 | 0 | 2.05e+03 | 1.081e+06 | 2.28e+03 | 648 | 34 | 3.15e+03 | 749 | 709 | 2.4 | 286 | 1.081e+06 | 2.94e+03 | 122 | 30 | 1.084e+06 |

### ratio to baseline `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` (lead / baseline; the repeated baseline row is the noise floor)

| configuration | columns | pts/col | trials/rel | rows | la_row_ops | canon | setup_adds | decomp_adds | gae_fb | gae_setup | gae_decomp | gae_la | gae_verify | gae_total |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk` | 1 | 1 | 0.969 | 0.964 | 0.935 | 0 | 21.6 | 1 | 1 | 21.6 | 0.997 | 0.935 | 1 | 21 |
| `binary-subspace[dimension=11] + mitm[negation_folded=1] + — + walk` | 23.1 | 0.0444 | 1.64 | 13.2 | 11.7 | 0 | 22.7 | 5.3 | 0.958 | 22.7 | 2.73 | 11.7 | 1 | 22.1 |

### B·D²/N (B = signed points, D = online GAE per target = decomp + la + verify, N = r)

- `koblitz-orbit[divisor=1] + mitm-frobenius + — + walk`: B=2.02e+03 D=1.12e+03 → B·D²/N = 1.21e+03
- `koblitz-orbit[divisor=1] + mitm[negation_folded=1] + — + walk`: B=2.02e+03 D=1.11e+03 → B·D²/N = 1.2e+03
- `binary-subspace[dimension=11] + mitm[negation_folded=1] + — + walk`: B=2.08e+03 D=3.09e+03 → B·D²/N = 9.5e+03
