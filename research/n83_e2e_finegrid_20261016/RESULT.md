# N83 adjacent-K cold screen: K1107 leads the exposed-point panel

The complete K1102–K1120 screen replayed 21 IC processes and one same-point
rho process on GitHub Actions run [38978215742](https://github.com/aburan28/crypto/actions/runs/38978215742).
All 21 bases, four-point relations, and recovered scalars passed independent
replay. The source, arm order, topology, CPU telemetry, resource envelope, and
the two A/B noise controls passed. This is a one-point exploratory cold-process
screen; its selected candidate remains null pending the frozen fresh-point
holdout.

| K | Actual usable points | Folded columns | Cold process (s) | One-target online (ms) |
| ---: | ---: | ---: | ---: | ---: |
| 1107 | 183,762 | 1,107 | 116.222 | 354.511 |
| 1112 | 184,592 | 1,112 | 117.262 | 59.898 |
| 1111 | 184,426 | 1,111 | 117.598 | 488.007 |
| 1103 | 183,098 | 1,103 | 117.637 | 43.962 |
| 1104 | 183,264 | 1,104 | 119.440 | 51.613 |
| 1120 | 185,920 | 1,120 | 119.632 | 23.300 |

All values are from one previously exposed public point; A/B values are means
of two fresh processes, and other values are single fresh processes. The
same-point rho online interval was 720.959 ms. The online numbers are stage
diagnostics from this cold experiment; the primary online comparison needs
its own matched fresh-point Tier 1 panel.

The raw provider ZIP SHA-256 is
`b07cef5779312fed6fb49c2af767d5d2c8cbdc1eca6c14de85ddd6753fb6290f`.
The local replay exited 0, wrote no stderr, and reproduced the archived
analysis exactly as parsed JSON. `ic1_rows_38978215742.json` binds the raw artifact to
19 exact candidate manifests, 21 run IDs, the workload ID
`eebf98ee1031`, stage and resource measurements, and correctness receipts.

The next screen compares root11 and root22 at K1107 and K1112. The sealed
sixteen-point fixture in `research/n83_e2e_final_holdout_20261016` is reserved
for candidate selection after the exposed-point screens are complete.
