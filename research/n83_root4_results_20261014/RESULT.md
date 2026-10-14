# N83 four-byte root slots: 596,936 usable factor-base points

The exact four-byte root-slot layout completed **K3596 with 596,936
subgroup-usable factor-base points** and 3,596 folded columns. Its
1,073,268,803 exact root keys occupied a fixed `2^31`-slot table. The
observed `wait4` process peak was 14,988,685,312 bytes (13.959 GiB),
43,700,224 bytes below the 14-GiB guard. All seven frozen arms, the
same-point strong-rho arm, and independent cloud and local base, relation,
subgroup, and scalar replay passed.

| Arm | Usable points | Root keys | Target probes | One-point online ms | Setup s | `wait4` peak GiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K3400 root4 | 564,400 | 959,460,734 | 49,703 | 37.263316 | 689.3 | 13.350 |
| K3400 root5 | 564,400 | 959,460,734 | 49,703 | 29.448567 | 570.8 | 10.350 |
| K3500 root4 | 581,000 | 1,016,729,255 | 44,677 | 29.173061 | 625.2 | 13.656 |
| **K3550 root4** | **589,300** | **1,045,986,026** | **7,503** | **5.762825** | **655.4** | **13.813** |
| K3596 root4 | 596,936 | 1,073,268,803 | 20,658 | 14.015392 | 665.7 | 13.959 |
| K2543 root4 | 422,138 | 536,737,856 | 16,174 | 8.146978 | 303.6 | 11.048 |
| K2543 root5 | 422,138 | 536,737,856 | 16,174 | 7.803710 | 300.1 | 8.050 |

The K3400 and K2543 root4/root5 controls each have identical complete
base hash, root-key count, target probe count, and recovered scalar; both
four-point relations replay independently. The four-byte slot retains a
two-bit fingerprint and checks the complete canonical root before accepting
a match. K3550 was
the fastest arm on this previously used public point; its one-point
interval is an exploratory target-search observation. The same-point rho
worker recovered and replayed the scalar in 1,249.096758 ms online.
The pre-frozen 100-point secondary panel will compare target-search cost
across the surviving layouts before a new five-point Tier 1 validation.

The 29-bit state ordinal imposes a sharp limit on this encoding:
K3596 has 536,647,262 retained states, below `2^29`, while K3597 has
536,945,772 states, above it. `RootFour.lean` checks the slot encode/decode
identities and both capacity inequalities. The [proof map](../n83_root4_20261014/PROOF_MAP.md)
identifies the remaining formal bridge from those arithmetic lemmas to
the Rust table transitions and canonical GF(2^83) roots.

## Reproducibility

[Cloud run 38530953283](https://github.com/aburan28/crypto/actions/runs/38530953283)
executed commit `278b1a784f6fb42e2085cd10b8d04d2166e3e481` with producer
SHA-256 `5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`
and public-point SHA-256
`21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
Artifact ID `11882662982` has raw ZIP SHA-256
`8c6a609445b8547957452fb9a2ae645fb34c303d18762f7da1539f86f97f7eba`;
`results/analysis.json` has SHA-256
`5787ddc706373530099889383df541ee1270349f117d9427e6e294890775774f`.
The independent cloud verifier and a second local run of
`python3 research/n83_root4_20261014/verify.py EXTRACTED_RESULTS`
both exited zero with that same analysis digest.

`materialize.py` checks the raw ZIP, analysis, executed source, exact
public point, and complete base dumps before regenerating the seven
canonical `IC1` candidate manifests, workload, and measurement rows.
The raw artifact retains each per-arm resource and correctness receipt
for 14 days.
