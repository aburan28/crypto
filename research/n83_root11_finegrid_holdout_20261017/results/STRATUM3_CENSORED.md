# N83 root11 finegrid: independently replayed censored stratum 3

The corrected independent verifier replayed all six sealed public points,
relations, and scalars from the preserved stratum 3 provider ZIP. It recorded
`PASS_INDEPENDENT_REPLAY`, `PASS` resources and event order, and empty local
stderr. The host is `NOT_ADMITTED` because two before/after calibration
pairs at point 18 exceeded the frozen 10% limit. The local verifier exited 2
after writing the canonical analysis. The full 24-point Tier 1 selector
therefore remains `BOUNDED_UNKNOWN` even though K3072 has the shortest
observed target interval in this stratum.

| Point 18 policy | Before calibration ms | After calibration ms | Relative spread |
| --- | ---: | ---: | ---: |
| K2816/root11 | 463.828635 | 525.773321 | 12.519% |
| K3072/root11 | 522.707942 | 464.747311 | 11.739% |

The relative spread is the absolute difference divided by the mean of the
two calibration intervals. All reported point measurements below passed
the independent relation and scalar replay, but the host gate excludes
them from a promoted full-panel comparison.

| Point index | K2543 online ms | K2816 online ms | K3072 online ms | Same-point rho for K3072 ms |
| ---: | ---: | ---: | ---: | ---: |
| 18 | 23.665896 | 10.388522 | 3.659816 | 1,714.959900 |
| 19 | 14.586919 | 39.611636 | 10.374900 | 2,256.691944 |
| 20 | 54.096871 | 15.655679 | 5.600986 | 2,449.279679 |
| 21 | 60.882844 | 63.124597 | 1.120223 | 1,213.422930 |
| 22 | 67.378728 | 21.500692 | 15.529024 | 1,383.403946 |
| 23 | 172.637994 | 36.956281 | 38.547192 | 524.461565 |

The descriptive geometric target-online means are 48.629761 ms for K2543,
26.195778 ms for K2816, and 7.228111 ms for K3072. The [analysis JSON](stratum3_analysis.json)
retains per-point cold costs, exclusive online phases, base hashes, rank,
resources, calibration, host topology, and scalar replay. The
[provider receipt](stratum3_provider_receipt.json) binds artifact
`12108510285`; its 202,741,349-byte ZIP has SHA-256
`5d6dc325968e356b489d6add655afe59a4e8cff8f5fc24178e6d0a4b3e5fec1d`.
The [local replay receipt](stratum3_local_replay_receipt.json) binds
corrected verifier SHA-256
`5b9ef05cabce00267239e83a9a7e1104b24305bedc970e69d60c2bf2e4e6b36e`
and canonical analysis SHA-256
`46d495c3349af105eb5babd5f1a42a90c1bdf5a322a5cc0acaa2e6472139948f`.
