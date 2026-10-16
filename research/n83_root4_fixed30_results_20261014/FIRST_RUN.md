# N83 fixed30 first execution: exact insertion fault and three replayed controls

The first six-arm execution completed three five-byte control arms and
isolated a direct-state four-byte insertion fault. At source commit
`7294b395961b8353c8be70e7873c5151c51c45be`, every new
`root4f30` arm reached index construction and exited 101 at
`examples/koblitz_orbit_dlp_fast.rs:1690`: the builder called the generic
`insert_if_absent` path on an `Indexed83Four` table, whose explicit guard
requires indexed insertion. The three controls independently replayed all
1,000 public targets each. This receipt remains separate from the corrected
rerun on the same frozen points.

| Arm | Producer/replay status | Mean target online ms | Peak RSS GiB | Recorded outcome |
| --- | --- | ---: | ---: | --- |
| K2940 root4f30 | `FAILED` | — | 7.087 | exit 101 during index insertion |
| K2940 root22 | `PASS_INDEPENDENT_REPLAY` | 48.615448 | 12.719 | 1,000 target scalars and relations replayed |
| K3200 root4f30 | `FAILED` | — | 8.386 | exit 101 during index insertion |
| K3200 root21 | `PASS_INDEPENDENT_REPLAY` | 38.605925 | 13.714 | 1,000 target scalars and relations replayed |
| K3300 root4f30 | `FAILED` | — | 8.916 | exit 101 during index insertion |
| K3300 root11 | `PASS_INDEPENDENT_REPLAY` | 33.085301 | 10.053 | 1,000 target scalars and relations replayed |

The [Actions run](https://github.com/aburan28/crypto/actions/runs/38596699717)
used source SHA-256
`68f9b22972ee09511e15ae74a998ba74b82e74eed36e3014bd0460c918b6cdcb`.
Its seven-step job passed source/fixture freeze, release example tests, all
bounded child executions, and raw artifact upload; the independent verifier
step exited 1 because the three new arms had producer failures. Artifact
`11917005743` has SHA-256
`bd270ee4609c503c99b07c43a186dd46a661219af53fd4c0adc45fddf7319c2f`.
The retained `results/analysis.json` has SHA-256
`eb541897edce4449411131d608f70b226c01a189802aa6cf0121508fdc0ebbd8`;
each failed arm keeps its stderr, exit status, wall time, peak RSS and
stop reason. The runner reported an AMD EPYC 7763 and 16,373,452 KiB
host memory.

The corrected producer selects `insert_indexed83` from the allocated table
variant for both exact four- and five-byte slots. Its source SHA-256 is
`6b22a065d95d5481342c1b1cc37dc9479200dd4f0832bcd53ae6ce4b8c9312d5`.
The second execution retains the same six arms, public points, limits and
resource guard and receives a distinct Actions run and artifact identity.
