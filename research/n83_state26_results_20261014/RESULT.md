# N83 26-byte states: 484,720 usable points within 14 GiB

Packing the two 83-bit roots and three state labels into 26 bytes completed
K2920 with **484,720 subgroup-usable factor-base points**, 2,920 folded
columns, 707,677,409 exact root keys, and a 13.936-GiB `wait4` peak. The
independent checker replayed every base point, subgroup membership,
four-point target relation, scalar, and exclusive online phase. All five
arms passed correctness replay.

| Arm | Usable points | Root keys | Target probes | Online ms | Setup ms | `wait4` peak GiB | 14-GiB envelope |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| K2940 state26 | 488,040 | 717,404,838 | 16,117 | 10.487786 | 520,121.517 | 14.054 | Exceeded observed peak |
| K2920 state26 | 484,720 | 707,677,409 | 240,882 | 132.934704 | 490,641.715 | 13.936 | Within |
| K2900 state26 | 481,400 | 698,016,386 | 116,155 | 63.162314 | 475,719.719 | 13.814 | Within |
| K2543 state26 | 422,138 | 536,737,856 | 16,174 | 8.590523 | 355,158.063 | 11.794 | Within |
| K2543 state28 | 422,138 | 536,737,856 | 16,174 | 8.525601 | 355,773.431 | 12.296 | Within |

The K2940 producer and independent replay completed, but its final
15,090,896,896-byte peak exceeded the 15,032,385,536-byte guard by
58,511,360 bytes. Its row retains verified mathematics and the observed
resource-envelope failure separately. K2920 is the largest completed arm
inside that envelope in this sweep. At K2543, the 26-byte and 28-byte
controls have identical base hash, exact root-key count, and scalar.

The same-point strong-rho interval was 1,874.384228 ms. The online numbers
above are individual cloud-instance observations on a previously used
public point. `measurements.json` retains exclusive phases, setup, process
usage, correctness, and canonical candidate IDs. Base selection will use
the archived secondary target panel; a new point and five alternating
Tier 1 pairs will then test the selected candidate.

## Evidence and replay

[Cloud run 38514041440](https://github.com/aburan28/crypto/actions/runs/38514041440)
executed commit `d29a9ca1cbb5ca8e0da09b57a8804664877160d7`, source SHA-256
`0d78b826504daa6f9d05579286241ce4bfa011ea3978b9c6a91d0cfad0000058`,
and public point SHA-256
`21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8`.
GitHub Actions artifact ID `11872642267` has SHA-256
`473b5ae6b35666413afd33fac74f6ec1e12bcb35cfd216b9d0a1de7c0cabf0f4`.
Its `results/analysis.json` SHA-256 is
`0dad359865c0240e83eaafa2bc8cbd0da36242c50567c6b4068f64c0e36d032b`.
Download that artifact, extract it, then run
`python3 materialize.py --artifact-zip PATH_TO_DOWNLOADED_ZIP PATH_TO_EXTRACTED_RESULTS`
to check the raw archive and regenerate the canonical `IC1` candidate and
one-point workload records. The complete cloud verifier passed before the
artifact was uploaded.

The subsequent 22-byte root-only sweep completed K3114 with 516,924
usable points under the same 14-GiB guard. The new 21-byte fixed-slot gate
tests a larger base with fallback arms.
