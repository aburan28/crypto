# N83 mirrored one-point online finalist screen

All eight frozen IC arms and the same-point rho arm passed independent base,
relation and scalar replay. The two exact finalist manifests and all eight
[IC1 run rows](ic1_rows_38905299528.json) are versioned. K3072/root11
measured a 6.594 ms mean one-target online interval over four fresh
processes; K2543/root11 measured 11.894 ms. The same-point rho interval was
483.986 ms. These are the recorded intervals for the previously exposed
screening point.

| Exact candidate | Four online intervals, ms | Adjacent A/A spreads | Admission |
| --- | --- | --- | --- |
| `IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0h074bc8345b97` (K3072/root11) | 7.364, 6.285, 6.214, 6.512 | 15.80%, 4.69% | Noisy first pair |
| `IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h67d9bea98295` (K2543/root11) | 12.087, 11.822, 11.652, 12.016 | 2.22%, 3.07% | Both pairs pass |

Each candidate's four arms had one identical base hash and exact-index entry
count: K3072 `c5511247735bc183359d5425c4dd794aaad49ad38d2ca964799e1caa76afaa08`
with 783,270,873 entries, and K2543
`eb0a315471ae0d01598fc24d10b8dae41873cc3fab556fc91f85a584880ccf58`
with 536,737,856 entries. All arm resource receipts passed. The K3072
first adjacent pair exceeded the predeclared 10% A/A gate, so
[`analysis_38905299528.json`](analysis_38905299528.json) records
`BOUNDED_UNKNOWN`, `selected_pipeline: null` and
`selected_candidate_id: null`. The measured interval difference does not
override that gate.

The frozen [cloud run](https://github.com/aburan28/crypto/actions/runs/38905299528)
used commit `6670a9aeee73ddf2e942927d6256970bdd4a1de6` and producer source
SHA-256 `3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`.
Artifact `12050698461` is archived at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-online-rescreen-38905299528.zip`,
SHA-256 `e33bcce6c7b794db1060f2a88afb8615a7d12d5f7ba4b7ef408d545bd5792212`.
The ZIP passed its integrity check. Local replay with the byte-identical
frozen verifier reproduced the cloud analysis byte for byte; its SHA-256 is
`205ebad3997139c2448d80d6a2328b60e609a7c8e9016391f14ca43c3625589f`.
The local replay log SHA-256 is
`132c23abcf8a4d01fe7a78866a4419b3a6d4fa973a7a3dd123d6d12b034800a5`.
The IC1 materializer checked the raw archive and replay hashes, base dumps,
source, arm order and exact candidate identities, then ran twice without
changing a versioned record.

The next selection experiment runs K3072 A/B and K2543 A/B on the one
approved `g5g.metal` Tier 1 allocation, seals an admitted candidate as a
versioned S3 object before the five fresh-point producer, and preserves this
screen as a separate noisy cloud-host result.
