# N83 exact-root online factor-base range screens

Two frozen same-point cloud screens independently replayed every completed
factor base, four-point target relation and scalar, plus the matched rho
answer. K3072 with one-root dynamic exact slots had the lowest completed
one-target online interval in both screens: 7.384 ms in the eight-arm screen
and 7.282 ms in the expanded twelve-arm screen. K2543 with the same one-root
layout measured 13.896 and 13.882 ms. The expanded screen included the
K1152/K1536 cold-shortlisted pipelines; their online intervals are retained
in the [expanded analysis](analysis_38896107047.json).

| Pipeline | Eight-arm online, ms | Expanded online, ms | Status |
| --- | ---: | ---: | --- |
| K3072/root11 dynamic | 7.384 | 7.282 | Independently replayed in both |
| K2543/root11 dynamic | 13.896 | 13.882 | Independently replayed in both |
| K1024/root22 dynamic | 47.880 | 45.289 | Independently replayed in both |
| K3450/root11 fixed | 78.146 | 79.161 | Independently replayed in both |
| K3584/root11 fixed | — | — | 900-second timeout in both |

The expanded screen's [IC1 rows](ic1_rows_38896107047.json) bind its 11
completed arms to 10 exact candidate manifests and retain the K3584 timeout
as a separate row with a null candidate identity. Both K2543/root22 A/B
controls resolve to the same candidate ID,
`IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0hcec5838fbe33`.
K3072/root11 resolves to
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hd45656791e9f`.
The materializer checked the exact raw ZIP and replay hashes, executed source,
base dumps, usable-point counts and arm order; it ran twice without changing
any versioned record. The selected candidate remains null under the A/A gate.
The original eight-arm screen has its own [IC1 rows](ic1_rows_38895743960.json):
seven completed arms bind to six exact manifests, and the K3584 timeout
retains a null candidate identity. Its A/B controls also share one exact
candidate ID. The original materializer checked the historical runner bytes
from commit `e83a78be0b74d22fb34e315006d47e6bd82a1e4c` and was likewise
run twice without changing the versioned records.

The predeclared K2543/root22 A/B control was noisy in both runs:
18.570 versus 14.487 ms (24.70%) in the eight-arm screen and 18.327 versus
14.120 ms (25.93%) in the expanded screen. Both exceed the 10% gate, so
their selected pipelines and candidate IDs remain `null` with
`BOUNDED_UNKNOWN` status. This preserves the directly observed timings
without promoting either screen to the fresh-point selector. The matched
rho online intervals were 723.096 and 647.238 ms; these exploratory
same-point screen rows do not replace the five-point Tier 1 comparison.

The [eight-arm run](https://github.com/aburan28/crypto/actions/runs/38895743960)
used source commit `e83a78be0b74d22fb34e315006d47e6bd82a1e4c` and
artifact `12046277866`. Its exact ZIP is
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-online-range-eight-38895743960.zip`,
SHA-256 `41d4aa0d623967ab9c62cc3ded653e01e9ec049e78596d691d87b7f40f4cea45`.
The [eight-arm analysis](analysis_38895743960.json) SHA-256 is
`5db680f0423f1e4b4c4138303135e826e543cb90f3b427147977ee9be192927c`.
Local replay with byte-identical `run.py` and `verify.py` from that commit
recreated the cloud analysis byte for byte; replay log SHA-256
`236aa527f3335d52f7083d2aed2539d175b32424edf9cc187922f3390279d5bb`.

The [expanded run](https://github.com/aburan28/crypto/actions/runs/38896107047)
used source commit `c001093fb47304814bdee2a6b3d31ea8389d9213` and
artifact `12047173997`. Its exact ZIP is
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-online-range-expanded-38896107047.zip`,
SHA-256 `53cdc2fa1ccc46a4617cf9ccc80609d0bc19a5e6d90b2b6a5d849d33ed8c0606`.
The expanded analysis SHA-256 is
`0e33e7ec12ad35f19ec15314455155a86f5137a998535d9144aa77d5c5859717`.
Local replay with the exact `c001093` producer and verifier scripts
recreated it byte for byte; replay log SHA-256
`96f5a578e01269015168b2b10b79a41a8a4f921769fd4746b38e02da2285a3e7`.

The [mirrored rescreen](../n83_online_rescreen_20261016/PROTOCOL.md) now
measures K3072/root11 and K2543/root11 four times each in adjacent A/A
pairs and mirrored order. It requires exact base/index identity, full
replay, pair spreads at most 10%, and separation exceeding observed A/A
noise before it can name the candidate for the sealed Tier 1 panel.
The recurring admission gap and its exact exit conditions are tracked in
the [N83 online A/A problem record](../../problems/n83-online-aa-noise.md).

The [mirrored screen result](../n83_online_rescreen_20261016/RESULT.md)
independently replayed all eight IC arms and rho. Its K2543 adjacent pairs
passed at 2.22% and 3.07%, while K3072's first pair had a 15.80% spread.
The predeclared gate therefore retained `BOUNDED_UNKNOWN` and a null
selection. Two exact candidate manifests and all eight run rows are versioned
with that result. The next selection screen is on the approved Tier 1 host
before the five fresh-point panel.

The [current-source expanded rerun](https://github.com/aburan28/crypto/actions/runs/38897449665)
used source commit `f0a5fce29a441566811f9afc260bc30a59f18144` and
artifact `12050929059`. Its exact ZIP is
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-online-range-current-38897449665.zip`,
SHA-256 `c6fbc8567699767a89095aae793ffe2394658cdfaaa33a76e10e8243ca85e594`.
The [analysis](analysis_38897449665.json) SHA-256 is
`49e71b6542b655b4097700fb6e8d98d595da23acac17061c9d1167f60b2f9965`.
Local replay recreated the cloud analysis byte for byte; its log SHA-256 is
`f872f2bfccefaec99097c5e817d0d4254add2fa45096d6b5edc178961277551e`.
All 11 completed IC arms and same-point rho passed independent replay; the
K3584 arm reached its 900-second cap. The K2543/root22 A/B control measured
18.355451 and 14.125281 ms, a 26.05% spread above the frozen 10% gate.
K3072/root11 measured 7.376757 ms and K2543/root11 14.076503 ms as stage
diagnostics, with rho at 654.993539 ms. Selection remains `null` with
`BOUNDED_UNKNOWN` status. The [12 IC1 rows](ic1_rows_38897449665.json)
bind the completed arms to 10 exact source-specific candidate manifests and
retain the timeout as a null-identity row. The materializer checked the raw
ZIP, replay analysis, historical source and runner hashes, exact base dumps,
usable point counts, and arm order; a second run left every record unchanged.

The same-instance Tier 1 selector was frozen in
[its protocol](../n83_tier1_online_20261016/PROTOCOL.md) after these screen
results. Its approved EC2 allocation attempt returned AWS account-level
`Blocked` at `RunInstances`, before any instance was created. The attempted
client token was `n83tier1-20261016-1854bb9d`; the EC2 lookup for that token
returned an empty instance set. The separate Tier 1 panel therefore awaits
an available authorized host.
