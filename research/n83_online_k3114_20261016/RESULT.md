# K3072 wins the N83 K3114 online frontier screen

The exact-root N83 one-point screen selected K3072/root11 over the proved
K3114 capacity frontier under its predeclared A/A and strict-separation
rules. K3072's two online intervals were 8.373915 and 7.606348 ms;
K3114's were 17.036886 and 17.051287 ms. Both bases were constructed
twice from empty processes, had identical base and index identities within
each A/B pair, passed resource admission, and independently replayed the
target relation and scalar. The same-point rho reference independently
replayed at 627.662905 ms online.

| Exact-root factor base | Actual usable points | Folded columns | Online A/B (ms) | A/A spread | Target probes | Peak RSS A/B (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K3072/root11 | 509,952 | 3,072 | 8.373915 / 7.606348 | 9.606% | 19,875 | 9.403 / 9.404 |
| K3114/root11 | 516,924 | 3,114 | 17.036886 / 17.051287 | 0.0845% | 45,668 | 9.519 / 9.518 |

K3072's slower interval was 8.662971 ms faster than K3114's faster
interval, exceeding the largest observed absolute A/A difference of
0.767567 ms. The target-dependent point-decomposition phase dominated:
K3072 required 19,875 exact-index probes on this point and K3114 required
45,668. The larger base therefore did not improve this frozen point's
probe trajectory. This is an exposed-point factor-base comparison on the
GitHub-hosted four-vCPU envelope; the five fresh-point Tier 1 panel is a
separate acceptance experiment.

The [run](https://github.com/aburan28/crypto/actions/runs/38916776487)
used source commit `9abf6afe8bc52ad4f136a19ff1603d60c344f30d` and
artifact `12054666847`. Its exact local ZIP is
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-online-k3114-38916776487.zip`,
SHA-256 `08f09b113c7256d829a606f2e4c630ccd77e09605060002e654a615ccd32b4bc`;
ZIP integrity passed. The [analysis](analysis_38916776487.json) SHA-256
is `7a5f69618b000dac4fdf2a94626ce0639620868cf0960c00bcab6dbbda9113de`.
Independent local replay on the executed source recreated the cloud
analysis byte for byte; replay log SHA-256 is
`c4130401e059e9b5a7285c60d832c1221bf4f77363626922ee756604529490fe`.
All four IC arms and rho passed independent replay; no cgroup throttle or
OOM events occurred. The host exposed CPUs `0,1,2,3` and 16,373,452 KiB
of memory. The [IC1 rows](ic1_rows_38916776487.json) bind all four
measured arms to two exact candidate manifests. K3072's selected ID is
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hafe77bc68a40`;
K3114's ID is
`IC1N83Ckb1a1fb516924PDP4rootRCguidedLAgaussTDdirectISO0h1103f7d44ce0`.
The materializer checked the raw ZIP and replay hashes, frozen source,
base dumps, actual usable points and A/A identities, then ran twice
without changing any versioned record.
