# N83 K1107/K1112 root-layout screen: noise gate inconclusive

[Cloud run 38993637162](https://github.com/aburan28/crypto/actions/runs/38993637162)
completed all eight IC arms and the same-point rho arm on source commit
`3ab9dc25ff6ae4139669676342e72391458fe608`. The cloud verifier reported
`PASS_INDEPENDENT_REPLAY` for every point, base, relation and scalar, plus
passing execution order, topology, host telemetry and resource gates. Two
K1107 A/B pairs exceeded the frozen 5% cold-process spread gate. The panel
therefore reports `BOUNDED_UNKNOWN`, `selection_admitted: false`, and no
selected candidate. The verifier's exit code 2 caused the workflow failure;
it represents the preregistered noise outcome.

| Pipeline | Actual usable points | A cold process (s) | B cold process (s) | A/B spread | Gate |
| --- | ---: | ---: | ---: | ---: | --- |
| K1107/root11 | 183,762 | 115.911 | 110.247 | 5.0087% | noisy |
| K1107/root22 | 183,762 | 114.822 | 108.801 | 5.3850% | noisy |
| K1112/root11 | 184,592 | 115.106 | 110.699 | 3.9041% | pass |
| K1112/root22 | 184,592 | 111.706 | 110.832 | 0.7855% | pass |

The exposed point and source were identical across the layout arms. At fixed
K, the root11 and root22 arms had matching independently replayed factor-base
hashes. The observed root22 minus root11 cold differences were -1.089 s and
-1.446 s at K1107, and -3.401 s and +0.133 s at K1112. These are diagnostic
paired values; the failed noise gate precludes ranking the four pipelines.
The rho online interval was 666.923 ms and independently replayed on the same
point.

The provider retained artifact `12087241838` under the name
`n83-e2e-root22-neighbors-38993637162-1` (35,059,707 bytes). Download,
provider-digest verification, independent local replay, and immutable IC1
row materialization remain pending. The next measured action is a complete
rerun of the frozen eight-arm schedule, retaining this failed pair as a row.

The repository also retains **configuration-only** candidate manifests for
K1107/root22
(`IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0hcb73613ebc5a`)
and K1112/root22
(`IC1N83Ckb1a1fb184592PDP4rootRCguidedLAgaussTDdirectISO0hbc5ec3072dc4`).
They use the independently replayed root11 base dumps, which the cloud
verifier found equal to the root22 bases at each K. As a construction check,
the same conversion reproduced the previously archived K1104/root22 and
K1120/root22 manifests exactly. The new IDs are expected identities for the
frozen algorithm; raw artifact binding and measured costs are still pending.
