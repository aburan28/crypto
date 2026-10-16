# N83 100-point factor-base selection

The preregistered 100-point panel selected **K2940 root22**, with 488,040
usable subgroup factor-base points and 2,940 folded columns. Its mean
single-target online interval was **22.548 ms** across 100 distinct archived
public points under the 14 GiB guard. All twelve candidate arms independently
replayed every recovered scalar, subgroup/base point and four-point relation
in the cloud artifact. The selected candidate is
`IC1N83Ckb1a1fb488040PDP4rootRCguidedLAgaussTDdirectISO0h751130834947`.

The closest measured arm was K3220 root21 at 22.882 ms, a paired mean
difference of 0.334 ms. Its preregistered 10,000-resample target-block
bootstrap interval for runner-up minus selected was **[-5.681, 6.004] ms**.
The ordering is descriptive on this sample. The separate five-point reserve
fixture was frozen before this panel and is the next same-point IC/rho
validation for K2940.

## Frozen method and allocation

The source commit was `93d19669bd6118f5448c37ed8001fdc8183f41ea`,
with IC source SHA-256
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
The 100 public target file SHA-256 was
`fa0044f7cad85185812ac4c79f4feabb300ddcdc1732b6a0890778d5b8e72af6`;
the independent scalar sidecar SHA-256 was
`fcb167978cda4a72732c322442a3b77d3f0a997b74676a2a4d5767eb6c5262bf`.
Each arm used the same N83 curve, rank seed 29, symmetric S3 traversal,
four rank threads, one target thread, and a fresh process with its own
base/index/log setup. The table below reports the arithmetic mean of the
100 target-dependent online intervals after setup, plus the full process's
observed peak RSS. The setup, rank, final LA, five online phases, each
point's probes and time, and failures are in `measurements.json` and the raw
artifact. This secondary panel selects a candidate; each target interval is
kept as an individual observation.

| Arm | Usable points | Mean online (ms) | Peak RSS (GiB) | Status |
| --- | ---: | ---: | ---: | --- |
| K3240 root21 | 537,840 | 26.174 | 13.928 | PASS |
| K3400 root4 | 564,400 | 31.161 | 13.351 | PASS |
| K3300 root11 | 547,800 | 24.335 | 10.052 | PASS |
| K3500 root4 | 581,000 | 28.563 | 13.656 | PASS |
| K3220 root21 | 534,520 | 22.882 | 13.821 | PASS |
| K3550 root4 | 589,300 | 25.063 | 13.813 | PASS |
| K3114 root22 | 516,924 | 28.832 | 13.637 | PASS |
| K3596 root4 | 596,936 | 23.171 | 13.959 | PASS |
| K3200 root11 | 531,200 | 25.685 | 9.762 | PASS |
| K3400 root5 | 564,400 | 30.215 | 10.351 | PASS |
| K3190 root21 | 529,540 | 24.817 | 13.660 | PASS |
| **K2940 root22** | **488,040** | **22.548** | **12.718** | **PASS; selected** |

The K3400 four-byte/five-byte controls had identical complete base hashes,
root-key counts and target probe counts. Candidate-specific manifests bind
each exact enumerated point-set digest and all algorithm flags; the compact
candidate names distinguish bases and layouts beyond `N83` and point count.

## Evidence and next gate

[GitHub Actions run 38543929816](https://github.com/aburan28/crypto/actions/runs/38543929816)
completed on one Ubuntu 24.04 allocation. Raw artifact `11888953582` had
downloaded ZIP SHA-256
`caea68c20f6890a45cd6415979ffd078b50d136f90f27b412f9e6e5612d72f9c`,
matching GitHub's artifact digest. The cloud independent analysis SHA-256 was
`e75000c34a3b363a5c2282a2bcd77765846fa7777c34ec78522f63ca32124f83`.
The versioned `materialize.py` checks those digests, the executed commit,
target sidecar and every complete base dump before generating twelve
canonical IC1 manifests, the 100-point workload ID
`f0caa6eaed63`, measurement rows, and `selected_candidate.json`.
The selected manifest SHA-256 is
`b6f0d5bf6301506ada0483d1aacdeeea4bc5c95bd1b2a272d59831f7e290b9c8`;
the selected base hash is
`12c164dd4bbc72df19d8b0c638f1d808c534bd464397c9033822dace74cca92b`.
An independent local execution of `verify.py` against the downloaded raw
artifact exited 0 after replaying all twelve complete bases and all 1,200
four-point target relations and scalars. It reproduced the cloud analysis
SHA-256 exactly:
`e75000c34a3b363a5c2282a2bcd77765846fa7777c34ec78522f63ca32124f83`.
The local run was a correctness replay; its elapsed time is outside the cloud
target-online intervals.

The reserve five public points have target SHA-256
`644c56a7bae2c92085704bb6f58ab70b38b278710de90152ece7631058c41dc9`
and answer-sidecar SHA-256
`b30206982765ba649dcdf6f1eb48ded2cbe0fd3c8ec6a226a1fdfb261e5873b4`.
They are disjoint from all 100 selection targets and earlier five-point
validation targets. The next run uses those points for five alternating
same-point K2940/root22 and strong-rho pairs under the Tier 1 protocol.
