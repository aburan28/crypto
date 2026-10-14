# N83 exact five-byte index on a 16 GB cloud VM

The exact five-byte root index completed three point-defined N83 factor-base
candidates on one public target in a fresh four-vCPU Linux VM. The largest
tested base, K1798 with 298,468 usable points and 1,798 folded columns, gave
the shortest measured target interval: 60.890432 ms. Its process peaked at
13,105,192,960 bytes RSS and its five-byte root table occupied 5 GiB.
The same-point strong-rho reference completed in 1,879.377195 ms, giving an
exploratory ratio of 30.8649 for this single execution.

| Candidate | Usable points | Folded columns | Verified target probes | Online ms | Setup ms | Peak RSS bytes | Rho / IC, exploratory |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `IC1N83Ckb1a1fb298468PDP4rootRCguidedLAgaussTDdirectISO0hc4398d2c7c5a` | 298,468 | 1,798 | 111,004 | 60.890432 | 440,245.479376 | 13,105,192,960 | 30.8649 |
| `IC1N83Ckb1a1fb249000PDP4rootRCguidedLAgaussTDdirectISO0h09ee3f7d63fb` | 249,000 | 1,500 | 223,726 | 110.369533 | 289,055.889317 | 10,785,693,696 | 17.0280 |
| `IC1N83Ckb1a1fb210986PDP4rootRCguidedLAgaussTDdirectISO0hafa8513db0d5` | 210,986 | 1,271 | 682,144 | 337.220614 | 206,943.017525 | 6,588,059,648 | 5.5731 |

The input, curve, source digest, ordered arms, limits, and interval boundary
are frozen in [the protocol](../n83_indexed5_tier1_20261013/PROTOCOL.md).
Each candidate's manifest is in `candidate_K*.json`; `workload.json` fixes the
single point and `measurements.json` retains the full run and phase rows.
The original [Actions run](https://github.com/aburan28/crypto/actions/runs/38486670301)
executed commit `8954fe0571bc5bb7c152c7f12d8ff387ae250ee0`.

The complete raw artifact, including base dumps, per-run stdout/stderr,
usage, host and cgroup receipts, build logs, and the verifier output, is
`raw-artifact-38486670301.zip`. GitHub artifact ID `11861225775` published
SHA-256 `f30d5f879dea2c7f7fd23c1d46dd0b1a2b8335b2985e766fbb1f40972ea2776b`;
the committed ZIP has the same digest. `results/analysis.json` inside it has
SHA-256 `ae0bb93ea9f2470d37c76b424a834b8f8afc790e038a2a9909a5c054854ea98e`.
An independent local replay of the extracted artifact reproduced that exact
analysis digest and returned `PASS_INDEPENDENT_REPLAY` for all three arms and
rho. The replay checks the complete base, subgroup, four-point relation,
scalar, phase sum, and archived first-witness parity.

The VM identified an AMD EPYC 7763, four visible CPUs, Linux x86-64, and
16,373,452 KiB host-reported memory. The K1798 process completed in
511.840448 s; its target-independent setup was 440.245479 s. The cgroup
recorded no OOM kills or CPU throttling. This is an instance-scoped Tier 1
correctness and capacity result with one run per arm. The current timing
ratios remain exploratory until the prescribed alternating paired-run noise
gate is executed. The existing 64-point panel and the independent 256-point
panel answer different empirical questions; this single target establishes
the best completed candidate only among these three frozen bases on this
point and VM.

The next experiment is the symmetric S3 state quotient. S3 is symmetric in
the two input x coordinates; paired state labels related by exchange and
Frobenius rotation may share one canonical-root orbit. A checked quotient
could reduce index states and permit larger base searches within the same
memory envelope. It requires exact-key and first-witness tests before any
timing comparison.
