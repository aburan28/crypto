# P-256 million-curve grid: verified result

Status: **complete; negative preregistered screen; durable publication blocked**

Date: 2026-10-06

## Outcome

The frozen 1,000 by 1,000 product grid completed and independently replayed
exactly **1,000,000 distinct canonical curves** in the registered P-256
isogeny class.  It contains 999,999 explicit parent isogenies and 999,999
kernel certificates.  Generation and replay agree on every count, histogram,
extremum, short-identity collision, terminal summary and rolling-chain value.

No preregistered structural screen fired.  This is a negative model-level
screening result, not an ECDLP experiment:

- `qr_prefix_64` ranged from 14 through 51; there were no values at the
  exploratory comparison threshold 52 or the primary threshold 54;
- the smallest canonical signed `b` length was 233 bits, above 224;
- among models without `a = -3`, the smallest signed `a` length was 236 bits,
  above 224; and
- no repeated j-invariant, canonical model, full ICV1 identity, 48-bit EC1
  alias or full curve UID occurred.

No factor-base relation test or whole index-calculus solve was selected by the
frozen rule.  **ECDLP speedup remains unset.**

## Requested scope and status

| requirement | status | evidence or remaining gap |
|:--|:--|:--|
| generate 1,000,000 P-256-isogenous curves locally | verified complete | exactly 1,000,000 unique canonical identities |
| independently replay construction and certificates | verified complete | all 999,999 parent edges and two fresh order-audit points per curve passed |
| search frozen structural anomalies and discrepancies | verified complete for this grid | no threshold hit or full-identity discrepancy; expected 32-bit display collisions recorded |
| test an index-calculus speedup | not selected by the frozen protocol | no screen crossed its threshold; relation yield and a full solver were not run |
| support local cache plus on-demand S3 synchronization | implemented and locally verified | durable upload blocked by absent ambient AWS identity and AWS CLI |
| coordinate the durable receipt through Cairn | implemented but not submitted | no Cairn node, objective or signing identity is configured |
| search `2^32`, `2^40`, or every P-256 isogeny | partial | this run covers 1,000,000 curves, or 0.023283064% of `2^32` |

## Population and correctness gates

| gate | generated | independent replay | result |
|:--|--:|--:|:--|
| unique curves | 1,000,000 | 1,000,000 | pass |
| parent edges | 999,999 | 999,999 | pass |
| kernel certificates | 999,999 | 999,999 | pass |
| prime-order audits | 1,000,000 at one point | 1,000,000 at two new points | pass |
| non-singular / valid generator | 1,000,000 / 1,000,000 | same | pass |
| unique j / model / full ICV1 / EC1 / UID | 1,000,000 each | same | pass |
| uncompressed bytes | 1,773,118,481 | 1,773,118,481 | pass |
| record-chain SHA-256 | `4de6e93e0d7a90ff53a26f8f1fd1cbd03ab55f701a25f8e49874cb8bead40302` | same | pass |

The verifier regenerated every modular-polynomial root choice, predecessor
exclusion, kernel polynomial, Velu codomain, target isomorphism, canonical
model, generator, identity and detector.  It used two audit points beginning
at `x = 7`, independent of generation's one audit point beginning at `x = 0`.

The typed transfer and its still-open solver obligations are recorded in
`TRANSFER_ASSESSMENT.json`.  In brief:

```text
E_(0,0)/F_p -- degree-13 spine --> E_(0,y)/F_p
             -- degree-11 row --> E_(x,y)/F_p
```

Every edge is a separable cyclic isogeny over `F_p`.  Its induced map on the
prime-order rational group is an isomorphism because 11 and 13 are coprime to
that order.  This establishes subgroup-preserving transfer, not a reduction in
discrete-log cost.

The shorter display slug is deliberately not a uniqueness gate.  There were
121 distinct two-way 32-bit suffix collisions and no triple collision, close
to the birthday expectation
`1,000,000 * 999,999 / (2 * 2^32) = 116.415`.  All longer identities remained
unique.

## Structural screens

| statistic | frozen threshold | observed | decision |
|:--|:--|:--|:--|
| `qr_prefix_64` primary | at least 54 | maximum 51; zero hits | no candidate |
| `qr_prefix_64` comparison only | at least 52 | zero hits | no candidate |
| canonical `b` signed bits | at most 224 | minimum 233 at index 4,289, coordinate `(293,3)` | no candidate |
| non-`a=-3` signed `a` bits | at most 224 | minimum 236 at index 143,043, coordinate `(186,142)` | no candidate |

The sole curve at `qr_prefix_64 = 51` is index 867,068 at coordinate
`(935,866)`, full UID
`urn:ec-record:1:sha256:0cab9d6f7b3526de6f2862a44dbadf621940adb5cbd28434a0f1edc015a089af`.
It did not cross a frozen selection threshold.

Exploratory distribution moments were mean `31.996739`, population variance
`16.2977284`, skewness `0.0011113`, and excess kurtosis `-0.0372237`.  The
simple independent `Binomial(64, 1/2)` reference has mean 32, variance 16,
skewness 0 and excess kurtosis -0.03125.  The earlier 65,536-curve multigraph
had mean `31.9990082` and variance `16.1231527`.  The million-grid variance is
1.86% above the simple reference, but the deterministic isogeny-grid samples
and the 64 probes within a curve are not independent.  This post-execution
moment comparison is not a registered anomaly and does not select a curve.

## Boundary accounting

| variant | unique curves | certified edges represented | fraction of `2^32` | ratio to `2^32` | ECDLP speedup |
|:--|--:|--:|--:|--:|:--|
| verified multi-prime reference | 65,536 | 302,653 | 0.001525879% | 1 / 65,536 | unset |
| verified degree-13/11 product grid | 1,000,000 | 999,999 | 0.023283064% | 1 / 4,294.967296 | unset |

The new population is 15.258789 times the reference curve count.  It is not an
exhaustive enumeration of the P-256 isogeny class and is only
`1,000,000 / 2^32` of the previously requested screening boundary.

## Execution receipt

The clean generation source was
`316f546373db8f42ac894ef7d8b4232e57d87f8e`.  The exact release binary was
3,522,416 bytes with SHA-256
`37ccab8a4076cf04d06ffd87d00eb9f1b76fbbc411f6e5bc688f675656231b71`.
The GitHub-published equivalent production tree is
`c6daafdf0dc8913bd72f0f0d8ba8ecaa8fc46ee9` at commit
`86cf65727adf50451b48008149098b470e038d77`; the different commit ID is solely
the consequence of publishing through GitHub's object API from a different
parent commit with an identical base tree.

| phase | threads | wall | user CPU | system CPU | peak RSS |
|:--|--:|--:|--:|--:|--:|
| generation, one audit point | 4 | 5,816.209 s | 22,687.881 s | 3.265 s | 1,010,831,360 B |
| replay, two audit points | 5 | 6,601.237 s | 23,164.372 s | 3.878 s | 1,011,425,280 B |

Host: Linux 6.18.44 x86-64 KVM, five exposed cores of an AMD EPYC 9V74,
18,882,699,264 bytes RAM; Rust `1.98.0 (88d9e12ae 2026-08-18)`.

```text
p256_isogeny_million --threads 4 generate --side 1000 --batch-rows 8 \
  --audit-points 1 \
  --source-commit 316f546373db8f42ac894ef7d8b4232e57d87f8e \
  --output p256-grid-1m.jsonl.gz

p256_isogeny_million --threads 5 verify \
  --input p256-grid-1m.jsonl.gz --audit-points 2 \
  --audit-seed-x 7 --batch-rows 8
```

## Artifact and local cache

| object | bytes | SHA-256 |
|:--|--:|:--|
| canonical gzip | 527,157,501 | `87522b6de08e64804dcee17bba058eb91922d6b8dd9853cd32523f67a3d0cbdc` |
| exact decompression | 1,773,118,481 | `2a1a27d06f48b22b85ce3c96ac4017f5af6a0256e9d1ad277e3fb1f516475454` |
| generation receipt | 24,346 | `2e330d01a1df8ba4e722e1ac5d8dc0d25e80ddc4d1c811c70e8db10d806cbba3` |
| replay receipt | 24,388 | `0cf69c6180e2f9759f88982e3cd05787f442772c28dd2ffadf7ab1a9962529df` |
| cache receipt | 570 | `55408727b68f9a39489f3b3e0ef728b5e03f2be5ae3d8c80c6c7f1d27a478f9e` |
| cache manifest | 92,736 | `4e4cceae95d7cc11632204877f29e2a608d2d4b60a211f177fd6868c884241c3` |

Post-run storage code is commit
`53aff5c2a4ccc1c924a83e0844a62c1cf082a6a7`; its release binary is 4,424,584
bytes with SHA-256
`1e334c63b8aac63948627049f43101aece2156cc7fb0a3b93afae9b4dbd7f6b7`.
Its GitHub-published equivalent tree is
`25fd5ae0304fad8127af2ccd420f5f7b07b9619b` at commit
`59a2915b9eed1012a4a34270bf6e3c62c03d1011`.
It derived one preamble, 125 eight-row members and one summary member.  An
external concatenation of all 127 decompressions reproduced the canonical
content hash and byte count above.  Row members contain 7,992 curve records
and average 4,213,504 stored bytes.  All cache members total 527,260,083 stored
bytes.

The canonical artifact and cache remain under
`/tmp/p256-isogeny-million-20261006/` for this task.  That path is local
evidence, not durable publication.

## S3 and Cairn status

The implemented hybrid path uploads the canonical artifact, both receipts and
all row members under `objects/sha256/<first-two>/<digest>`, with a declared
SHA-256 and `If-None-Match: *`.  It creates
`runs/p256-grid-1m-20261006/complete.json` last.  Fetches validate both stored
and decompressed hashes before an atomic rename.  A compact Cairn
`coordination-receipt` names the S3 marker and canonical hashes through the
existing signed, persisted commit-reveal transport; it is explicitly labelled
`s3-content-address-only`.

No upload or Cairn submission was attempted.  This managed environment has no
configured AWS outbound identity, no AWS CLI, and no Cairn node, objective or
identity.  Chat-provided access keys were not used.  Durable publication is
therefore **blocked on ambient worker identity and destination configuration**,
not on artifact preparation.

## Validation and retained failures

- 3 by 3 generation/replay, deterministic-gzip and all registered tamper
  cases passed before production.
- The relevant isogeny-walk tests, release build, `cargo check`, and clippy
  with warnings denied passed.
- Post-run storage validation passed incremental SHA known-answer/boundary
  tests, exact slice reconstruction, malformed-partition and byte-tamper
  rejection, S3 location checks, compile and clippy.
- One broader pre-run `cargo test --test curve_id` attempt exhausted the
  disposable debug link workspace and the linker exited with a bus error
  before tests ran.  It was an infrastructure failure, not recorded as a
  passing test; the failed attempt was reported and the relevant focused
  tests subsequently passed.

## Conclusion and next achievable goal

Within this fixed million-curve degree-13/11 grid, there is no structural
discrepancy and no preregistered candidate for an index-calculus experiment.
The result does not support the claim that index calculus is faster on an
isogenous P-256 model.

The next scientifically achievable step is a separately frozen control study
of the slight QR-count overdispersion, using disjoint walks and translated or
pseudorandom probe sets.  Only if that study pre-registers and reproduces a
selection effect should it advance candidates to matched factor-base yield
tests and then cold, full-cost one-target IC-versus-rho solves.
