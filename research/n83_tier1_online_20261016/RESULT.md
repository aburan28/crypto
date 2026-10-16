# N83 Tier 1 same-instance selector allocation

The four-arm selector, exact candidate materializer, sealed S3 selection
bundle, and five-point IC/rho panel are frozen at source commit
`4753f58df6d3e9e8672806453c1b8cf15de36bb3`. The selector preflight
checked its source and input hashes; Python syntax, shell syntax, and the
local noisy-control and positive-selection rehearsals passed. These are
implementation checks, not cloud measurements.

The approved attempt requested one `g5g.metal` Spot instance in `us-east-1b`
with a $0.50/hour bid ceiling, 165-minute watchdog, encrypted 50 GiB root
volume, one-time request, and client token
`n83tier1-20261016-1854bb9d`. The live preflight at
`2026-10-16T14:24:20.866330+00:00` quoted $0.4118/hour and confirmed the
approved subnet and ARM64 Amazon Linux image. The exact EC2 dry run passed.
`RunInstances` returned AWS account-level `Blocked`: the account was not
recognized as a valid account. The exact client-token instance lookup returned
an empty set. No selector process or five-point producer ran on AWS.

The local raw receipt at
`/Volumes/SSD990-2/n83-archive-20261016/aws/n83-tier1-host-selector-20261016/failed-allocation.json`
retains the quote, source commit, client token, planned launch and provider
error. The temporary security group, empty versioned S3 bucket, IAM role and
instance profile were deleted; subsequent checks returned `InvalidGroup.NotFound`,
S3 `404`, and IAM `NoSuchEntity`. The previously configured RunPod research
pod also refused its SSH connection during this attempt, so it supplied no
replacement Tier 1 measurement. The candidate selection and five-point
online speedup fields remain unknown until a suitable authorized host runs
the frozen protocol and its archive passes independent replay.

The first GitHub-hosted Tier 1 route ran on source commit
`e0a1c7bc5987a126865ba5f0399084c2940fd288` as
[run 38917594127](https://github.com/aburan28/crypto/actions/runs/38917594127).
Its full release library suite passed 4,154 tests with 117 ignored and no
failures; the IC example passed 20 tests and rho example four. The selector
stopped before constructing any base because the hosted runner did not expose
four distinct physical cores on one visible NUMA node. The seal and five
fresh-point stages were skipped. Exact artifact `12054875724` is archived at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-tier1-github-38917594127.zip`,
SHA-256 `8a043d3065af58c94866918b47da43f870026a03e2082ae6bef7a220009127d5`;
ZIP integrity passed. This is a host-capacity failure, with selector choice
and panel result unknown. The separately frozen
[two-core hosted route](GITHUB_PROTOCOL.md) retains the same A/A and replay
gates on a feasible cloud-runner resource envelope.

## Two-core hosted selector and sealed candidate

The separately frozen two-core GitHub-hosted selector selected K3072/root11
for the five fresh-point panel. [Run 38921627515](https://github.com/aburan28/crypto/actions/runs/38921627515)
used source commit `4238a20bacba174b634cd9efdf88d048e61a77ff`,
two distinct physical CPUs `0,2` on visible NUMA node `node0`, and the
same previously exposed public point for both bases and rho. Independent
replay passed all four IC arms and rho; no arm recorded cgroup throttling,
OOM or steal time.

| Base | Actual usable points | Online A/B (ms) | A/A spread | Outer cold A/B (s) | Peak RSS A/B (bytes) |
| --- | ---: | ---: | ---: | ---: | ---: |
| K3072/root11 | 509,952 | 7.070100 / 6.689632 | 5.53% | 384.329 / 374.825 | 10,097,074,176 / 10,094,952,448 |
| K2543/root11 | 422,138 | 13.234984 / 12.659500 | 4.44% | 266.742 / 266.743 | 8,643,420,160 / 8,642,473,984 |

Both adjacent A/A spreads passed the predeclared 10% gate, with identical
base hashes and index-entry counts within each pair. K3072's slower online
interval beat K2543's faster interval by 5.589400 ms, exceeding the largest
observed absolute A/A difference of 0.575484 ms. The same-point rho online
interval was 468.241890 ms. The selector sealed exact candidate
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0ha2d30fb17dd0`;
the five fresh-point panel began only after the same job confirmed the
selector artifact ID, digest, run, source and host.

Exact selector artifact `12058065194` has GitHub SHA-256 digest
`fa3f471254ed965ad3f89904dbf3f7e7dfd8f920893217053412da5211c9a030`.
The downloaded ZIP at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-tier1-selector-38921627515.zip`
has that same SHA-256 and passed ZIP integrity. Its analysis SHA-256 is
`a0ec9b360442081d2c00f209bff910edbc3da5fe419df2f98cf38b34489f3e70`.
Independent local replay on the executed source recreated the analysis byte
for byte; replay log SHA-256 is
`4834db768b3888e60b6fcfe0e7c57aa5cfa77538161d2ce2cd722e31770ac637`.
The selector measures an exposed point. The five fresh-point same-host panel
is still running and will supply the separate Tier 1 online comparison.
