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
