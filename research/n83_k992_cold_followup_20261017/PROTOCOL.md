# N83 K992 cold-frontier follow-up on ten sealed public points

The fresh-point question is whether the K992/root11 factor base lowers the
complete one-target cold process time against K1107/root11, the sixteen-point
Tier 1 cold-frontier size, while retaining a verified target-online solve.
K1024/root11 is a prespecified intermediate-size control. The primary IC
comparison remains one-target **online** wall time against an independently
solved, same-point rho target on the same host; cold process cost, including
base, root index and factor-base log preparation, is a separately reported
factor-base objective. Neither objective borrows the other's timing interval.

The [sealed fixture](targets.jsonl) contains ten public points on
`EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)` with modulus
integer `9671406556917033397649557`, curve `y²+xy=x³+x²+1`, prime
subgroup order `8569786107849059`, and cofactor `1128547018`.
`generate.py --verify-sealed` independently regenerates every point and
checks subgroup membership and non-overlap with all earlier precision,
finalist, and five-arm points. The scalar fixture law is
`1 + (little_endian(SHA256(domain || u32_le(index))) mod (r - 1))`, using domain
`n83-k992-cold-followup-ten-20261017-v1`. The public-point file SHA-256 is
`89c7208dd4fc7a5621cf727b3fac9a6714f04f66a312e516b418c31764e7547e`;
the verifier-only [answer sidecar](answers.jsonl) SHA-256 is
`363adc7f9bded3f4a25a80c7b27da5a2305afb30ee6741932cbe76f455c2865b`.
These bytes were sealed after the exposed-point K992 screen and before any
fresh-point K992 measurement or the outcome of the queued five-arm panel.

| Fixed arm | Usable points | Folded columns | Exact current-source candidate |
| --- | ---: | ---: | --- |
| K992/root11 | 164,672 | 992 | `IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h829ed476f094` |
| K1024/root11 | 169,984 | 1,024 | `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h3eb055926289` |
| K1107/root11 | 183,762 | 1,107 | `IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h82329aa062fa` |

The three configurations use the same source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and current-source commit `416489e09c3d8885b6c12f262935e6ba29373460`.
All use the symmetric S3 state quotient, exact five-byte indexed root slots,
root11 reconstruction, dynamic 75%-load request, sparse canonical map,
four guided rank threads and one target thread. Their enumerated base digests
and exact implementation flags belong in immutable IC1 manifests before
dispatch. The [local screen](local_screen/RESULT.md) is a one-point L0
feasibility diagnostic, not an input to the ten-point
selection beyond the three fixed candidate choices.

Run points 0–4 and 5–9 on independent Tier 1 cloud instances. Each instance
runs each of the three arms once per point in a fresh process with empty
base, index and log caches. Give every IC process its own same-point rho
solve. Rotate arm order by point index so each arm occupies early, middle
and late positions; alternate IC then rho and rho then IC by fixed arm
index and point index. Before and after each IC–rho pair, run the fixed
exposed-point rho calibration. This is fifteen IC–rho pairs and thirty
before/after calibrations per host, thirty pairs and sixty calibrations
overall. Record source and binary hashes, exact commands, all exclusive
online phases, cold preparation, operation and rank counts, peak RSS,
topology, cgroup limits, affinity, throttling, steal and event order.

The IC child cap is 480 seconds, each rho/calibration child cap 120 seconds,
and the per-host producer cap 120 minutes. Preserve timeout and failure
records. A host stratum is admitted only when all fifteen IC solves, rho
solves and their arithmetic replays pass, all calibrations pass a 10% paired
spread gate, the frozen order and identities match, and exposed affinity,
cgroup and telemetry receipts satisfy Tier 1. Pool only two independently
admitted strata. Physical-host isolation is an additional Tier 0 property;
it is not required here.

For each arm, report ten paired rho/IC online ratios and the host-stratified
median and uncertainty, plus complete cold process time. Test the three
unordered candidate pairs on the same ten points with exact two-sided
`2^10` sign flips: mean log IC time ratios for online, mean cold millisecond
differences for cold. Use Bonferroni `0.05/3` per pair. Select an exact
online or cold configuration only if it is faster than both others and
passes both corresponding tests. Otherwise preserve a null selection and
the bounded data. The same-point online result is the primary rho
comparison. A fresh-point rank/full-solve failure or noisy stratum remains
in the evidence and makes its selection unknown.

The fixture and method are frozen independently of the pending five-arm
outcome. A later dispatch must bind completed candidate manifests, binary
hashes, runner/verifier source, resource envelope and a fresh cloud run;
no cloud allocation starts from this fixture commit alone.
