# N83 five-arm current-source challenger on ten sealed public points

K1280/root11 dynamic reached **4.73325 ms** for one exposed N83 public point,
with 212,480 usable factor-base points, 1,280 folded columns and a 1.34 GB
packed root index. Its current-source arithmetic replay and the K1248
control are [versioned](../n83_smallk_current_source_20261017/RESULT.md).
K2543 and K3072 are [current-source controls](../n83_k2543_current_source_20261017/RESULT.md)
on the same point. The prospective question is which factor-base size gives
the shortest verified **one-target online** interval over unused points under
matched same-point rho and Tier 1 cloud resources. Cold process cost is a
separate supplementary metric.

The four fixed candidate IDs are:

| Arm | Actual usable base points | Folded columns | Exact candidate ID |
| --- | ---: | ---: | --- |
| K1248/root11 dynamic | 207,168 | 1,248 | `IC1N83Ckb1a1fb207168PDP4rootRCguidedLAgaussTDdirectISO0hc9eda193920f` |
| K1280/root11 dynamic | 212,480 | 1,280 | `IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dac` |
| K2543/root11 dynamic | 422,138 | 2,543 | `IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55` |
| K3072/root11 dynamic | 509,952 | 3,072 | `IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a` |

The fifth arm is the frozen eight-arm finalist protocol's selected **online**
candidate if both Tier 1 strata are admitted and `online.status == SELECTED`.
Otherwise it is the predeclared K1120/root11 fallback,
`IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234`.
A fallback measurement leaves the finalist selection null. The same IC source
SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
must execute every arm. The curve is `EC1N83Ckb1a1haacb69ee9f37`, GF(2^83)
with polynomial integer `9671406556917033397649557`, equation
`y² + xy = x³ + x² + 1`, and prime subgroup order `8569786107849059`.

The ten public points were sealed in the K3072 followup fixture **before**
the K2543, K1248 and K1280 exposed-point results. This panel copies those
target and answer bytes exactly. `targets.jsonl` SHA-256 is
`9bf9bfca72d4e82a5780759d52ac000b812ca25e831af04238983c8b72d2d777`;
the answer sidecar SHA-256 is
`3f7a75dc0134f883a77fa2b87eb6839b6df804c53950d8dc77e7f29b2500284c`.
The answer sidecar is used only for independent replay. The inherited
`generate.py --verify-sealed` must reproduce the fixture and its original
manifest before measurement. This successor panel supersedes the dormant
three-arm K2543 challenge; neither fixture has been dispatched. The finalist
result and provider archive binding must be committed and independently
replayed before a later `prepare_dispatch.py` creates `dispatch.json`. Only a
push of that file arms the cloud workflow.

Run public points 0–4 and 5–9 on separate Tier 1 cloud allocations. Each
host runs all five candidates on its five points, each candidate-point arm
in a fresh process with empty base, index and log caches. Each arm gets an
independent strong rho solve of **the same point**. Rotate candidate order
cyclically by global point index, placing every candidate in each position
once per host. Alternate `IC,rho` and `rho,IC` by point and fixed arm index.
Before and after each pair, run the fixed exposed-point rho calibration using
the earlier Tier 1 calibration point and seed. The paired rho seed is
`2026102500 + 5*point_index + candidate_index`, where fixed indices 0–4 are
K1248, K1280, K2543, K3072 and finalist/fallback. This makes 25 IC–rho
pairs and 50 calibration processes per host, 50 pairs over the ten points.

Preserve exact source and binary hashes, commands, timestamps, affinity,
topology, cgroup limits, frequency, steal time, throttling, OOM events and
peak memory. Cap K1248, K1280 and finalist IC processes at 480 seconds,
K2543/K3072 at 900 seconds, and each rho or calibration process at 120
seconds. The producer cap is 190 minutes per host; the job cap is 260
minutes. Preserve every timeout, failed replay, noisy calibration and
global-cap stop. A host stratum is admitted only if all 25 IC–rho pairs,
their before/after calibrations, base and relation replay, scalar replay,
source and candidate IDs, resource envelope, event order, and Tier 1
telemetry pass. A calibration pair passes when its two online times differ by
at most 10% of their mean. The two strata may be pooled only after each
independently passes Tier 1; physical host isolation is not required for
this tier.

For each candidate, report all ten paired `rho_online_ms / IC_online_ms`
values, their median and host-specific means, and the paired 95% interval
for rho minus IC online time. A host-specific Tier 1 speedup requires IC
online time below rho on all five of that host's public points. Compare all
ten unordered candidate pairs on their ten matched points: log ratios for
online intervals, ordinary differences for supplementary cold process cost,
with means, descriptive 95% paired intervals, and exact two-sided sign-flip
p-values over all `2^10` sign assignments. Select a candidate only if it is
faster in the observed direction than each of the other four and all four
corresponding p-values are at most `0.05/10 = 0.005`. Otherwise record a null
selection. Never replace the primary online choice with the cold choice.
