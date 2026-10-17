# N83 current-source online leader comparison, frozen before the finalist result

The exact K3072/root11 candidate is
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a`.
Its local arithmetic and candidate identity replay is recorded in
`RESULT.md`. The same IC source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
must execute both comparison candidates. The second candidate is the frozen
N83 finalist protocol's selected **online** candidate if its two Tier 1
host strata are admitted and `online.status == SELECTED`. Otherwise use the
fixed fallback K1120/root11,
`IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234`.
A fallback measurement preserves the finalist panel's null selection.
Commit and independently replay the finalist selection and its provider
artifact binding before dispatching this comparison.

The ten sealed public points are in `targets.jsonl`, SHA-256
`9bf9bfca72d4e82a5780759d52ac000b812ca25e831af04238983c8b72d2d777`.
The answer sidecar, SHA-256
`3f7a75dc0134f883a77fa2b87eb6839b6df804c53950d8dc77e7f29b2500284c`,
is used only for independent replay. `generate.py --verify-sealed` must
pass before measurement. Each point is a previously unused, one-target
workload. Run points 0–4 and 5–9 on separate dedicated Tier 1 cloud
allocations. Within each allocation, run both candidates on every point in
fresh processes with empty base, index and log caches. Each candidate-point
arm has its own strong rho solve of the **same** public point. Run the two
candidates in opposite orders on successive local point indices. Alternate
`IC,rho` and `rho,IC` for each candidate across the five points. Run a fixed
exposed-point rho calibration before and after each IC–rho pair, using the
same seed `2026101601` and public point as the earlier Tier 1 online panel.
The paired rho seed is `2026101900 + 2*point_index + candidate_index`, where
candidate index 0 is K3072 and 1 is the selected finalist or fallback.

Archive the exact binaries, source and input hashes, commands, timestamps,
cgroup limits, affinity, topology, CPU frequency, steal time, throttling,
OOM events and peak memory. Cap K3072 IC processes at 900 seconds, finalist
IC processes at 480 seconds, and every rho or calibration process at 120
seconds. Cap the producer at 130 minutes and the entire cloud job at 180
minutes. Preserve every timeout, failed replay, noisy calibration and
global-cap stop. A host stratum is admitted only if all ten IC–rho pairs,
their before/after calibrations, source and candidate IDs, relation/base
replay, scalar replay, resource envelope and event order pass. A calibration
pair passes when its two online times differ by at most 10% of their mean.
The two strata may be pooled only after each independently passes Tier 1.
Long-horizon repeats may be retained as diagnostics without deciding
admission, provided the required pairs and calibrations are complete.

For each candidate, report all ten `rho_online_ms / IC_online_ms` values,
their median and host-specific means, and the paired 95% interval for rho
minus IC online time. A Tier 1 speedup for a host requires IC online time
below rho on all five of that host's public points. To choose between the two
IC candidates, calculate ten matched-point log time ratios
`log(finalist_online_ms / K3072_online_ms)`, their mean and descriptive 95%
paired interval, and the exact two-sided sign-flip p-value over all `2^10`
sign assignments. Select K3072 only if its mean log ratio is positive and
the exact p-value is at most 0.05; select the finalist only if the mean is
negative and the same gate passes. Otherwise record a null selected
candidate. Report complete cold-process costs and paired differences as a
supplementary objective, with no substitution for online selection.
