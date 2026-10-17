# N83 K2543 current-source challenger on ten prospectively sealed points

The exact K2543/root11 dynamic candidate is
`IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55`.
The exact K3072/root11 dynamic candidate is
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a`.
Their local arithmetic and candidate identity replays are recorded in
[`K2543`](../n83_k2543_current_source_20261017/RESULT.md) and
[`K3072`](../n83_k3072_current_source_20261017/RESULT.md). The same IC source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
must execute all three comparison candidates. The exposed-point K2543
control completed in 401.115 seconds with an 18.983 ms online interval,
making it the next distinct base size to compare with K3072 on fresh
points. The third candidate is the frozen
N83 finalist protocol's selected **online** candidate if its two Tier 1
host strata are admitted and `online.status == SELECTED`. Otherwise use the
fixed fallback K1120/root11,
`IC1N83Ckb1a1fb185920PDP4rootRCguidedLAgaussTDdirectISO0hc922b499a234`.
A fallback measurement preserves the finalist panel's null selection.
Commit and independently replay the finalist selection and its provider
artifact binding before dispatching this comparison.
The branch workflow is armed only by a later committed `dispatch.json` with
the SHA-256 digests of those two versioned finalist files and their bound
IC1 run rows. No dispatch file
is present in the frozen fixture commit.

The ten sealed public points were fixed and committed before the K2543
control completed. They are byte-identical to the unused
[`K3072 followup fixture`](../n83_k3072_current_source_20261017/fixture_manifest.json),
which has not been dispatched. They are in `targets.jsonl`, SHA-256
`9bf9bfca72d4e82a5780759d52ac000b812ca25e831af04238983c8b72d2d777`.
The answer sidecar, SHA-256
`3f7a75dc0134f883a77fa2b87eb6839b6df804c53950d8dc77e7f29b2500284c`,
is used only for independent replay. `generate.py --verify-sealed` must
pass before measurement. Each point is a previously unused, one-target
workload. Run points 0–4 and 5–9 on separate dedicated Tier 1 cloud
allocations. Within each allocation, run all three candidates on every point in
fresh processes with empty base, index and log caches. Each candidate-point
arm has its own strong rho solve of the **same** public point. Rotate the three
candidates in cyclic order by global point index. Alternate
`IC,rho` and `rho,IC` by global point and fixed candidate index. Run a fixed
exposed-point rho calibration before and after each IC–rho pair, using the
same seed `2026101601` and public point as the earlier Tier 1 online panel.
The paired rho seed is `2026102300 + 3*point_index + candidate_index`, where
candidate indices 0, 1 and 2 are K2543, K3072 and the selected finalist or
fallback, respectively.

Archive the exact binaries, source and input hashes, commands, timestamps,
cgroup limits, affinity, topology, CPU frequency, steal time, throttling,
OOM events and peak memory. Cap K2543 and K3072 IC processes at 900 seconds,
finalist IC processes at 480 seconds, and every rho or calibration process
at 120 seconds. Cap the producer at 150 minutes and the entire cloud job at 210
minutes. Preserve every timeout, failed replay, noisy calibration and
global-cap stop. A host stratum is admitted only if all fifteen IC–rho pairs,
their before/after calibrations, source and candidate IDs, relation/base
replay, scalar replay, resource envelope and event order pass. A calibration
pair passes when its two online times differ by at most 10% of their mean.
The two strata may be pooled only after each independently passes Tier 1.
Long-horizon repeats may be retained as diagnostics without deciding
admission, provided the required pairs and calibrations are complete.

For each candidate, report all ten `rho_online_ms / IC_online_ms` values,
their median and host-specific means, and the paired 95% interval for rho
minus IC online time. A Tier 1 speedup for a host requires IC online time
below rho on all five of that host's public points. For the three IC
candidates, calculate ten matched-point log time ratios for each of the three
unordered pairs, their means and descriptive 95% paired intervals, and exact
two-sided sign-flip p-values over all `2^10` sign assignments. Select a
candidate only if it is faster in the observed direction than both
alternatives and both p-values are at most `0.05/3`. Otherwise record a null
selected candidate. Apply the same three-pair gate to supplementary complete
cold-process costs without substituting that selection for the online result.
