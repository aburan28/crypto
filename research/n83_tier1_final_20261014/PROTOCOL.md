# N83 fresh-point Tier 1 final comparison

Five new public points are frozen here before the secondary factor-base
selection result is observed. `generate_fixtures.py` derives each scalar as
`1 + little_endian(SHA256("n83-tier1-final-20261014-v1" || u32_le(i))) mod (r-1)`
for `i=0..4`, then multiplies the exact N83 subgroup generator.
`targets.jsonl` has SHA-256
`7d80c73fa7930ae0c22113bad3f253d8dfb5c9c71b1ce8fbeedd06f45649556b`;
the independent answer sidecar has SHA-256
`03e299325a2e65e353dc5128d7785617d70a4841092a4b9a45903448a5b60b42`.
The generator verifies subgroup membership and excludes all five archived
selection points and the earlier capacity point. The IC producer receives
one public point per process and does not read the scalar sidecar.

The exact candidate will be frozen in `selection.json` after the secondary
panel completes. The selected manifest must bind its canonical `IC1` ID,
factor-base size `K`, state layout, source SHA-256
`07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae`,
14-GiB observed peak guard, child cap at most 1500 seconds, and the sealed
selection-analysis digest. The selected arm is the lowest arithmetic mean
of five verified online intervals among panel arms inside the 14-GiB
envelope. A failed or censored panel arm cannot be selected.

## Paired execution

Run all five pairs on one dedicated cloud allocation. For point indices
0, 2 and 4, run IC then rho; for 1 and 3, run rho then IC. Each IC arm is a
fresh process that constructs its exact factor base, index and factor logs
before its single point arrives at the target stage. Record setup separately.
The strong rho arm uses the archived `fold-L64-DP6-E0` configuration, one
thread and a 120-second cap. Both arms receive the same public point,
source-pinned binaries, declared four-CPU affinity, cgroup resource envelope
and independent scalar replay. All failed attempts and raw timing receipts
remain in the artifact.

The runner pins itself and children to four visible CPUs, preferring
distinct physical cores as exposed by the instance. It records allocation
metadata, visible topology, cgroup limits, raw `/proc/stat` and interrupt
snapshots, CPU and memory events, `wait4` peak, timestamps and hashes.
Three fixed SHA-256 calibration passes run immediately before and after
each pair. The verifier requires no cgroup throttling or OOM events, an
observed peak at or below 14 GiB, mean steal time per execution CPU at most
1% of the pair wall time, and maximum/minimum calibration median at most
1.15 across the panel. These thresholds were frozen before the run.

The independent checker replays each target's subgroup membership and
scalar, each IC process's complete base and four-point relation, and rho's
same-point answer. It checks the five exclusive IC online phases against
the producer interval and requires the same base hash across processes.
The Tier 1 speedup gate requires all five pairs to pass correctness and
noise gates, with IC online time below same-point rho on all five points.
Report every ratio, the median ratio, a deterministic 10,000-resample
target-block percentile interval for that median (seed `2026101405`),
arithmetic mean online times, setup and memory peaks, and the cloud instance
class. The paired result describes this frozen target law and instance tier.

## Selected candidate bound after panel replay

The independently replayed secondary panel selected K3240 root21 under the
rule above. `selection.json` has SHA-256
`fbe9757f05b661ae4a2b1fefe161a8ccf35cc51c9c22934003f7bd0619fd7094`
and binds canonical candidate
`IC1N83Ckb1a1fb537840PDP4rootRCguidedLAgaussTDdirectISO0hff882a640d8e`.
The panel selection record SHA-256 is
`d8f5b10961a77127f3f43fbe35ae3004a8c5ea14202c8e950ddc6b63078aa240`;
its analysis SHA-256 is
`e23c04abc39a2c517860c325ec66ab767f86e11b9ac41aa988a875d1d26b9663`.
The selected candidate manifest SHA-256 is
`8fc0c9ecdda418171afee38616720ac8d4d5ce9784b5dd9290d6d42444a276d3`.
The runner checks all three files and the exact enumerated-base hash before
starting; the independent checker compares every fresh process's base hash
to that selected manifest. The selected child cap is 1500 seconds per point.
