# N83 five-point factor-base selection: K3240, 537,840 usable points

The frozen selection rule chose **K3240 with the 21-byte retained root state**
from seven complete four-summand factor-base pipelines. Its 3,240 folded
columns reached full rank, and all five archived public targets passed
independent base, relation, subgroup, and scalar replay. The arithmetic mean
of the five target-online intervals was **14.4220794 ms** on the cloud
instance. The observed `wait4` process peak was 14,955,081,728 bytes
(13.928 GiB), 77,303,808 bytes below the 14-GiB guard.

| Arm | Usable points | Folded columns | Five-point mean online ms | Setup s | `wait4` peak GiB | Status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| K3400 root11 | 564,400 | 3,400 | 31.9082072 | 639.3 | 10.351 | replay pass |
| **K3240 root21** | **537,840** | **3,240** | **14.4220794** | **616.3** | **13.928** | **selected; replay pass** |
| K3200 root11 | 531,200 | 3,200 | 51.0535448 | 518.2 | 9.762 | replay pass |
| K3114 root22 | 516,924 | 3,114 | 31.3006784 | 523.1 | 13.638 | replay pass |
| K2940 root22 | 488,040 | 2,940 | 37.8488706 | 449.6 | 12.719 | replay pass |
| K2543 root22 | 422,138 | 2,543 | 71.9605518 | 333.1 | 10.797 | replay pass |
| K2543 root11 | 422,138 | 2,543 | 71.5846808 | 323.2 | 8.048 | replay pass |

Each row uses one fresh process with target-independent setup followed by
the same five archived target points. Setup is recorded separately from the
online intervals. The K2543 root11/root22 controls have identical base hash,
exact root-key count, and target probe counts on all five points. The
selected candidate is
`IC1N83Ckb1a1fb537840PDP4rootRCguidedLAgaussTDdirectISO0hff882a640d8e`;
[`selected_candidate.json`](selected_candidate.json) binds its canonical
manifest, source-base hash, workload, and artifact digests.

The selected mean has a five-target bootstrap interval of
7.7648412–22.242011 ms. The next fastest measured arm was K3114 root22;
its paired mean time minus the selected mean was 16.878599 ms, with a
target-block bootstrap interval of −0.2540498–34.0112478 ms. This overlap
makes the pre-frozen fresh-point Tier 1 same-point IC/rho experiment the
next measurement. The selected K3240 arm has only 73.7 MiB of observed
headroom, so that run must retain its own resource-gate result.

## Reproducibility

[Cloud run 38528070197](https://github.com/aburan28/crypto/actions/runs/38528070197)
executed commit `6f070ffd3b65d11549afd8253aeb4d62efc193ce` and producer source
SHA-256 `07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae`.
Its five-point input SHA-256 is
`b96672b4d75df4b90a337c1e471d54dc147edbf236ef8d0054adb716d2314871`.
Artifact ID `11879829587` has raw ZIP SHA-256
`ef8b0221257789d8e45672c6125d102e6f11374947869ee03a8f4bd88c2fc469`;
`results/analysis.json` has SHA-256
`e23c04abc39a2c517860c325ec66ab767f86e11b9ac41aa988a875d1d26b9663`.
The independent cloud verifier and a second local run of
`python3 research/n83_factor_base_panel_20261014/verify.py EXTRACTED_RESULTS`
both exited zero and produced that same analysis digest.

`materialize.py` checks the raw ZIP, analysis, executed source, workload,
and complete factor-base dumps before regenerating the seven immutable
`IC1` manifests, workload, measurement rows, and selection record. The
GitHub artifact retains every per-arm raw stdout, target result, base dump,
resource receipt, build/test log, and host receipt for 14 days.
