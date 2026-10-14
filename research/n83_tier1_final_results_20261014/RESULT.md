# N83 five-point Tier 1 paired result

The selected K3240 root21 index recovered each of five previously unseen
public-point scalars faster online than the same-point strong-rho reference on
one GitHub Actions allocation. The median paired rho/IC ratio was **54.10**;
the frozen 10,000-resample target-block bootstrap interval was
**[17.18, 58.90]**. All five four-point relations and both recovered scalars
per point passed independent replay. The calibration maximum/minimum ratio was
1.0182 against the preregistered 1.15 bound; each pair passed the recorded
CPU-throttle, OOM, steal-time, memory and calibration gates.

## Exact candidate, input and allocation

- Candidate:
  `IC1N83Ckb1a1fb537840PDP4rootRCguidedLAgaussTDdirectISO0hff882a640d8e`.
  Its complete manifest is the hash-pinned
  `research/n83_factor_base_panel_results_20261014/candidate_K3240_root21.json`.
  The base contains 537,840 usable subgroup points before sign/Frobenius
  folding and 3,240 effective columns; its complete base hash is
  `54847e372e3d68a2509bb16227aef70e09ceb9948d9e73149b6abc21355a2c0e`.
- Curve: `EC1N83Ckb1a1haacb69ee9f37`, the binary Koblitz model and prime
  subgroup specified in the candidate manifest. The five public targets were
  generated from the frozen SHA-256 scalar law in
  `research/n83_tier1_final_20261014/generate_fixtures.py`; public target
  SHA-256 is `7d80c73fa7930ae0c22113bad3f253d8dfb5c9c71b1ce8fbeedd06f45649556b`.
  They are disjoint from the factor-base selection targets. Each row in
  `workloads.json` fixes one public point and has its own IC1 workload ID.
- Allocation: GitHub Actions run
  [38538483710](https://github.com/aburan28/crypto/actions/runs/38538483710),
  `ubuntu-24.04`, x86-64, AMD EPYC 7763, image `20261004.327.1`, instance
  `gh-38538483710-1-five-fresh-pairs`. The same pinned exposed CPU set
  `{0,1,2,3}` and 14 GiB RSS guard applied to all alternating pairs. The
  cgroup reported `cpu.max = max 100000` and `memory.max = max`; exact
  affinity, timestamps, counters, limits, binary hashes and raw failures are
  in `host.json`, `source.json` and each pair's `usage.json` in the artifact.

## One-target online results

Each row is one fresh public-point workload. The IC clock begins with its
target-dependent query after base, index and factor logs are ready and ends
after producer recovery checking. The rho clock covers its same-point walk
through producer scalar verification. Independent replay is recorded outside
both producer intervals. IC setup is reported separately and was rebuilt in
the fresh process for each point. `measurements.json` retains all five
exclusive IC online phases, target points, hashes, memory, noise and process
receipts.

| Point | IC online (ms) | Rho online (ms) | Rho / IC | IC setup (s) | IC peak RSS (GiB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 55.927 | 2,051.409 | 36.68 | 639.554 | 13.928 |
| 2 | 35.793 | 1,936.234 | 54.10 | 605.382 | 13.928 |
| 3 | 85.190 | 5,018.004 | 58.90 | 589.510 | 13.929 |
| 4 | 50.084 | 2,819.416 | 56.29 | 584.095 | 13.928 |
| 5 | 67.550 | 1,160.618 | 17.18 | 580.409 | 13.928 |

The arithmetic means were 58.909 ms IC and 2,597.136 ms rho. The observed
maximum IC peak was 14,956,396,544 bytes, below the 15,032,385,536-byte
guard. Setup times above are target-independent preparation and are not part
of the displayed online ratio. The row-level phase sums, full process wall
times and resource counters remain in the sealed measurement file.

## Evidence and replay

The producer source commit was
`04e6a33bb67f37177b7d251982676e9dc7959682`; the IC source SHA-256 was
`07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae`.
GitHub artifact `11884791327` had downloaded ZIP SHA-256
`fb2fad6fa3994a9aed37d8854fafff20ef803bcc522d220d8b78fe65e6e1a4bc`,
matching GitHub's artifact digest. The cloud `analysis.json` SHA-256 was
`66395ad10b13cfa237f1ba17ac0180aa45ce4790107093215d373d4ca3644cb6`.
An independent local invocation of
`python3 research/n83_tier1_final_20261014/verify.py <extracted-results>`
exited 0 and reproduced that analysis digest exactly. The versioned
`materialize.py` checks those digests plus the selection, manifest, source,
target and answer hashes before writing five canonical run rows; a second
invocation reproduced the same rows without changes. The raw artifact is
retained by GitHub for 14 days and the downloaded ZIP plus replay log are in
`/private/tmp/n83-final-38538483710` on this host.

An automatic branch rerun after main acquired the four-byte producer
([run 38582925768](https://github.com/aburan28/crypto/actions/runs/38582925768))
stopped at the frozen-source preflight. Its failed receipt is preserved. The
versioned workflow now uses manual dispatch and checks out the measured
`04e6a33...` source commit explicitly for a repeat allocation.

The later 100-point secondary panel selected K2940/root22 as its descriptive
mean minimum under the same 14 GiB guard. Its separate five-point reserve
fixture is the next validation of that changed candidate; this K3240 result
remains the completed fresh-point reference for that decision.
