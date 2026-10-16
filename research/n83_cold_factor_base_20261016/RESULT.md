# N83 cold factor-base screen: K1024 leads the completed root22 arms

On the frozen public point, K1024/root22 completed cold construction, rank acquisition, target handling, and producer verification in **162.363 s**, **163.404 s**, and **150.173 s** on three cloud allocations. It used 169,984 distinct usable base points. K1536/root22 followed at 179.871 s, 180.796 s, and 171.348 s. K256/root22 reached its 900 s cap in the first two allocations and completed in 842.986 s in the third. Independent replay accepted each completed factor base, relation set and answer, plus each allocation's same-point rho answer.

The first allocation's K2543/root22 control differed by 11.53% of its mean, above the frozen 10% noise gate, so that screen retains `BOUNDED_UNKNOWN`. The second allocation's control differed by 9.86% and passes. The corrected verifier ranks completed resource-admitted arms while preserving K256's timeout, yielding `SCREEN_ADMITTED_WITH_TIMEOUT` and advancing K1024/root22 and K1536/root22 to the frozen reserve-point comparison. The third allocation's A/A spread was 20.87%, so its result remains `BOUNDED_UNKNOWN` despite all arms completing and passing independent replay. The selected factor base remains null until the reserve-point comparison. The first two allocations exposed four AMD EPYC 7763 CPUs; the third exposed four AMD EPYC 9V74 CPUs. These allocations provide descriptive cold-stage evidence, not a Tier 1 online speedup admission.

| Arm | Usable base points | Folded columns | Cold process, s | One-target online, ms | Peak RSS, GiB | Replay status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| K2543/root22 A | 422,138 | 2,543 | 346.955 | 25.173 | 10.80 | pass |
| K256/root22 | unverified (42,496 planned) | 256 planned | unknown (900 s outer cap) | unknown | 0.10 sampled | timeout |
| K512/root22 | 84,992 | 512 | 378.198 | 321.549 | 0.44 | pass |
| K768/root22 | 127,488 | 768 | 230.107 | 1,416.595 | 0.90 | pass |
| K1024/root22 | 169,984 | 1,024 | **162.363** | 72.512 | 1.63 | pass |
| K1536/root22 | 254,976 | 1,536 | 179.871 | 451.149 | 3.43 | pass |
| K2048/root22 | 339,968 | 2,048 | 230.955 | 109.908 | 6.30 | pass |
| K3450/root11/wide3 | 572,700 | 3,450 | 605.336 | 74.066 | 11.51 | pass |
| K2543/root22 B | 422,138 | 2,543 | 309.133 | 21.597 | 10.79 | pass |

The paired rho process recovered and independently replayed the same point in 721.293 ms online. The IC online intervals above begin after reusable base, root-index, and factor-log preparation; the cold intervals include those stages. The rho pair was run once after the IC arms, so this screen does not establish a controlled IC/rho speedup ratio. The K256 timeout preserved its exit status, 900.422 s outer wall, 3,334.782 user CPU seconds, and resource receipt.

The [IC1 row record](ic1_rows_38863274637.json) binds the eight completed rows to seven immutable candidate manifests and workload `eebf98ee1031`, including the frozen generator domain and SHA-256 seed. The K1024/root22 candidate is `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0hcc9ec7592620`. The two K2543 controls share one candidate ID and have distinct `R1`/`R2` run IDs. K256's base dump was not produced before its cap, so its candidate and run IDs remain null while its timeout receipt remains in the raw artifact.

`materialize.py` checks the archived artifact hash, source and workload identities, independent replay digest, each complete base dump, and each canonical candidate hash before writing those manifests. It is idempotent on the saved artifact:

```sh
python3 research/n83_cold_factor_base_20261016/materialize.py \
  /path/to/extracted-artifact/results \
  --archive /path/to/n83-cold-factor-base-38863274637.tar.zst
```

The [original analysis](analysis_38863274637.json), [executed-source receipt](source_38863274637.json), and [host receipt](host_38863274637.json) are checked into this directory. A [second local replay](analysis_replay_38863274637.json) against the same raw artifact completed with `evidence_status: PASS_INDEPENDENT_REPLAY` and `status: BOUNDED_UNKNOWN`; this verifier revision makes evidence integrity distinct from factor-base selection. The full raw artifact is [GitHub Actions run 38863274637](https://github.com/aburan28/crypto/actions/runs/38863274637), artifact ID `12032697449`, including all stdout, stderr, base dumps, usage and verifier records. A local copy is archived at `/Volumes/SSD990-2/n83-archive-20261016/raw/n83-cold-factor-base-38863274637.tar.zst`; its SHA-256 is `57b4e0391ca86dc2af40d0b4a00ac642d41069e563b91c34a041022b5506b6fa` (50.0 MiB compressed). The original analysis SHA-256 is `d8c24607bf77cad5c446e70e9ec202e0226a9c031fc73a71232d52f1247d633a`; the local replay analysis SHA-256 is `609b86d76e6f50054711a3ab452f55f319e8b17706a3d12eed76fafc6cf60e06`.

## Second allocation and timeout decision

[GitHub Actions run 38873482289](https://github.com/aburan28/crypto/actions/runs/38873482289) completed all nine bounded attempts and the rho control. Its original cloud verifier reported `BOUNDED_UNKNOWN` because it required even the capped K256 arm to complete before ranking any arm. The frozen protocol instead calls for ranking **completed** resource-admitted pipelines while retaining timeout rows. The corrected verifier was replayed independently against the full downloaded artifact. Its [replay analysis](analysis_replay_38873482289.json) records `PASS_INDEPENDENT_REPLAY`, `SCREEN_ADMITTED_WITH_TIMEOUT`, `bounded_arms: [k256root22]`, and confirmation candidates `k1024root22` and `k1536root22`. No raw producer receipt or original cloud analysis was changed.

| Arm | Second cold process, s | One-target online, ms | Replay status |
| --- | ---: | ---: | --- |
| K2543/root22 A | 346.785 | 25.825 | pass |
| K256/root22 | 900 s cap | unknown | timeout |
| K512/root22 | 382.783 | 327.998 | pass |
| K768/root22 | 226.281 | 1,420.160 | pass |
| K1024/root22 | **163.404** | 72.304 | pass |
| K1536/root22 | 180.796 | 452.993 | pass |
| K2048/root22 | 233.318 | 111.596 | pass |
| K3450/root11/wide3 | 620.368 | 75.173 | pass |
| K2543/root22 B | 314.188 | 21.238 | pass |

The second same-point rho online interval was 722.126 ms. The original analysis SHA-256 is `9e2c36063e5cce484e82d5ff65f41cc09947c9efb0d033329bb4128d3a1375f9`; the corrected replay analysis SHA-256 is `c49d7d2e11f77a6be8f514f06eb45a3046157845da778d4e5be09469cee79709`. Artifact ID `12038165977` is archived byte-for-byte at `/Volumes/SSD990-2/n83-archive-20261016/raw/n83-cold-factor-base-38873482289.zip`, SHA-256 `367079f7a5ce915b2bef01eaa2665af872efc1d5bd0323cb921e6bd3e3d37a46`.

## Current-head third allocation

[GitHub Actions run 38882773376](https://github.com/aburan28/crypto/actions/runs/38882773376) completed successfully at source commit `644fabdfca2423609c6166a445ed95bf2acca132`. Its verifier independently replayed all nine IC arms and same-point rho. Unlike the earlier allocations, K256/root22 completed its factor base, rank and target in 842.986 s, with 42,496 usable points and 256 folded columns. K1024/root22 remained the fastest cold arm at 150.173 s; K1536/root22 took 171.348 s and K3450/root11/wide3 took 594.241 s. The same-point rho online interval was 634.032 ms. All completed arms passed the resource checks, with zero recorded cgroup throttle or OOM events.

The K2543/root22 A/B cold intervals were 368.583 s and 298.934 s, a 20.87% relative spread. The frozen A/A gate therefore records `BOUNDED_UNKNOWN` and no confirmation candidates from this allocation. This run exposed four AMD EPYC 9V74 CPUs. Its [original analysis](analysis_38882773376.json), [executed-source receipt](source_38882773376.json), and [host receipt](host_38882773376.json) are checked in. Raw artifact ID `12040969445` is archived byte-for-byte at `/Volumes/SSD990-2/n83-archive-20261016/raw/n83-cold-factor-base-38882773376.zip`, SHA-256 `db9429e3d880f808e392d1954e9cf059a464f4103ac0ab8b40614dbf59be6e43`. A fresh local invocation of the independent verifier on the extracted artifact produced a byte-identical analysis with SHA-256 `eb04275ea56bb083e54a5e6d931189bbda32d8ad090ee608374b9e0b284a3b32`.

Next, run the frozen four reserve points with K1024/root22, K1536/root22, K1152/root22, and K1152/dynamic root11 in rotated order, retaining A/A controls and same-point rho receipts. The K1152 finalists come from the subsequent seven-arm neighborhood screen. The search covers these specified contiguous-orbit bases and index variants; a wider factor-base search remains open.
