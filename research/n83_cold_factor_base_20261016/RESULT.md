# N83 one-point cold factor-base screen: K1024 leads the completed root22 arms

On the frozen public point, the K1024/root22 pipeline completed cold construction, rank acquisition, target handling, and producer verification in **162.363 s**. It used 169,984 distinct usable base points and 1.751 GB peak RSS. K1536/root22 followed at 179.871 s; K256/root22 reached the 900 s cap. Independent replay accepted all eight completed factor bases, relation sets, and answers, as well as the same-point rho answer.

The repeated K2543/root22 control took 346.955 s and 309.133 s, a difference of 11.53% of their mean. This exceeds the frozen 10% noise gate. The screen therefore records `BOUNDED_UNKNOWN`, leaves the selected factor base null, and carries the complete rows into the frozen confirmation campaign. These are measurements on a GitHub Actions cloud instance, with four visible AMD EPYC 7763 CPUs and 15.61 GiB host memory; the instance does not have a Tier 1 speedup admission.

| Arm | Usable base points | Folded columns | Cold process, s | One-target online, ms | Peak RSS, GiB | Replay status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| K2543/root22 A | 422,138 | 2,543 | 346.955 | 25.173 | 10.80 | pass |
| K256/root22 | 42,496 | 256 | unknown (900 s outer cap) | unknown | 0.10 sampled | timeout |
| K512/root22 | 84,992 | 512 | 378.198 | 321.549 | 0.44 | pass |
| K768/root22 | 127,488 | 768 | 230.107 | 1,416.595 | 0.90 | pass |
| K1024/root22 | 169,984 | 1,024 | **162.363** | 72.512 | 1.63 | pass |
| K1536/root22 | 254,976 | 1,536 | 179.871 | 451.149 | 3.43 | pass |
| K2048/root22 | 339,968 | 2,048 | 230.955 | 109.908 | 6.30 | pass |
| K3450/root11/wide3 | 572,700 | 3,450 | 605.336 | 74.066 | 11.51 | pass |
| K2543/root22 B | 422,138 | 2,543 | 309.133 | 21.597 | 10.79 | pass |

The paired rho process recovered and independently replayed the same point in 721.293 ms online. The IC online intervals above begin after reusable base, root-index, and factor-log preparation; the cold intervals include those stages. The rho pair was run once after the IC arms, so this screen does not establish a controlled IC/rho speedup ratio. The K256 timeout preserved its exit status, 900.422 s outer wall, 3,334.782 user CPU seconds, and resource receipt.

The [IC1 row record](ic1_rows_38863274637.json) binds the eight completed rows to seven immutable candidate manifests and workload `ee4634a4e091`. The K1024/root22 candidate is `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0hcc9ec7592620`. The two K2543 controls share one candidate ID and have distinct `R1`/`R2` run IDs. K256's base dump was not produced before its cap, so its candidate and run IDs remain null while its timeout receipt remains in the raw artifact.

`materialize.py` checks the archived artifact hash, source and workload identities, independent replay digest, each complete base dump, and each canonical candidate hash before writing those manifests. It is idempotent on the saved artifact:

```sh
python3 research/n83_cold_factor_base_20261016/materialize.py \
  /path/to/extracted-artifact/results \
  --archive /path/to/n83-cold-factor-base-38863274637.tar.zst
```

The [original analysis](analysis_38863274637.json), [executed-source receipt](source_38863274637.json), and [host receipt](host_38863274637.json) are checked into this directory. A [second local replay](analysis_replay_38863274637.json) against the same raw artifact completed with `evidence_status: PASS_INDEPENDENT_REPLAY` and `status: BOUNDED_UNKNOWN`; this verifier revision makes evidence integrity distinct from factor-base selection. The full raw artifact is [GitHub Actions run 38863274637](https://github.com/aburan28/crypto/actions/runs/38863274637), artifact ID `12032697449`, including all stdout, stderr, base dumps, usage and verifier records. A local copy is archived at `/Volumes/SSD990-2/n83-archive-20261016/raw/n83-cold-factor-base-38863274637.tar.zst`; its SHA-256 is `57b4e0391ca86dc2af40d0b4a00ac642d41069e563b91c34a041022b5506b6fa` (50.0 MiB compressed). The original analysis SHA-256 is `d8c24607bf77cad5c446e70e9ec202e0226a9c031fc73a71232d52f1247d633a`; the local replay analysis SHA-256 is `609b86d76e6f50054711a3ab452f55f319e8b17706a3d12eed76fafc6cf60e06`.

Next, run the frozen four reserve points with repeated K1024 and K1536 arms plus K1271's dynamic root11 index, alternating order and same-point rho pairs. Keep resource failures and A/A noise in the resulting rows. The present screen does not select a global optimum from the unbounded space of factor-base constructions; the immediate search is over the specified contiguous-orbit bases and index variants.
