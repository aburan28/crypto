# n37 native online wall screen: both bases remain live

The [frozen protocol](PROTOCOL.md) asked whether the K16 target-only instruction
lead survives as native wall time on new public points. On 16 one-target n37
Koblitz workloads, every one of the 384 executions recovered and verified its
scalar. The 320 measured executions were independently replayed on a different
environment class. The hosted Linux run earned **L1**, below the required L2
isolation level, and its identical K16 control had a **12.85%** maximum paired
online deviation, above the frozen **5%** noise gate. The preregistered decision
is `carry_both_host_noise_exceeds_gate`. The admitted one-target
`rho_online_ms / IC_online_ms` speedup is **unknown**.

Target index zero, `Wd281dbafa0cb`, was designated the primary one-target
workload before execution. Each row below is the mean of five paired measured
rounds on that **same public point**; all online ratios are descriptive L1
host-screen diagnostics. The online interval includes every target-dependent
attempt and verification, while reusable preparation stays outside it.

| Arm | Exact method or candidate | Usable base points | Folded columns | Verified rounds | Mean online wall (ms) | Rho / IC online |
|:--|:--|--:|--:|--:|--:|--:|
| Strong signed-Frobenius rho | `ECM1h1d9961ee4601` | — | — | 5/5 | 0.821473 | reference |
| IC K8 | `IC1N37Ckb0fb592PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0hd7ebd078ae0c` | 592 | 8 | 5/5 | 1.224381 | 0.670929, exploratory |
| IC K16 | `IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0hbbfdf029e5a2` | 1,184 | 16 | 5/5 | 0.114444 | 7.177938, exploratory |
| Identical K16 control | same K16 candidate | 1,184 | 16 | 5/5 | 0.114470 | A/A control |

The other 15 frozen one-target workloads are separate replications, not a
batch DLP or a replacement for target zero. Their individual rows and all
five rounds are in [DECISION.json](DECISION.json). Across 16 targets and 80
measured rounds per arm, the preregistered secondary ratio of summed online
nanoseconds and target-block bootstrap intervals were:

| Same-point online comparison | Descriptive ratio of sums | Target-block 95% interval |
|:--|--:|:--|
| K8 / K16 | 7.1103 | [3.6693, 16.9818] |
| Rho / K16 | 3.5462 | [2.0282, 8.0584] |
| Rho / K8 | 0.4987 | [0.3621, 0.7246] |
| K16 / identical K16 | 0.9991 | [0.9940, 1.0024] |

The bootstrap interval does not cure host contention: the largest *single
pair* A/A deviation was 12.85%. Host preflight recorded CPU PSI some avg10
20.37, above its 5.00 limit, and the run used `--allow-busy`. Every measured
row remained L1. On the preregistered target zero, the complete cold solve
means were 1.097646 ms for rho, 9.889181 ms for K8, and 10.172437 ms for
K16: descriptive cold IC/rho is 9.0094 for K8 and 9.2675 for K16. Across all
16 separate targets, complete cold solve means were 1.058 ms for rho, 10.210
ms for K8, and 10.278 ms for K16. These cold wall figures are also exploratory;
they are reported because the reusable work dominates this implementation's
solve, not as a controlled crossover claim.

| Native IC stage diagnostic, mean or fixed count | K8 | K16 |
|:--|--:|--:|
| Actual usable base points / folded columns | 592 / 8 | 1,184 / 16 |
| Pair-table entries | 2,344 | 9,412 |
| Target-blind rank trials / verified useful rows | 40 / 8 | 25 / 18 |
| Mean target PDP attempts, including failed attempts | 8.450 | 1.325 |
| Mean base construction wall (ms) | 0.326 | 0.604 |
| Mean pair-table setup wall (ms) | 1.167 | 3.920 |
| Mean relation/rank preparation wall (ms) | 6.907 | 4.993 |
| Mean complete cold solve wall (ms) | 10.210 | 10.278 |
| Maximum recorded RSS (KiB) | 8,248 | 8,612 |

The [sealed session](sessions/hosted_ubuntu_01) retains all 384 records, raw failures (none),
phase timings, complete resource envelopes, target seeds, source and binary
hashes, point and scalar checks, and the CPU isolation record. The hosted
[local audit](LOCAL-AUDIT.json) reproduced all 320 measured children. The
independent macOS ARM64 [receipt](INDEPENDENT-AUDIT-MAC.json) reproduced the
same 320 children exactly from a different environment class; the
[evidence manifest](EVIDENCE.json) pins their hashes and the [hosted CI run](https://github.com/aburan28/crypto/actions/runs/37201442645).
The native [Rust analyzer](../../examples/ecbench_n37_native_online_wall_analyze.rs)
checks the frozen source/plan, target identities, phase sums, scalar
verification, equal-point pairing, independent replay and decision gate, then
reproduces `DECISION.json` byte for byte. From the repository root:

```sh
cd research/ecbench_n37_native_online_wall_20261004
shasum -a 256 --check SHA256SUMS
cd ../..
cargo build --release --example ecbench_n37_native_online_wall_analyze
target/release/examples/ecbench_n37_native_online_wall_analyze \
  research/ecbench_n37_native_online_wall_20261004 \
  research/ecbench_n37_native_online_wall_20261004/sessions/hosted_ubuntu_01 \
  research/ecbench_n37_native_online_wall_20261004/INDEPENDENT-AUDIT-MAC.json \
  /tmp/n37-native-wall-decision.json
cmp /tmp/n37-native-wall-decision.json \
  research/ecbench_n37_native_online_wall_20261004/DECISION.json
```

This screen validates native correctness and exposes a sizeable K16 online
lead worth testing, but **does not select K16 under the frozen wall gate**.
The next experiment must repeat the same one-target comparison on an auditable
L2 host without relaxing the noise rule. Both K8 and K16 remain controls in
n41/n53; no ECC2K-130 transfer or isolated online speedup follows from this
n37 hosted screen.
