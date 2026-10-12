# Direct six-summand n83 F6 system: admitted, root gate inconclusive

The [preregistered protocol](PROTOCOL.md) and its
[memory-limit](AMENDMENT_1.md), [column-cap](AMENDMENT_2.md), and
[affine-row](AMENDMENT_3.md) amendments preceded their respective runs.
This branch is stacked on PR #1409. The native 512-bit monomial path
constructs a *coupled* six-summand S3 chain on public K0 curve
`icv1-f2m83-tm6151469093347-debefd74` (field modulus
`z^83+z^45+z^2+z+1`, exact subgroup order
`2417851639230796216685689`, cofactor 4). The dimension-16 standard
source subspace has 64,907 geometric points and its frozen projected
inventory has 64,904 usable points in 32,452 signed columns. The direct
chain has 428 Boolean variables (96 source-coordinate and 332
intermediate bits) and 415 cubic coordinate equations.

## Exactness and bounded ordinary query

The final binary replayed source indices `[0,2,4,6,8,10]`: all 415
Boolean equations evaluated to zero; six full curve points summed to
`(1df49627af40feb6fa421,7da13c02c90c95ffe206)`. Cofactor-four
projection gave
`(6c3da8931c913ae9934bb,15ea852eb0d5db8fc71bb)`; the subgroup
preimage and exactly one of four verified rational 4-torsion offsets
(offset 3) reconstructed the source sum. This is a planted correctness
control, **not ordinary relation yield**. The focused 512-bit boundary
and exact affine-row tests passed 2/2. The final planted [raw
row](planted_final.jsonl) and [status](planted_final.status) are retained.

The ordinary query was the previously pinned public target T001,
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`. Its subgroup
preimage and all four rational torsion offsets each constructed 415
equations, with 1,402,606–1,412,608 term occurrences. The offset-0
own-degree reduction had **384,895 distinct monomial columns, rank
415, no contradiction, and one affine row**:

`v80 ⊕ v345 ⊕ v426 = 0`.

Here `v80` is a sixth-summand coordinate and `v345,v426` are free
intermediate coordinates. This single row can be satisfied for either
value of `v80`; it does not itself prune any source-summand choice.
The reduction produced no witness. It does not prove that T001 has no
six-summand decomposition, nor that stronger algebraic reductions
cannot find one. The ordinary relation gate therefore **fails to
establish a usable solver** under this protocol. The direct-system
admission succeeds, but its own-degree root step is insufficient.

| Run | Input and revision | Outcome | Result receipt |
| --- | --- | --- | --- |
| `planted_1` | Planted, original source | Runner failed before binary output: sandbox denied `ps` | Empty [JSONL](planted_1.jsonl), no status |
| `planted_2` | Planted, original source | Exact equations and group/cofactor replay pass | [JSONL](planted_2.jsonl), [status](planted_2.status) |
| `ordinary_0` | T001 offset 0, 300,000-column cap | Stopped at column 300,001; inconclusive | [JSONL](ordinary_0.jsonl), [status](ordinary_0.status) |
| `ordinary_1` and `ordinary_2` | T001 offsets 1 and 2 | Both constructed; reduction not attempted | [JSONL 1](ordinary_1.jsonl), [JSONL 2](ordinary_2.jsonl) |
| `ordinary_3` | T001 offset 3, original runner | Construction row written; runner lost exit status in `ps`/`pipefail` race | [JSONL](ordinary_3.jsonl), no status |
| `ordinary_3_retry` | T001 offset 3, fixed runner | Construction complete, exit 0 | [JSONL](ordinary_3_retry.jsonl), [status](ordinary_3_retry.status) |
| `ordinary_0_extended` | T001 offset 0, 1,500,000-column cap | Root rank 415, one linear row, no contradiction | [JSONL](ordinary_0_extended.jsonl), [status](ordinary_0_extended.status) |
| `ordinary_0_affine` | T001 offset 0, exact row receipt | Same rank and column count; row support `[80,345,426]` | [JSONL](ordinary_0_affine.jsonl), [status](ordinary_0_affine.status) |
| `planted_final` | Planted, final binary | Exact replay pass, exit 0 | [JSONL](planted_final.jsonl), [status](planted_final.status) |

All stderr files, including the empty ones, and every raw/status SHA-256
are in [SHA256SUMS](SHA256SUMS). The failed/partial rows are retained and
are not successes. The ordinary affine-row run's reported peak RSS was
329,236,480 bytes; the monitor sampled at most 321,520 KiB. The
120-second time guard and sampled 7-GiB RSS guard were not hit. macOS
rejected a kernel `ulimit -v`/`-d` limit, so the RSS guard is sampled,
not a hard memory bound.

## Reproduction and claim boundary

Source revision for the final replay:
`d3285c0bc413c2ff6e5a5713a36669ebf4afa8b8` (the result files
follow in a later commit). SHA-256: native solver source
`8a21ad401b21378e059f7e48459ffdaffeeb070ba8361937dc63cd10ac21e094`,
probe source
`4c87984a35c2999163fca11bbf3a85892fb5df62ac8eaf0ffc627a09d8749520`,
protocol
`0e838bb3066c609bf6a3c6f2ab988eb91042c0de5ecc053ac1317b0917d6dd3c`,
final native executable
`1314c31268cc2e0ed72ff326d84d21b4355a4af752407a7a41b9b98c8e13f663`.

Host: physical Apple M4 Pro, 14 CPU cores (10 performance, 4
efficiency), 48 GB RAM, macOS Darwin 25.6.0 arm64, Rust 1.93.1.
One process was run at a time. Rebuild and replay:

```sh
CARGO_INCREMENTAL=0 CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target cargo build --offline --release --example f6_n83_sixsum_wide_probe
./research/f6_n83_sixsum_wide_20261005/run.sh ordinary_0_affine ordinary 0
./research/f6_n83_sixsum_wide_20261005/run.sh planted_final planted
```

This host has no auditable exclusive CPU partition. The final planted
setup time ranged from 0.92 seconds in `planted_2` to 11.86 seconds
in `planted_final`, illustrating contention; **no CPU speedup ratio is
reported**. The combinatorial six-sum capacity ratio 42.950379 is an
upper bound before duplicate sums, not a measured success rate. This
experiment is an algebraic feasibility diagnostic, not an IC candidate:
`candidate_id`, online single-target cost, `S`, relation yield/rank,
final sparse linear algebra, scalar recovery, matched rho reference,
and F4/F5/F6 or IC speedup are all **unknown**. No claim is made that
the direct 512-bit path beats F4, F5, or rho.

Decision: stop this root-only direct chain as a search strategy. A
follow-on solver must communicate a restriction on the *summand set*
across the chain before branching over 64,907 source points, and must
return a full-group-verified ordinary relation with all failed attempts
charged before any throughput comparison. This is a negative search
gate and a positive exact system-admission result, not a lower bound on
other F6 algorithms.
