# n83 five-summand degree-four source prolongation: column cap before pruning

The [preregistered protocol](PROTOCOL.md) used one selective Boolean
Macaulay prolongation of the exact dimension-18 five-summand system in
#1417. The planted `[0,2,4,6,8]` assignment satisfied every original
and added equation for every executed `k`. The public ordinary T001
source preimage at torsion offset zero produced no contradiction or
source-only affine row at either completed reduction. The `k=8` row
reached the frozen 1,500,000-column cap; `k=16` was therefore not run.
This is a **bounded structural negative**, not a proof that T001 has no
five-summand decomposition or that a higher-degree F6 method cannot
work.

The registered K0 curve is
`icv1-f2m83-tm6151469093347-debefd74`. The exact standard
dimension-18 base has 261,447 geometric points, 261,444 distinct
subgroup-usable points and 130,722 sign-folded columns. The ordinary
input was the same public T001 subgroup point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)` as #1417. The
five-summand system has 339 Boolean variables, 90 of which are source
coordinate bits, and 332 original cubic equations. This experiment
multiplied each **original** equation once by the first `k` variables
of the frozen source sequence; it did not recursively multiply added
rows. Boolean idempotence and XOR cancellation were exact.

| Source multipliers `k` | Equations | Term occurrences | Distinct columns | Rank | Affine rows | Source-only rows / contradiction | Reduction word XORs | Status |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |
| 1 | 664 | 2,209,164 | 596,027 | 664 | 1 | 0 / no | 107,255,002 | complete |
| 4 | 1,660 | 5,513,470 | 1,463,392 | 1,660 | 1 | 0 / no | 622,318,443 | complete |
| 8 | 2,988 | 9,922,598 | 1,500,001 observed | unknown | unknown | unknown | unknown | column cap |
| 16 | — | — | — | — | — | — | — | not run after cap |

The one affine row at `k=1` and `k=4` is exactly the inherited
`v72 ⊕ v256 ⊕ v337 = 0`: one source bit and two free intermediate
bits. It can be satisfied for either source-bit value, so it prunes no
factor-base point. Full rank in the two completed matrices means the
selected added rows produced no further linear dependency at this
degree. `k=8` is inconclusive because the column cap stopped index
construction before elimination. Since the variable sequence is nested,
the `k=16` column set contains the `k=8` set and would hit the same cap.

All six executed processes exited zero. The planted controls at
`k=1,4,8` reported both `all_equations_zero: true` and
`group_sum_replayed: true`; their source sum and cofactor/torsion
certificate are inherited from the exact #1417 control. The focused
release test `selective_source_prolongation_preserves_planted_zero`
passed. Every raw JSONL, stderr and exit status is retained here. The
largest reported peak RSS was 1,205,534,720 bytes at `k=4`, below the
7-GiB observation gate. The `k=1` and `k=4` reduction intervals were
1.217 and 5.233 seconds; these are **contended feasibility timings**
on an Apple M4 Pro, not speedup evidence.

The release binary SHA-256 was
`7a3332fd73e530c5c110f6e0df918ece891f727da5af5b4f926e2aa1f0f14556`.
The helper source SHA-256 was
`999e5a5e4240d69c422c37c3b9a81cdec920b9dd341b7a2b18647b684e060d6f`;
the probe source SHA-256 was
`e43479be234ddf034a1056bf2e84c17527cf341d118c01d753296e7ea2a72a6f`.
The [release build](build.log), [focused test](helper_test.log),
[runner](run.sh), [status table](status.tsv), raw outputs and
[`SHA256SUMS`](SHA256SUMS) provide the reproducibility receipt. A
future gate needs a different column representation or a better
selection of prolongation rows before increasing the degree-four cap.

No ordinary full-group relation was recovered, so there is no complete
F6 point-decomposition result, factor-base relation yield, one-target
IC online interval, target scalar recovery or same-point rho reference.
The IC candidate and speedup remain unset (`candidate_id: null`,
`IC_online_ms: null`, `rho_online_ms: null`). The nominal five-multiset
capacity of 4.21016 per subgroup element from #1417 is a counting
bound, not a measured natural relation rate.
