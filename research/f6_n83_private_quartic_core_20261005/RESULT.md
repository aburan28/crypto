# n83 five-summand private degree-four certificate: a smaller unresolved core

The [preregistered first gate](PROTOCOL.md) completed on the exact K0
five-summand system. A two-pass scan selected 32 exact degree-four
monomials per prolonged row and counted each selected monomial in
**every** prolonged row. A count of one certifies that its owner cannot
participate in any linear combination whose degree-four terms cancel.
The filter is exact: a missed private monomial only leaves a row
unresolved. It does not create a false certificate.

The registered curve is `icv1-f2m83-tm6151469093347-debefd74`.
The standard dimension-18 base has 261,447 geometric points, 261,444
distinct subgroup-usable points, and 130,722 sign-folded columns.
The ordinary input was the public T001 subgroup point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)` at torsion
offset zero. The five-summand Boolean system has 332 original cubic
equations, 339 variables and 90 source bits. The frozen variable order
and row IDs are in the protocol.

| Source multipliers | Prolonged rows | Rows with degree four | Degree-four term occurrences | Candidate columns | Rows certified private | Rows unresolved | Peak RSS (MB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 5,312 | 3,984 | 14,131,299 | 119,343 | 3,384 | 1,928 | 356.254 |
| 90 | 29,880 | 22,410 | 79,326,072 | 666,283 | 14,058 | 15,822 | 491.651 |

At `k=16`, the 1,928 unresolved rows are below the protocol's
2,000-row core threshold, but that first protocol only authorized core
reduction for **the full `k=90` certificate**. At `k=90`, 15,822 rows
remain unresolved, so no core reduction was run. The 7,470 `k=90`
rows without degree-four terms necessarily remain unresolved by a
degree-four certificate; other unresolved rows may have private terms
outside the 32 selected candidates or may genuinely share all such
terms. The certificate alone does not decide either case.

All four processes exited zero. The planted `[0,2,4,6,8]` control
checked all 5,312 and 29,880 generated products at `k=16` and `k=90`,
respectively, with every equation zero and the full curve sum replayed.
The focused synthetic test distinguished private and shared degree-four
terms and passed. The ordinary certificate witness digests (BLAKE3 over
ordered row IDs and full monomials) are
`ba12cc6bb138b135f514f6bb6198db21286101b4a55eeb6f42a5af12e4a49d74`
for `k=16` and
`153342de74d4e40a2ea5f96ba35d24b7bdf2726700ed9b04714c9c89ca49dbdd`
for `k=90`. Every unresolved row ID is in the raw JSONL.

The certificate intervals were 1.507 and 9.902 seconds on a contended
Apple M4 Pro. They are feasibility diagnostics, not controlled CPU
speedups. The release executable SHA-256 was
`77807a743bd1ffd72f5c57775c9e68630cd542350932bb774235e7c3e3e15a68`.
The helper source SHA-256 was
`26b681034a6bed583ccd23423c1b91abd55b444787984b72338acb0d1a80d186`;
the probe source SHA-256 was
`d57c893a8abd8d08733a6afe67d102c519d2669ac105167d1e35fd4e77b54660`.
The [build log](build.log.gz), [focused test](helper_test.log.gz),
[runner](run.sh), [status](status.tsv), planted and ordinary JSONL,
stderr, and [SHA-256 manifest](SHA256SUMS) preserve the evidence.

## Amendment 1: exact k=16 core reduction

[Amendment 1](AMENDMENT_1.md) was committed before the new probe was
built or run. Its frozen guard reproduced the 1,928 unresolved rows
and the original `k=16` witness digest. The 1,928 rows plus all 332
original equations produced a 2,260-row core with 6,376,371 term
occurrences. The ordinary reducer stopped at its registered
1,500,000-column cap (`column_limit` at 1,500,001), before rank or
affine consequences could be computed. Peak observed RSS was
746,438,656 bytes. Both the repeated planted control and ordinary
process exited zero; the ordinary cap is an inconclusive structural
outcome, not evidence that no source constraint exists.

The amended probe source SHA-256 was
`775dd89f1c060d08b403871117750e4d94d2ace3fcc96d4643ff23240dec6f7c`;
the frozen release binary SHA-256 was
`37326e268752341ebac1d361cd7614336b7e039de52f6b20e882bc46a8def441`.
The [amended runner](run_core16.sh), [status](core16_status.tsv),
raw `core16_*` JSONL/stderr and build log preserve the new run.

## Amendment 2: the same core completes

[Amendment 2](AMENDMENT_2.md) was committed before adding the explicit
6,500,000-column cap to the reducer. The default cap remains 1,500,000.
The repeated private certificate passed the frozen row-count and
witness-digest guards. The same 2,260-row, 6,376,371-term core had
**1,881,674 distinct columns** and reduced to **rank 2,260**. It left
one affine row, the inherited mixed intermediate/source row, and zero
nonconstant source-only rows or contradictions. Thus this particular
`k=16`, one-round source-variable prolongation did not improve source
pruning on ordinary T001 offset zero. It does not rule out more
multipliers, more rounds, or a different solver.

The planted control and ordinary process both exited zero. The live
RSS sampler reached 1,834,272 KiB, below the registered 7-GiB kill
threshold; the process-reported peak was 1,886,830,592 bytes. The
reduction counted 702,364,047 64-bit XOR operations. Its 1.427-second
interval is a contended-host diagnostic, not a CPU speedup. The
amended helper source SHA-256 was
`9ea54011648e6f3eadb6b18b20b46cf09ce6bf8e573b9d77dc64668b81a50eff`;
probe source SHA-256 was
`ae2057b0c9c6ee5380177cf0f7494b3bdcfb7705774c132fcf9ddf3c9ac73b86`;
frozen binary SHA-256 was
`0f661231c038e7edeff566e172debb817914c2325107fbb22885c23ba3808fa4`.
The [runner](run_core16wide.sh), [RSS samples](core16wide_rss.tsv),
[status](core16wide_status.tsv), build log, and raw JSONL/stderr retain
both processes. This is a bounded negative structural result.

This result narrows the exact degree-four search space but has no
ordinary full-group relation, complete F6 decomposition, one-target IC
online interval, target scalar recovery, or same-point rho reference.
The complete candidate ID and speedup remain unknown.
