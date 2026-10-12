# Fused sigma batch geometry: frozen screen

The counter-free fused sigma schedule removes one state load pass. This screen
tests whether that changes the old batch-16 optimum by amortizing inversion over
32 or 64 slots while keeping the number of live scalar walks and completed
updates identical.

This is a bounded screen, not a promotion run. It keeps the sigma walk,
arithmetic, witness-off mode and one RTX PRO 6000 fixed. No search, solver,
collision recovery or multi-GPU aggregation is run.

## Arms and equal work

All arms use `SIGMA_FUSED=1`, `WITNESS=0`, inline polynomial products and the
same arithmetic as `HEADLINE-PROTOCOL.md`.

| arm | batch | block threads | min blocks | benchmark workers | live walks |
|---|---:|---:|---:|---:|---:|
| b16t256 | 16 | 256 | 2 | 385,024 | 6,160,384 |
| b32t256 | 32 | 256 | 2 | 192,512 | 6,160,384 |
| b32t512 | 32 | 512 | 1 | 192,512 | 6,160,384 |
| b64t256 | 64 | 256 | 2 | 96,256 | 6,160,384 |
| b64t512 | 64 | 512 | 1 | 96,256 | 6,160,384 |

Each timing row uses 1,024 steps and 32 launches, exactly 201,863,462,912
complete scalar updates.

## Correctness and decision

For correctness, every arm runs the same 1,540,096 seeded walks through seven
odd 95-step launches at DP weight 48. Each must replay 300/300 reports with
zero drops and produce the identical sorted headerless-v1 corpus. The native
comparator retains one canonical sorted stream and its SHA-256; the five
redundant unsorted corpora are removed only after identity passes.

One warmup per arm is excluded. Three screen rounds use forward, reverse and
forward arm order. The native summarizer requires three valid equal-work rows
per arm. A non-baseline arm qualifies for a confirmation experiment only when
its median is at least 1.01 times the maximum b16t256 screen sample and every
one of its samples exceeds the minimum b16t256 sample. Otherwise batch 16 is
retained. No screen-only row is promoted into a preset.

