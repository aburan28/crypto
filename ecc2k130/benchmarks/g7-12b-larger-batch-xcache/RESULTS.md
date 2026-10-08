# Larger-batch shared-X cache result on G7

The selected B16/cache4/min3 implementation remains production. Increasing the
shared-X cache recovers about 1.5--1.75% inside B32, but B32 still trails B16.
The B24 cache increase is flat. No arm qualifies for confirmation and the
12 billion complete scalar updates/s goal remains unmet.

## Initial screen

Each row contains three interleaved samples of 34,359,607,296 complete scalar
updates on the same AWS g7.2xlarge / RTX PRO 4500 at 165 W.

| Variant | Median B updates/s | Paired speedup vs selected | 95% CI | Rate / 12 B/s | Decision |
|---|---:|---:|---|---:|---|
| selected B16/cache4/min3 | 6.178327 | 1.000000 | [1.000000, 1.000000] | 0.514861 | reference |
| B24/cache4/min3 | 6.013335 | 0.970669 | [0.962303, 0.979108] | 0.501111 | regression |
| B24/cache4/min2 | 5.903623 | 0.950172 | [0.938710, 0.961774] | 0.491969 | regression |
| B24/cache6/min2 | 5.982853 | 0.964314 | [0.946721, 0.982234] | 0.498571 | regression |
| B32/cache4/min3 | 5.909547 | 0.945873 | [0.922701, 0.969628] | 0.492462 | regression |
| B32/cache4/min2 | 5.891344 | 0.942877 | [0.911779, 0.975035] | 0.490945 | regression |
| B32/cache8/min2 | 6.015860 | 0.959373 | [0.924547, 0.995511] | 0.501322 | regression |

The direct cache contrasts separate cache size from batch size and launch bound:

| Contrast | Paired speedup | 95% CI | Conclusion |
|---|---:|---|---|
| B24 cache6 / cache4, min2 | 1.014883 | [0.999140, 1.030875] | unconfirmed |
| B32 cache8 / cache4, min2 | 1.017495 | [1.008476, 1.026596] | local gain |
| B24 min2 / min3, cache4 | 0.978883 | [0.970612, 0.987225] | min2 regresses |
| B32 min2 / min3, cache4 | 0.996832 | [0.979984, 1.013970] | unresolved |

## Follow-up with the selected compiler schedule

The follow-up requested min3 for the larger caches. Shared memory still limits
the candidates to two active blocks per SM, while ptxas returns to the selected
80-register schedule. Each row again contains three interleaved equal-work
samples.

| Variant | Median B updates/s | Paired speedup vs selected | 95% CI | Rate / 12 B/s | Decision |
|---|---:|---:|---|---:|---|
| selected B16/cache4/min3 | 6.208908 | 1.000000 | [1.000000, 1.000000] | 0.517409 | reference |
| B24/cache4/min3 | 6.056122 | 0.970602 | [0.960041, 0.981280] | 0.504677 | regression |
| B24/cache6/min3 | 6.045993 | 0.967811 | [0.951462, 0.984442] | 0.503833 | regression |
| B32/cache4/min3 | 5.950580 | 0.949334 | [0.927919, 0.971242] | 0.495882 | regression |
| B32/cache8/min3 | 6.050341 | 0.963305 | [0.935844, 0.991572] | 0.504195 | regression |

B24 cache6/cache4 is 0.997125 with 95% CI [0.988729, 1.005591]. B32
cache8/cache4 is 1.014717 with 95% CI [1.008190, 1.021287], confirming that the
larger B32 cache itself helps while remaining insufficient to beat B16.

## Correctness and resources

All six new binaries reproduce the selected built-in log. The initial four
variants pass exact normalized checkpoint and DP-multiset comparisons at run
IDs 0 and 139; the two follow-up variants pass the same four comparisons. The
timing-eligible cache6/cache8 variants pass their actual-walk shared-X probes,
memcheck, initcheck and synccheck. The initial actual-walk probes cover 1,088
launches per candidate; all outputs match.

The min2 B24 binaries use 128 registers and the min2 B32 binaries use 126, with
zero spills. Cache6 allocates 38,272 bytes of kernel shared memory; cache8
allocates 46,976 bytes. The min3 follow-ups use 80 registers and an 8-byte
stack/spill load/store allocation, matching the selected register schedule.
Their shared allocations remain 38,272 and 46,976 bytes. The measured SM shared
capacity is 102,400 bytes, so both remain at two resident blocks.

Generic work remains `sqrt(n/262)` with ratio 1. Full-DLP S remains null. This
is an engineering comparison of the complete scalar rho walk.

[Frozen protocol](README.md) · [Follow-up protocol](FOLLOWUP.md) ·
[Initial timings](comparison.json) · [Follow-up timings](comparison-followup.json) ·
[Initial state checks](state-validation.json) · [Follow-up checks](followup-validation.json) ·
[Initial sanitizers](sanitizers.json) · [Follow-up gates](followup-gates.json)
