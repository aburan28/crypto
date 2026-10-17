# Per-family external comparison: isogeny_algos vs PARI/GP 2.15.4 (2026-10-17)

`OUT=results/p18-families external/families.sh`: inputs from `examples/ext_families.rs gen`
(`cases.gp`), PARI side `external/pari_families.gp`. Supersingular E: y² = x³ + x over
p = 4ℓⁿk − 1, a rational kernel point of order ℓⁿ. Two rounds, tool order reversed; every run on
one isolated CPU (`tools/isolated_bench.py`): 12 runs, none contended (`isolation.jsonl`).
Times are the best per call over both rounds (ours: wall per call; PARI: CPU time, 1 ms
resolution amortised over repeated calls). Every case gives identical outputs (codomain a4, a6;
x-map at r = 123456789): 26/26, 0 mismatches.

| workload | ours | PARI | ours / PARI |
|:--|:--|:--|--:|
| kernel → codomain (ours: √élu / x-only Vélu, faster of the two; PARI `ellisogeny(E, P, 1)`), ℓ = 101, 1009, 10007, 61–511 bits | 5.0 µs – 2.55 ms | 40 µs – 26.3 ms | 0.027 – 0.263 |
| kernel → rational maps (ours: x-multiples + Kohel; PARI `ellisogeny(E, P)`), ℓ = 101, 1009 | 28 µs – 12.7 ms | 6.3 ms – 4.31 s | 0.0013 – 0.0127 |
| ℓⁿ chain (ours: optimal strategy; PARI: `ellisogeny` + `ellisogenyapply` per step), ℓ = 2, 3, 2ⁿ, 3ⁿ ≈ √p, 127–511 bits | 0.63 – 30.0 ms | 2.45 – 92 ms | 0.224 – 0.338 |

Largest ratio: 0.338 (chain, 511 bits, ℓ = 3, n = 152). No case where PARI is faster.

Per case (µs):

| family | bits | ℓ | n | ours | PARI | ratio |
|:--|--:|--:|--:|--:|--:|--:|
| codomain | 61 | 101 | 1 | 5.0 | 40 | 0.124 |
| codomain | 61 | 1009 | 1 | 26.8 | 425 | 0.063 |
| codomain | 61 | 10007 | 1 | 273.7 | 4500 | 0.061 |
| codomain | 127 | 101 | 1 | 5.1 | 85 | 0.060 |
| codomain | 127 | 1009 | 1 | 24.5 | 925 | 0.027 |
| codomain | 127 | 10007 | 1 | 295.7 | 10250 | 0.029 |
| codomain | 255 | 101 | 1 | 16.4 | 130 | 0.127 |
| codomain | 255 | 1009 | 1 | 78.1 | 1475 | 0.053 |
| codomain | 255 | 10007 | 1 | 754.7 | 16000 | 0.047 |
| codomain | 511 | 101 | 1 | 55.3 | 210 | 0.263 |
| codomain | 511 | 1009 | 1 | 292.4 | 2750 | 0.106 |
| codomain | 511 | 10007 | 1 | 2545.2 | 26250 | 0.097 |
| maps | 61 | 101 | 1 | 28.4 | 6300 | 0.0045 |
| maps | 61 | 1009 | 1 | 1155.7 | 576000 | 0.0020 |
| maps | 127 | 101 | 1 | 61.2 | 12600 | 0.0049 |
| maps | 127 | 1009 | 1 | 1855.1 | 1424000 | 0.0013 |
| maps | 255 | 101 | 1 | 152.4 | 23400 | 0.0065 |
| maps | 255 | 1009 | 1 | 4478.9 | 2149000 | 0.0021 |
| maps | 511 | 101 | 1 | 508.6 | 40000 | 0.0127 |
| maps | 511 | 1009 | 1 | 12680.6 | 4310000 | 0.0029 |
| chain | 127 | 2 | 60 | 817.4 | 3650 | 0.224 |
| chain | 127 | 3 | 38 | 628.9 | 2450 | 0.257 |
| chain | 255 | 2 | 120 | 4816.6 | 14500 | 0.332 |
| chain | 255 | 3 | 76 | 3413.6 | 10750 | 0.318 |
| chain | 511 | 2 | 240 | 30021.2 | 92000 | 0.326 |
| chain | 511 | 3 | 152 | 19577.9 | 58000 | 0.338 |

Not compared, for want of an external tool on this host: Sage (passagemath) does not fit in the
session's disk allowance; PARI has no F_p equivalent of rational ℓ-isogeny enumeration with
kernels, BMSS, duals, the CSIDH action, KLPT, Richelot or isogeny-graph path finding. Earlier
records cover Φ_ℓ and SEA against PARI (`results/p17-sea/`) and polynomial arithmetic against
FLINT (`results/p16-external/`).
