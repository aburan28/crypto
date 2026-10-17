Research cycle 2026-10-16 v1

Hypothesis: reuse the complete three-source table across targets while preserving all original eight-coordinate roots; measure squarefree support cost rather than assuming inversions are cheap.

- batch-separator-manyfield-exact-coldpairs-20261016-v1: **passed**, worker cold 37.082s. Shared complete leftseparator and certified-edge joins matched every original root against prior cached separator across the manyfield panel.
  Whole-cold gain against strongest prior cached separator: 1.3830x; 20 fields, 38 contexts; context geometric gain 2.4514x.

| Case | Median old/new gain |
|---|---:|
| prime-3 | 3.7430x |
| prime-5 | 8.4530x |
| prime-7 | 5.9313x |
| prime-11 | 7.8510x |
| P192 | 1.9955x |
| P224 | 0.9990x |
| P256 | 1.6805x |
| P384 | 1.4667x |
| P521 | 1.1792x |
| Koblitz-2-a0 | 5.9968x |
| Koblitz-2-a1 | 7.1977x |
| binary-2 | 7.2352x |
| Koblitz-3-a0 | 5.9499x |
| Koblitz-3-a1 | 5.2524x |
| binary-3 | 5.6515x |
| Koblitz-4-a0 | 6.0375x |
| Koblitz-4-a1 | 6.1011x |
| binary-4 | 5.0914x |
| Koblitz-83-a0 | 1.3707x |
| Koblitz-83-a1 | 1.4539x |
| binary-83 | 1.3714x |
| Koblitz-163-a0 | 1.6231x |
| Koblitz-163-a1 | 1.5858x |
| binary-163 | 1.3285x |
| Koblitz-233-a0 | 1.4415x |
| Koblitz-233-a1 | 1.4875x |
| binary-233 | 1.2234x |
| Koblitz-283-a0 | 1.5528x |
| Koblitz-283-a1 | 1.6224x |
| binary-283 | 1.1840x |
| Koblitz-409-a0 | 1.1907x |
| Koblitz-409-a1 | 1.0876x |
| binary-409 | 1.2289x |
| Koblitz-571-a0 | 1.1356x |
| Koblitz-571-a1 | 1.1261x |
| binary-571 | 1.1321x |
| oddext-3-2 | 6.2141x |
| oddext-5-2 | 6.7883x |

- support-squarefree-cost-manyfields-20261016-v1: **passed**, worker cold 5.783s. Exact nilpotence and squarefree-shortcut support criteria compared with all repeated roots retained.

| Case | Median old/new gain |
|---|---:|
| prime-11 | 0.9338x |
| P192 | 0.7034x |
| P224 | 0.6788x |
| P256 | 0.5948x |
| P384 | 0.6633x |
| P521 | 0.6131x |
| Koblitz-3-a0 | 0.8071x |
| Koblitz-4-a0 | 0.9266x |
| Koblitz-83-a0 | 1.0060x |
| Koblitz-163-a0 | 1.0077x |
| Koblitz-233-a0 | 1.0157x |
| Koblitz-283-a0 | 1.0563x |
| Koblitz-409-a0 | 1.0319x |
| Koblitz-571-a0 | 1.0270x |
| oddext-3-2 | 0.5450x |
| oddext-5-2 | 0.5463x |

- support-squarefree-cost-actual83-20261016-v1: **passed**, worker cold 1.472s. Exact nilpotence and squarefree-shortcut support criteria compared with all repeated roots retained.

| Case | Median old/new gain |
|---|---:|
| actual83-fixedh-k-relaxed16 | 1.1409x |

The 15 separator and 4,448 support oracle cases passed, including repeated and inseparable quintics in characteristics2,3,5. Sole-model review with independent oracle algorithms is not an external independent review. Full-root output equality is checked in separator trials; support trials classify supplied quintics and do not search or count original roots.

The actual GF83 domain2^18 has a left-table bound2^56, explicitly refused by preflight. These measurements do not establish hard-system gains over genuine F4/F5/F6, or all-scale coverage. The original release gate and additional requested5x target remain unmet; the program stays active. Components are classical; no primary-source novelty claim is made.

Next experiments:

- batch-certified-output-native-arithmetic-followup-20261016: Profile every retained field, then bound unique-edge verification/output cost and test native carryless arithmetic or representation changes against the strongest cached separator. Full original tuples remain mandatory; batch reuse is a tiny-domain calibration, not a 2^18-domain solution.
- support-cost-dispatch-20261016: Use measured squarefree GCD overhead, especially actual GF83, to decide whether to reject blanket promotion. Freeze any cost dispatch as a new variant; compact nilpotence remains exact on repeated roots.
- walsh-six-to-eight-place-rankwidth-20261016: Before adding two F8 places (moduli11,13), profile sparse contraction/rank width of the augmented local-image constraints. The signature grows from96 to120bits with rank at most76: a dense 2^44 transform is forbidden. Budget90sec,4GiB,onecore; return an exact width bound or an explicit refusal, not enumerate the billion surviving parameters.
- hard-matched-groebner-5x-panel: Keep genuine source-identified matched all-roots F4/F5/F6 interface work independent. Require immutable20-case gate,3cold repeats,fastestbaseline>=1s,A/A<=10%,worstslowdown<=1.10,qualified isolation and requested field coverage; user target5x is additional. No current calibration qualifies.
