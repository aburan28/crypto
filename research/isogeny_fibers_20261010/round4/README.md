# Round 4: the §16 leads

Results are summarised in `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §17. Each lead
directory holds `RESULTS.md` (question, exact commands, tables, verdict,
limits), raw outputs, graphs with their plotting data, and `report.pdf`.
`BRIEF.md` is the shared protocol every lead followed.

| lead | directory | example | report |
|:--|:--|:--|:--|
| 16.1 feasibility crossover in m | `16_1_m_crossover/` | `examples/r4_16_1_m_crossover.rs` | `16_1_m_crossover/report.pdf` |
| 16.2 fast 4E subgroup check | `16_2_fast_4e/` | `examples/r4_16_2_fast_4e.rs` | `16_2_fast_4e/report.pdf` |
| 16.3 vertical security margin | `16_3_vertical_margin/` | `examples/r4_16_3_vertical_margin.rs` | `16_3_vertical_margin/report.pdf` |
| 16.4 ζ-fold over F_{p²} | `16_4_zeta_fold_fp2/` | `examples/r4_16_4_zeta_fold_fp2.rs` | `16_4_zeta_fold_fp2/report.pdf` |
| 16.5 class-wide multi-target rho | `16_5_class_multitarget_rho/` | `examples/r4_16_5_class_multitarget_rho.rs` | `16_5_class_multitarget_rho/report.pdf` |
| 16.6 genus-2 saturation law | `16_6_genus2_law/` | `examples/r4_16_6_genus2_census.rs` | `16_6_genus2_law/report.pdf` |
| 16.7 isogeny oracles in the AGM | `16_7_agm_isogeny_oracles/` | `examples/r4_16_7_agm_pushforward.rs` | `16_7_agm_isogeny_oracles/report.pdf` |

Host: a shared 4-core cloud container. Runs were single-threaded and
reported in counted operations; any wall time is labelled contended.
