#set page(paper: "a4", margin: 16mm)
#set text(font: "Arial", size: 8.5pt)
#set heading(numbering: none)
= First native CPU timing campaign
2026-10-15 · Run 38766895864 · AMD EPYC 7763

Eight BMSS methods, four field widths (192/224/256/511), degrees 7/31: 64 timed cells. Five A/A pairs and five interleaved A/B pairs, five constructor samples per arm per round. All twenty AMD commands pass native guest-isolation checks and independent artifact audit. Full 96-case coefficient comparisons pass on all four CPU profiles at degrees 7/31/101.

== 224-bit fixture / degree 31
#table(columns: (1.7fr, 0.8fr, 1.3fr, 1.3fr), inset: 3pt, [Method], [Median ratio], [Five-round range], [A/A range],
[linear_algebra], [16.919], [16.859–16.962], [0.995–1.003],
[stark1972], [1.053], [1.046–1.067], [0.977–1.003],
[atkin1992], [1.197], [1.182–1.201], [0.996–0.998],
[atkin_modcomp], [1.116], [1.111–1.152], [0.996–1.007],
[elkies1992], [1.303], [1.282–1.324], [0.979–1.002],
[elkies1998], [1.105], [1.089–1.181], [0.985–1.022],
[fast_elkies], [1.074], [1.065–1.190], [0.998–1.044],
[fast_elkies_prime], [1.126], [1.109–1.151], [0.997–1.005],
)
#image("RATIOS.svg", width: 100%)
#image("PIPELINE.svg", width: 100%)

Constructor-only measurements on synthetic supersingular field-width fixtures. Preparing sigma, fixture/reference generation, fingerprints and verification are excluded. These are not named P192/P224 whole-search or replay speedups. The AMD guest has four logical CPUs; the timed child is reserved and pinned by the native tool. Guest contention checks pass for all twenty commands. VM host neighbours and frequency remain unobserved. A/A envelopes and min/max ranges are descriptive, not confidence intervals.

All 64 cells remain in RESULTS.csv and the raw analysis. At 256 bits/degree 7, Atkin 1992 has a paired median ratio 0.984707 and Atkin modular composition 0.970765; their A/A ranges are 0.994689–1.001896 and 0.993324–1.016337. These slower rows are retained. Several other small-case changes are inside the observed A/A envelope. No blanket per-method or cross-CPU win is inferred.

The ARM first A/A run was refused before execution: runner housekeeping consumed 0.29 CPU-seconds in a two-second settling window (limit 0.20), with PSI below its limit. The producer completed nineteen subsequent admitted runs. The requested complete campaign remains censored; no partial ARM timing ratio is promoted. Both macOS profiles have exact repeated counts and outputs but no admitted wall samples.

Baseline 79006eb85a40 → candidate a21306ca2041. Raw ZIPs, all 64 rows, medians, source/binary/isolation audits and failures remain in the linked evidence. Analysis SHA256:
#text(size: 6.5pt)[c483b4b738910f7a62b436b155acc140d76d82c0f7e1f39bcc34b55169be2c4b]
