#set page(paper: "a4", margin: 22mm)
#set text(font: "Libertinus Serif", size: 10pt)
#set heading(numbering: "1.")

= DP32 ordinary-sigma table screen

*Run:* `20261012T1815Z`
*Provider:* RunPod, pod `l1of60e0sm1spt`
*GPU:* NVIDIA RTX PRO 6000 Blackwell Server Edition
*Driver / CUDA:* 580.82.07 / 13.3.73
*Source:* `615381ba5c5fa89236a228cd75338d636fa81903` plus uncommitted experiment changes

== Result

The exact polynomial-coordinate sigma table is correct for the ordinary witness-carrying DP32 walk, but this screen does not justify a production rollout or the 20 billion iterations/s target. The candidate median was *9.254351 B/s*: +0.7947% over the matched 512-thread control, but -1.0330% against the live 256-thread geometry. The live arm was fastest at *9.350950 B/s* on this card.

This is a same-GPU comparison. Existing workstation and separate table-walk numbers are not interchangeable with this campaign result.

#image("dp32-rates.svg", width: 145mm)

== Protocol

All arms use the ordinary sigma walk (`WALK_TABLE=0`, `SIGMA_FUSED=0`, `WITNESS=1`), DP32, 1,024 steps/launch, and 385,024 benchmark threads. The panel used one warmup per arm and three rotating timing pairs; each arm has three recorded samples. The candidate changes only the weighted-prefix sigma evaluation to the exact polynomial-coordinate table.

== Correctness

- CUDA arithmetic suite passed for all arms.
- Each arm produced 5/5 CPU-reference-verified points with zero drops.
- Sorted v2 corpora were identical across live, geometry, and candidate.
- Candidate receipt: 126 registers/thread, 68,352 dynamic shared bytes, one block/SM.

This is a screen, not full campaign qualification. A prior 300-report attempt hit the explicit 45-minute limit while CPU replay was the bottleneck; its raw receipt is preserved with the report.

== Conclusion

The table path is a valid optimization candidate, but it is not the route to 20 B/s on this card and is not being rolled out. Further work should target register/shared-memory pressure or the existing 256-thread geometry.
