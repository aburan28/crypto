#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261011T2258Z-glvgls-lattice]

*What ran.* `cargo run --release --example isogeny_loop_glv_gls_lattice -- --bits 127 --samples 10000 --seed 20261011 --max-tries 2000000 --discs -3,-4,-7,-8,-11,-23,-71,-431,-619,-3019 --json research/isogeny_loop_glv_20261011/glvgls_lattice.json`. Base commit db0393b4 plus the example source with SHA-256 `e789d26c…4c01` (`../../revision.txt`); rustc 1.97.0; Linux 6.18 x86_64 cloud container, not isolated. Start 2026-10-11T22:58Z (`../../glvgls_run_start_utc.txt`); wall time 7.6 s (`time`, recorded in `../../glvgls_lattice.log`). Seed 20261011 (per-discriminant streams derived from it).

*Arms.* Per discriminant: certified cheapest loop (exhaustive to the norm bound); GLS-2 lattice; GLV-GLS-4 lattice for both ring embeddings ($s$, $-s$). 10,000 seeded scalars per lattice.

*Runs requested / completed.* 10 discriminants requested; 9 produced an instance; $D = -4$ produced none within $2 dot 10^6$ draws (recorded; the cause is derived in the study README §3). *Earlier runs:* a 22:51Z run of the same source before formatting and clippy-only edits produced identical rows (checked field by field). *Failed earlier attempt, kept:* the first full run, without the draw cap, was stopped by hand after 13.5 CPU-minutes stuck in the $D = -4$ search; it wrote no records. Checks asserted per instance: $mu^2 = -1$, $mu^2 - t mu + p = 0$, $s^2 = D$, $lambda^2 - "Tr"(alpha) lambda + N(alpha) = 0 space (mod r)$; every LLL vector in its lattice; every decomposition congruent to its scalar. All passed (exit 0).

*Recorded values* (unit: bits, max / mean over 10,000 scalars of $max_i "bits"(|k_i|)$): GLS-2 124–126 / 122.8–125.5; GLV-GLS-4 62–63 / 61.2–62.6, against $log_2 r slash 4$ = 62.0–63.5. Cheapest loop costs (M, projective): 1, 1, 15, 15, 22.5, 37.5, 60, 97.5, 127.5, 210. Model ratio GLS-2 / GLV-GLS-4 (derived operation count, not measured): 1.158–1.357. Full table: `../../glvgls_lattice.md`; record: `../../glvgls_lattice.json`.

#figure(image("../../figures/glvgls_bits.svg", width: 80%), caption: [Recorded mean coefficient sizes per discriminant (run run-20261011T2258Z-glvgls-lattice).])

*Observations only.* Interpretation is in `../../README.md` §3 and `../../REPORT.pdf`.
