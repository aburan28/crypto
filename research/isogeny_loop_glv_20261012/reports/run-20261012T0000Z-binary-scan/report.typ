#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261012T0000Z-binary-scan]

*What ran.* `./target/release/examples/binary_gls_scan --m M --dmax X --keep 6 --samples 5000 --lmax 61 --json …` with these parameters:

#table(columns: 3,
  [*M*], [*X*], [*start*],
  [127], [4,000,000], [00:00:27Z],
  [83], [8,000,000], [00:00:37Z],
  [131], [8,000,000], [00:00:46Z],
)

Final source, before the checkpoint commit after 8b58692c5. rustc 1.97.0, release, cloud container, not isolated, seed 20261012. Wall times: 5.2 s, 9.6 s and 10.2 s (`../../binary_scan_m*.log`).

*Method.* For each fundamental D ≡ 1 (mod 8) the scan requires:
- the class of 𝔭#sub[2] has order exactly m;
- Cornacchia gives 4·2#super[m] = t² + |D|v² with t odd;
- r = ((2#super[m] − 1)² + t²)/2 is a probable prime.

It then searches odd prime ℓ ≤ 61 for [𝔩][𝔭#sub[2]]#super[k] = 1. It keeps F-oriented α = (A + B√D)/2 of norm ℓ·2#super[min(k, m−k)] and ranks them by model cost 7.5ℓ + 2·(F steps) M. GLS-2 and GLV-GLS-4 decompositions are measured on 5,000 seeded scalars.

*Runs requested / completed.* 3 / 3, exit 0. An earlier m = 127 run of the same scan printed "kept 2 with prime r" while listing 3 prime-r rows; the wording conflated two counts. The summary line was corrected to report both counts, and the run was repeated. Only the repeat is recorded.

#table(columns: 6,
  [*m*], [*scanned*], [*ord = m*], [*prime r*], [*with loop*], [*cheapest (D: ℓ, F steps, M)*],
  [83], [810,576], [440], [4], [0], [—],
  [127], [405,292], [319], [3], [2], [−1532327: 3, 53, 128.5],
  [131], [810,576], [459], [4], [3], [−283583: 13, 16, 129.5],
)
GLV-GLS-4 max coefficient bits: 64 (m = 127) and 66 (m = 131), against GLS-2 at 126 and 130.

#figure(image("../../figures/binary_scan.svg", width: 95%), caption: [Recorded funnel counts per m (log scale), and the cheapest loop's model cost per kept discriminant (run-20261012T0000Z-binary-scan).])

*Observations only.* Interpretation: `../../README.md` §4.
