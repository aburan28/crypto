#set page(paper: "us-letter", margin: 0.65in)
#set text(size: 10pt)

= Bounded structural run: torsion fibers
2026-10-12T19:47Z · work item T-2070 · base commit `e3090cbc98c657106fef10020da8ccc19be49a3f`

The frozen census scanned all 1,014 curves specified by `PROTOCOL.md`:
primes p ≤ 199 with p ≡ 1 mod 3, and all nonsquare b in y² = x³ + b.
Direct point enumeration and quadratic-character counts agreed in every
cell; every trace passed the Hasse check. There were 349 prime-order rows.
The named F₁₃ map check counted 19 rational points on each side, 19
distinct images, and 361/361 ordered-pair homomorphism matches.

The runs used native Rust (`rustc 1.93.1`, Homebrew arm64), single process
in a local interactive workspace. No host isolation or timing inference is
claimed. Both commands exited 0; there were no failed or timed-out cells.
Raw records: `census.tsv`, `summary.tsv`, `map_check.txt` in the parent
study. The source and output hashes are in this run's `README.md`.

#figure(image("../../fiber.svg", width: 100%), caption: [The named map-control fiber.])

#figure(image("../../census.svg", width: 100%), caption: [Exact census populations by base prime.])

The figure compares exact population counts, not algorithm costs. There
was no target-logarithm run, matched rho run, or end-to-end S value.
