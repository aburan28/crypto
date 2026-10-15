#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
= Run receipt: quadratic-extension torsion-fiber DLP

*2026-10-14 · Darwin 25.6.0 arm64 · shared/L0 · Rust 1.93.1*

The frozen native run used four public instances over $F_(p^2)$, subgroup orders 19, 61, 211 and 811, three rounds per size, and seven complete-cost arms. The fixture SHA-256 is `1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`; the spec SHA-256 is `acab50556e889563b7352799647f027f9554985673504169101196c126021f82`. The V3 solver code is recorded at commit `4596c3ab59` with source SHA-256 `034c05398c4652453f74e2a1eb675dee9c474e496582a0be339dcc9e35608500`.

The V3 record SHA-256 is `ca9a0ec7e31de59c58a78ae5c269350b37963b62ff9de815c805c36c2ec26693`. The independent field/curve arithmetic implementation read those records in a second process and checked 75 recovered scalars by $[k]G=Q$, while preserving nine cap-exhausted raw-pair rows. Its receipt SHA-256 is `785c5e615a18742302f05f6c2a219b260d48d5e2b63e984cec945222bd51ac5c`. Every relation used for recovery passed an exact source group check.

Cold current-thread CPU intervals include line lifting, cofactor projection, table build, unsuccessful queries, rank, recovery and internal verification. The isogeny and six-class rho arms charge map and eigenvalue work respectively. The timing clock stops before serializing the report. Per-arm peak RSS and L2 isolation were not captured. The direct Semaev and $chi$ polynomial oracles remain unimplemented.

At order 811, median cold CPU in locally calibrated group-addition equivalents per $sqrt(r)$ was 202.56 for the projected pair table, 284.96 for the $chi$ index, 217.45 for the isogeny control, 7.30 for negation rho and 9.53 for six-class rho. The $chi$/negation-rho paired total-cost ratio was 51.89. The three-round ranges and all arms appear below.

#image("../../GRAPH.svg", width: 100%)

The first run produced 84 rows but its system-time wrapper failed under sandboxed `sysctl`; the second run included JSON serialization inside its timed interval. Both remain preserved as V1 and V2. V3 supplies the reported comparison. Source paths, exact commands, all raw records, phase counters, finite-support census, the full table and limitations are in `RUN_REPORT.md` and `../../REPORT.md`.
