#set page(width: 210mm, height: 297mm, margin: (x: 17mm, y: 16mm))
#set text(font: "Helvetica Neue", size: 9pt)
#set heading(numbering: none)

= N83 current-source control: run record

2026-10-17 · one previously exposed public point · shared macOS host

K1280/root11 dynamic recovered the public point in 4.73325 ms online with
212,480 distinct subgroup-usable factor-base points, 1,280 folded columns,
and 1,342,177,280 bytes of packed root allocation. Its receipt wall was
146.960 s. Independent arithmetic replay passed for all eight current-source
controls in the accompanying results report.

The exact curve is `EC1N83Ckb1a1haacb69ee9f37`, over GF(2^83) with
defining-polynomial integer `9671406556917033397649557`, equation
`y² + xy = x³ + x² + 1`, and prime subgroup order `8569786107849059`.
The previously exposed public point came from the frozen 2026-10-16 cold
factor-base fixture; the known-answer scalar was used for independent replay.

The source SHA-256 was
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
The local release binary SHA-256 was
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
The K1248 and K1280 processes ran sequentially after K2816 exited, each
with a 600 s cap; their receipt walls were 143.224 s and 146.960 s.
No local CPU isolation admission was recorded. Each online interval begins
after factor-base, index and log preparation; receipt wall includes the
fresh producer process.

The raw small-base archive is pinned by SHA-256
`f2b9cbc454705add4b8691785fefc8a4498764cc8066e0e6e294ce931f16f685`;
the K2816 archive by
`61daed9917be546f46bf72cc8d634525fefbce5a27a01cf72fd56fc5f6cfde63`.
Versioned independent checkers reconstruct every base point, sum each
four-point witness, replay the scalar, check online phases, and derive the
exact candidate identity. The archives, receipts, replay summaries, protocols
and checker source are committed in the study directories.

This run is a candidate-selection diagnostic. The next experiment will
compare K1248, K1280, K2543, K3072 and the selected finalist on ten
previously unused public points, in alternating Tier 1 same-point IC–rho
pairs.

Source report: `REPORT.md`; plotted rows and record paths: `plot_data.json`.
