#set page(paper: "a4", margin: 15mm)
#set text(font: "Arial", size: 8pt)
= Exploratory native family coverage
2026-10-15 · Three native profiles · 395 records each

#image("COVERAGE.svg", width: 100%)

#table(columns: (1.5fr,1fr,1fr,1fr), inset: 2pt, [Reported group], [AMD true / rows], [ARM true / rows], [M4 true / rows],
[big], [82 / 82], [82 / 82], [82 / 82],
[bmss], [42 / 42], [42 / 42], [42 / 42],
[chain], [6 / 6], [6 / 6], [6 / 6],
[char2], [13 / 13], [13 / 13], [13 / 13],
[csidh], [9 / 9], [9 / 9], [9 / 9],
[endo], [3 / 3], [3 / 3], [3 / 3],
[find], [12 / 12], [12 / 12], [12 / 12],
[genus2], [6 / 6], [6 / 6], [6 / 6],
[kernel], [28 / 28], [28 / 28], [28 / 28],
[kernel2], [9 / 9], [9 / 9], [9 / 9],
[kernel2_even], [4 / 4], [4 / 4], [4 / 4],
[kernel2_mont], [6 / 6], [6 / 6], [6 / 6],
[models], [18 / 18], [18 / 18], [18 / 18],
[newcomp], [14 / 14], [14 / 14], [14 / 14],
[p1kernel], [36 / 36], [36 / 36], [36 / 36],
[path], [15 / 16], [15 / 16], [15 / 16],
[phi], [46 / 46], [46 / 46], [46 / 46],
[quat], [3 / 3], [3 / 3], [3 / 3],
[radical], [9 / 9], [9 / 9], [9 / 9],
[relation], [2 / 2], [2 / 2], [2 / 2],
[theta], [12 / 12], [12 / 12], [12 / 12],
[theta4], [5 / 5], [5 / 5], [5 / 5],
[twopow], [6 / 6], [6 / 6], [6 / 6],
[walks], [8 / 8], [8 / 8], [8 / 8],
)

The 24 explicit command groups emit these 24 reported labels; aliases and split subgroups prevent a one-to-one command mapping. All rows, raw timings and verification flags remain preserved. The one false flag on every profile is ghs_with_kohel_volcano on a 17-bit fixture; its note describes split primes that may not generate the needed class-group orbit. A false flag is retained, not converted into a successful path.

AMD and ARM use original source 79006eb85a40 from hosted run 38766895864. Their outer native isolation receipts pass, but those artifacts did not retain the legacy bench executable. The new M4 run uses fourth-candidate source 9a55ce6d22b9473c1d2afbf398bfdf51ae22b544 and retains executable, source, compiler and shell/resource receipts. These different revisions are not pooled or compared for speed.

Each record is a legacy stage measurement, not an exact-output attestation. Neither complete input/output fingerprints nor matched A/A and interleaved A/B controls exist for this harness. The shared M4 supplies no qualified wall timing. Rows marked verified refer only to their existing harness-specific check, not independent replay. Named P192/P224 whole-search cost and the remaining exact family adapters stay open.

Census SHA256 receipts:
#text(size: 5.5pt)[ubuntu-24.04: 3b7a7339c8f1c6d175ec67ccee6892c4f9a60ee62eb6d5c5a18bd85757e8b85b\
ubuntu-24.04-arm: 3cda02896896bbceda877a3afe01f7c6d8400b5c15ffc281d19fda305ca603b2\
apple-m4-pro: 3381aa7c27b8d124170eb5d9b7dd80080ebfd04f4e82f0475358f96447b7f4a4\
]
