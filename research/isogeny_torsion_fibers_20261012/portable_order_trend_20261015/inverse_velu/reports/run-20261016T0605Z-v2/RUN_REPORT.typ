#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Direct inverse-Vélu fiber oracle: frozen V2 run
*Public `p=1237, r=1201` subgroup · 2026-10-16 · native Rust.*

The inverse x-coordinate cubic selects source witnesses for a
different-`j` degree-three Vélu destination. The oracle scanned 29
destination points, made 29 polynomial calls and 22 split attempts,
checked 28 complete maps, and selected the same eight factor-base
columns as a 1,531,275-point source census.

#image("PIPELINE.svg", width: 100%)

Five frozen arms ran at 1/4/16/64 public targets over three paired
rounds. All 60 children exited zero. An independent Rust checker
replayed 1,275 scalar outputs and 6,126 relation witnesses. Across
all 255 paired inverse/full target positions, support keys, trials,
rank, pair entries and recovered scalars matched. The 368-file raw
manifest passed SHA-256 verification. The largest inverse-child RSS
at 64 targets was 9,125,888 bytes under the 512 MiB cap.

At 64 targets, counted group-cost lower bounds per target were
`S_group=151.595` inverse, `1532.302` full census, `4.071` signed rho
and `1.315` BSGS. Paired full/inverse was 10.108, while inverse/rho
was 37.235. Complete portable `S` remains unset because native work
outside the counted group axis is not fully priced.

#image("GRAPH.svg", width: 100%)

The graph reports paired whole-child process CPU: full/inverse
43.664, inverse/rho 9.159, inverse/BSGS 11.380, and rho A/A 0.863
at 64 targets. The host was shared L0; its CPU-model `sysctl` was
denied, and the host receipt preserves that limit. No L2 wall result
is inferred. `REPORT.md` links the frozen spec, raw outputs, replay,
analysis, source, identity proof and next-test gates.
