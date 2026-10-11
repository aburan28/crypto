#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm), numbering: "1")
#set text(size: 9.5pt)
#set table(stroke: 0.4pt, inset: 4pt)
#show heading: set block(above: 1.2em, below: 0.6em)

#align(center)[
  #text(size: 15pt, weight: "bold")[Run report: run-20261011T1813Z-bench] \
  Eigen-Séta v2 Rust KEM: benchmark run (descriptive record) \
  #text(size: 8.5pt)[Study: research/eigen_seta_rust_20261011 · sources: commit 19c641d121097e80f727bcc650a08632e204ba1d]
]

= What ran

#table(columns: (auto, 1fr),
  [Start (UTC)], [2026-10-11T18:13:51Z],
  [Sources], [`eigen_seta/` at commit 19c641d12. The binary was built from that tree just before the commit, and a cargo rebuild after the commit was a no-op.],
  [Toolchain], [rustc 1.97.0 (2d8144b78 2026-07-07). Release: opt-level 3, thin LTO, 1 codegen unit, default target-cpu.],
  [Host], [Intel Xeon Processor \@ 2.80 GHz, 4 logical cores, 15 GiB, Linux 6.18.44 x86_64 (adx, avx2, avx512f, bmi2, popcnt)],
  [Isolation], [Shared cloud VM. No pinning and no frequency control. One benchmark process at a time.],
  [Command], [`bench_field`, then `run Toy16 5`, `run Toy32 3`, `run L1 3`, `run L3 2`, `run L5 2`],
  [Workload], [One keygen per set, then the stated number of encaps and decaps trials. Keys and ciphertexts come from OsRng.],
  [Output pin], [Each trial checks that decaps returns the encapsulated key. Result: 15/15 correct.],
  [Records], [`records/bench_output.txt`, `records/host_and_command.txt`],
)
Requested: 5 sets, 15 trials, and the field micro-benchmark. Completed: all. Failures: none.

= Measured values (wall clock, one thread)

#table(columns: 5,
  [*Set*], [*keygen (s), n=1*], [*encaps (s), per trial*], [*decaps (s), per trial*], [*correct*],
  [Toy16], [0.381], [0.0055, 0.0070, 0.0070, 0.0084, 0.0052], [0.131, 0.159, 0.186, 0.126, 0.118], [5/5],
  [Toy32], [2.341], [0.0095, 0.0116, 0.0091], [0.510, 0.446, 0.591], [3/3],
  [L1], [32.279], [0.2448, 0.2588, 0.2394], [14.223, 14.997, 13.424], [3/3],
  [L3], [118.589], [1.0256, 1.0304], [51.834, 53.940], [2/2],
  [L5], [271.791], [2.2045, 2.1651], [123.953, 120.180], [2/2],
)
Public key = ciphertext sizes: 216, 376, 1344, 1992 and 2632 bytes, from Toy16 to L5. Secret-degree bit lengths in this run: 51, 80, 278, 409 and 536.

#figure(image("kem_times.svg", width: 88%), caption: [Wall time per operation, every trial shown, log scale. Measured; records/bench_output.txt.])

#table(columns: 5,
  [*Set*], [*limbs*], [*mul (ns), min of 15*], [*inv. Fermat (ns), mean of 200*], [*inv. safegcd (ns), mean of 200*],
  [Toy16], [4], [54], [18,403], [4,741],
  [Toy32], [6], [72], [67,657], [10,316],
  [L1], [21], [611], [1,426,022], [35,608],
  [L3], [32], [1,516], [6,708,390], [109,506],
  [L5], [42], [2,449], [13,453,052], [96,327],
)
No spread is recorded beyond the table. The L3 safegcd figure is above the L5 one. That is what was recorded, and this run does not explain it.

#figure(image("field_ops.svg", width: 88%), caption: [F_p multiplication and inversion against limb count, log scale. Measured; records/bench_output.txt.])
