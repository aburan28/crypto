#set page(paper: "a4", margin: 16mm)
#set text(font: "Arial", size: 8.5pt)
#set heading(numbering: none)
= Seven kernel constructor baselines
2026-10-15 · Two repeated native 84-case runs

Seven methods at 192/224/256/511 field bits, degrees 7/31/101. All 84 paired output fingerprints and all four field-call counters repeat exactly. Synthetic supersingular fixtures; these are not the named NIST curve models.

== Measured calls at 224 bits / degree 101
#table(columns: (1.6fr, 1.3fr, 1fr, 0.8fr), inset: 4pt, [Method], [Supplied input], [Mul calls], [Inv calls],
[velu_from_point], [affine generator], [454], [50],
[velu_cyclic], [affine generator], [454], [50],
[velu_xonly], [generator x], [863], [49],
[velu_xonly_fast], [generator x], [1044], [1],
[sqrt_velu], [affine generator], [957], [15],
[sqrt_velu_fast], [generator x], [2532], [2],
[kohel_given_kernel], [kernel polynomial], [7999], [4],
)
#image("COUNTS.svg", width: 100%)
#image("PIPELINE.svg", width: 100%)

Generic Field API calls, inversion internals opaque; n=2 identical recorded counts. Different supplied-input boundaries remain separate. Reference/fixture generation, fingerprints and verification are excluded; owned representation allocation is included. These baseline counts establish no kernel speedup.

Degree, codomain and full point transports at 2ell+1 distinct non-kernel abscissae agree with both Kohel and point-enumerating Vélu. Counts are saved before verification. These fixture fingerprints are distinct from the BMSS coefficient comparison.

Shared Apple M4 Pro, 14 logical CPUs, 48 GiB, macOS ARM64. Zero timing samples requested or qualified. Other physical CPU baselines and optimization rounds remain open.

Source: #text(size: 7pt)[5521605f86ae444a9f4a433ca9b2f80d541e8470]
Binary SHA256: #text(size: 6.5pt)[e22a305acbd0ab58be40a7fe387867f10598255b3cf5ec03532d98fbb1b50c6e]
Comparison SHA256: #text(size: 6.5pt)[80dc2ebcf465fade24a71e55c47f366feed3ab200f1136196572b7ee8034c496]

Read REPORT.md, COUNTS.csv and KERNEL_BASELINE_PROTOCOL.md with the frozen raw records and source/binary archive. Both complete matrices, including all other widths and degrees, remain in the raw comparison.


Kohel's result contains all rational-map coefficients. The other constructor results store compact evaluation representations. Their build counts include different output materialization costs.
