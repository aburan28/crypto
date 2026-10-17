#set page(paper: "us-letter", margin: (x: 0.8in, y: 0.7in))
#set text(font: "New Computer Modern", size: 9.5pt)
#set par(leading: 0.55em)

= P-192 Type I conditioned three-sum support

The native finite census found *224 accepted targets among 816 distinct fixed-base three-sum targets*. The four character cells are ++: 224, +-: 211, -+: 203, and --: 178. A known scalar a gives P = aG and the quotient query q = -2P = a(-2G). The two-character filter changes eligible fixed-base support from 816 to 224 targets. A native replay checked the 64 frozen known-scalar equations and found zero matches; the separate relation verifier records zero useful rank.

#image("support.svg", width: 100%)

The instance is standard NIST P-192 over p = 6277101735386680763835789423207666416083908700390324961279 with prime subgroup order r = 6277101735386680763835789423176059013767194773182842284081. The 16 quotient points are frozen from the earlier verified base (fixture SHA-256 3984b55b793305b30beaf18918209282bce80d444c96217cd6143c7d9470d409). The character polynomial and accepted-set lower bound are frozen in the density receipt (SHA-256 fd0a47cf4aecf736b1852d133fd384e30d9584114761bac6495f803d4210ca29).

The Rust program verifies curve membership and r-torsion for every base point, enumerates all unordered triples, checks all target points are distinct, evaluates the sextic character and its cubic factors, and checks q = -2P. It replays 64 frozen accepted scalars, checks q = a(-2G), and finds zero hits. An earlier independent direct fiber-character replay obtained the same four cells. The source, JSON receipt, raw log, and replay records sit in the parent study directory.

With the certified accepted-scalar lower bound L = 1569275433846670190958947353987675545133315652035359680560, uniform conditioned-query hit probability is at most 224/L < 2^-182. Lean checks the exact integer inequality; finite-field enumeration is native and independently replayed. Degree-three divisor pairs are verified for all 16 base classes and all 64 accepted query classes. A unique fixed-origin Jacobian normal form and a larger relation collector remain separate tasks.

Run: 2026-10-15 08:17 UTC; Darwin arm64; Rust 1.93.1; git baseline ecb238fd9eef5ed6f699a609c803e41dc74e4ef9. Native run and two scoped tests exited 0. The repository-wide release library test stopped during compilation because the volume ran out of space. This is a finite support census; CPU timing and end-to-end attack cost were not measured.
