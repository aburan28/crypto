# Second candidate: multiplier inlining

The first candidate's P224/1471 verification profile records 49,944,361,778 instructions.
Of those, 39,279,219,914 (78.65%) are attributed to independent field multiplication and
9,533,786,830 (19.09%) to the polynomial product routine. The first candidate changes
verification instructions only slightly relative to the frozen baseline; retain that
result rather than treating the screening improvement as a verifier improvement.

Test a second candidate that forces inlining of the independently implemented field
multiplier and its tiny operations. Keep its arithmetic, input/output contract, independent
implementation, catalogue and source registry unchanged. Cross-check field arithmetic
and all exact-map tests before measurements. Compile with the same Rust version and
release flags as the first candidate; no global CPU feature flags.

Compare this candidate on the same five-round instruction panels, with particular focus
on the full P224 verification command. If forced inlining regresses instruction counts,
binary size or correctness, retain that evidence and do not promote it as a gain. These
measurements remain ARM Linux instruction counts, with native Mac throughput unset.
