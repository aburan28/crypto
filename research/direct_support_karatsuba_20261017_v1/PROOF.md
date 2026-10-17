# Adaptive direct source-product arithmetic

Keep the exact source-product/coefficient equality contract of the prior direct
compiler. Norm generation, solving and original root outputs are omitted from
this representation subsystem, and no component count is an8x solver result.

All polynomial coefficients here are Boolean wires over F2. For balanced
products split x=x0+T^m*x1, y=y0+T^m*y1 and compute
z0=x0*y0, z2=x1*y1, z1=(x0+x1)(y0+y1)+z0+z2.
The exact result is z0+T^m*z1+T^(2m)*z2. Addition is XOR, so no sign or
carry assumption is imported. If one input is more than twice the other's
length, divide the longer into short-length blocks and sum shifted block
products. Induction on maximum length proves exactness. Every recursive
balanced child has length at mostceil(max/2), and every unbalanced block is
strictly shorter than its parent; zero and tiny inputs use schoolbook.

The monic-quintic coefficient circuit, unique global lift q, root ordering,
multiplicities and exact equality constraints are unchanged. Only arithmetic
implementation differs. Each arm has a separate identical gate arena; no
cross-arm caching removes charged work. AND/XOR counts include Boolean
constant/idempotence simplification and identical hash-consing. Standard
Tseitin clauses retain the same3/4/unit model, not a SAT-runtime prediction.

Freeze cutoffs0(reference),8,16,32 on six binary field/source widths. The
prediction is that at least one FIXED preregistered nonzero cutoff reduces
AND gates by at least25percent at every width, with standard clauses no
more than1.10times reference at any width. Preserve every arm and regression;
selecting an exposed width-specific winner is calibration, not holdout success.

Independent controls cover27648 literal coefficient assignments/mutations,
plus288 separate convolution vectors at18x18,35x18,52x18,69x18,140x140,
and557x140 for all three cutoffs. The wide vector oracle is plain Boolean
convolution, not recursive multiplication. Tiny coefficient cases alone do
not exercise recursion; the extra vectors are required and were run.
90seconds/4GiB/oneCPU. Serialized profile completion and supervisor lifetime
remain charged. No solver, original root count, external-fastest, novel
priority or qualified end-to-end comparison is asserted.
