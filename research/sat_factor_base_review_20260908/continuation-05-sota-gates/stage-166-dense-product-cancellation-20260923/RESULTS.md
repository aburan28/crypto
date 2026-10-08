# Stage 166: dense cancellation before F4 product sorting

Strong single-target internal engineering on one finite public toy PDP. This is not a full solver panel, full index-calculus/rho crossover, independent reproduction, novelty review, or SOTA.

For Boolean systems with at most 20 variables, monomial multiplication now toggles mapped output masks in one generation-tagged dense array. Each first-seen mask stores a 32-bit key encoding descending degree and ascending mask order. The solver sorts only odd-parity survivors; `PQ_F4_DISABLE_DENSE_MUL=1` restores the prior map-all, sort-all, then-cancel path.

The clean blind F4 process takes 19.887753 wall seconds, 19.858236 core-seconds, and 900169728 bytes peak RSS. Across 818,171 products it maps 651,502,431 input terms to 549,909,017 survivors and cancels 101,593,414 terms before sorting.

The exact clean disabled control takes 21.582929 wall seconds. The selected clean speedup is 1.085x; three interleaved development pairs give 1.071x. Equation fingerprint, all algebraic and matrix work counters outside the new product counters, roots, and the exact curve-verified witness agree.

Direct MITM remains 6.63x faster by wall, 6.64x cheaper by CPU, and 19.64x smaller by RSS.

Licensed Magma F4, the complete same-instance solver panel, full index-calculus/rho crossover, and unaffiliated reproduction and novelty review remain open. This is not a Koblitz index-calculus SOTA result.
