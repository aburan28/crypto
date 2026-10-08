# Fixed-quadratic affine batches: does Boolean F5 keep its selected row skeleton?

## Mathematical question and reference

The [merged graded-block study](https://github.com/aburan28/crypto/pull/1207)
verified an exact cubic-block reuse mechanism on generated quadratic Boolean
matrices, but its strongest-control discovery passed only 2/8 groups and
stopped before holdouts. It did not run the repository's inherited Boolean
matrix-F5 call. At degree four, the quartic projection of `t*(q_j+a_j)` is
again independent of the affine tail `a_j` when each generator has fixed
quadratic part `q_j` and multiplier `t` has degree at most two. That identity
alone does **not** make the F5 matrix reusable: the Boolean-specific criterion
may prune different multiplier rows as the affine coefficients change.

This experiment measures the exact **selected row-label signature** before
attempting an F5 batch cache. Its reference is the current native
`F5Criterion::new` implementation in
`src/cryptanalysis/matrix_f5_f2.rs`, with no criterion change. It uses the
repository's rewrite-only Boolean-safe criterion and does not invoke an
ordinary comparable-signature rule. No solver timing, relation yield,
curve equation, scalar or key input is included. A signature match is a
necessary reuse condition, not a speedup.

## Inputs and exact signature

Use n=12/16/20/24, m=n, degree bound D=4 and full n-variable multiplier
mask. For each n and seed, initialize the same SplitMix64 recurrence and
2n-distinct-quadratic-term generator used by
`research/boolean_graded_tail_reuse_20261003/PROTOCOL.md`. Continue the PRNG
stream to form affine tails. Each batch starts with all affine coefficients
zero. In `independent_affine`, every later assignment draws one new constant
and n linear coefficient bits per generator from successive low bits. In
`walk_affine`, every later assignment toggles one affine slot per generator;
slot zero is the constant. Quadratic coefficients, generator order and
degrees remain fixed. All inputs are generated public Boolean systems.

For a generator `j`, enumerate all squarefree multipliers of degree at most
two in ascending numeric order. Record the Boolean selection bit
`!F5Criterion::prunes(j,t)` for every `(j,t)` label. Concatenate in generator
then multiplier order, and SHA-256 hash the complete bitstream. Also retain
the exact bitstream, selected/pruned counts, Koszul/Frobenius prune counts,
lower-level row/zero counts and criterion GF(2) word-XOR count. The native
verifier must regenerate every polynomial and selection bit from source,
check hashes and counts, and reject a mutated result. A zero row is not
silently counted as a selected nonzero basis row; this is a **row-label**
screen, not a matrix rank or relation count.

For each batch, report whether every changing assignment has exactly the
base assignment's signature, the number of distinct signatures, the largest
signature class fraction, and the fraction of adjacent pairs with identical
signatures. Different signatures prohibit replay of one fixed selected-row
matrix without fallback. Equal signatures do not prove that the high-block
elimination is fast or that full F5 output can remain packed.

## Frozen grid and stop rule

`protocol.json` fixes n=12/16/20/24, batches 2/8/32, both families,
discovery seeds 20261005/3141601/2718293 and unused holdout seeds
20261012/4242433/5772169. This gives 72 cells per phase and 1,008 system
evaluations per phase. Apply a 600-second worker cap, 8,192 selected-row
budget per system for the screen, and a 64 MiB output-evidence cap per run.
Any cap or incomplete criterion returns CENSORED, never a pass. The frozen
source, protocol, binary hash, host, raw bitstreams and native replay receipt
are sealed before interpreting results.

Advance to holdouts only if **every** n>=16, batch-32, family/seed discovery
cell has a largest exact signature class covering at least 80% of its 32
assignments and every oracle/record check passes. This is a conservative
structural feasibility gate, not a timing gate. If it fails, retain all cells
and stop without consuming holdouts. If it passes, run the untouched holdout
seeds with identical source and require the same all-cell condition. A later
separate protocol would charge F5 criterion recomputation, high-block setup,
per-assignment low work, output unpacking, memory and fallbacks against fresh
same-binary complete F5 calls. It would need a full-method relation/rank and
rho comparison before any cryptanalytic claim.

Use Rust for generator, criterion calls, verifier, analysis and evidence
generation. Thin shell build/run commands may call the repository's required
CPU-isolation controller. Historical Python files are context only and must
not become an execution path. No timing or result is asserted by this plan.
