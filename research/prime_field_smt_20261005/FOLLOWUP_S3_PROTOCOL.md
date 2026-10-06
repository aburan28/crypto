# Follow-up: native finite-field Semaev S3 encoding

Status: preregistered after both full-coordinate frozen encodings returned
`unknown` at 120,000 ms and before any S3-backend solver run.

The frozen cryptanalysis input, target, bound `B=4096`, cvc5 1.4.1
CoCoA-enabled artifact, `--ff-solver=split`, 120,000 ms cap, and all claim
exclusions remain unchanged.

## Fixed formulation

The exporter enumerates every `x` in `[0,B)` and retains exactly those for
which `x^3+a*x+b` is a quadratic residue in `F_p`. A deterministic Boolean
trie over the little-endian field bits constrains each of `x1,x2` to this
liftable set. The query then asserts only

```text
S3(x1,x2,xR) = 0
```

in cvc5's native prime-field theory, plus exact bitwise `x1 <= x2`
symmetry breaking. The emitted receipt records the allowed-abscissa count
and digest. No `y` or slope variable is sent to SMT.

For every SAT pair, the native checker computes both square roots for each
abscissa and tries all four sign combinations. It accepts only a pair whose
independently recomputed affine sum is the original signed target; the slope
and full witness are then replayed through the existing verifier. Because
the membership trie admits only `F_p` lifts and both signs are available,
the defining property of S3 makes this lifting complete.

## Runs and stop rule

Run the planted SAT and exhaustive UNSAT toy controls, then the frozen
20-bit case once. A timeout, `unknown`, malformed model, unliftable model,
or failed verification is retained. The trie, variable order, S3 formula,
solver option, target, bound, and cap are not changed after the frozen run.

This is still a two-summand decomposition-stage diagnostic. The direct
`R-P` factor-set reference is expected to dominate it. Pollard-rho ratio,
end-to-end `S`, and speedup remain `null`; the scoreboard is unchanged.

