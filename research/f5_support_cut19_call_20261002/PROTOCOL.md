# Cut-19 exact matrix-F5 complete-call candidate

## Source and hypothesis

Start from the exact cut-20 full-rank implementation at `9a5c835af`.
The separate cut screen (`codex/f5-support-cut-screen-20261002`, commits
`2f43c977e` and `57740ed4d`) showed cut 19 certifies all four frozen
matrices and uses 83.1–83.3% of cut 20's counted reduction XORs. Fix
the production support cut to variables 0–18, without the exploratory
environment selector. Keep the exact fallback and all other output
forms unchanged. The hypothesis is that this lower rank work improves
the complete n24 degree-4 call enough to pass the unchanged further-2×
gate against selective echelon.

## Frozen local screen

Use the seven F5 cases and seed XORs `0`, `badc0de1`, `5eed2026`,
`f5c02a28`, one release binary, one Rayon thread, and the same direct
packed-row, GF(2), and unpack settings as the cut-20 complete-call
protocol. Compare form 2 selective echelon with the new form 5 cut 19
in separate processes. Per seed use one warmup per arm, five
reference/reference A/A pairs, and five alternating reference/candidate
pairs. The native Rust runner must preserve every status, full output,
source/binary hash, exact rank, canonical row-space fingerprint, route,
counted work, exclusive phases, and complete-call timing.

The primary may return different raw rows from selective echelon, but
must match exact rank, canonical row space, F5 criterion, built/pruned
row counts, and columns. All six smaller cases must match every
nontiming output and route field. The cut-19 certificate must succeed
on all four seeds and use at most 35% of selective-echelon counted
reduction XORs on each. Apple ARM64 wall measurements are nonpromoting;
report paired medians, exact five-pair bootstrap intervals, and A/A
ranges without timing-based seed selection. Advance to Linux only if
exactness, route, and counted-work conditions all pass.

## Final claim gate

On physical isolated Linux x86-64, recompile the exact branch and run
the same four-seed one-thread paired complete-call protocol. Each seed
must have median **and** exact five-pair bootstrap 95% lower bound both
above 2.00× against selective echelon. Every smaller-case median must
avoid regression beyond its A/A lower bound, and a separate two-thread
job must meet the same nonregression control. Preserve the first clean
attempt per seed without selecting by timing, plus all failed,
contended, timed-out, or unavailable attempts. This is a solver-stage
claim only; it does not establish an IC online or rho speedup.
