# Width 6 versus width 8: paired complete-call follow-up

Width 6 passed the predeclared ≤90% counted-work gate on all four seeds
in `RESULT.md`. Before paired timing, freeze the same four systems and
seven F5 cases, one release binary, one Rayon thread, and the exact
direct packed-row, support cut-19, and unpack settings. The reference
sets the rank-table maximum to 8; the candidate sets it to 6. Both
arms request F5 form 5. Use separate processes, one warmup per arm,
five width-8/width-8 A/A pairs, and five alternating width-8/width-6
pairs per seed. Preserve every call, including failures and timeouts,
with source/binary hashes, host, rank, canonical and raw fingerprints,
route, counted word XORs, phase and complete-call times.

The primary must match every nontiming output field except reduction
word XORs; the original rows, term counts and raw fingerprint must
remain identical. All smaller cases must match every nontiming output
field. Width 6 must keep ≤90% counted work on every seed. Report the
paired complete-call median, exact five-pair bootstrap interval and
A/A range even if wall time regresses. Apple ARM64 timing is
nonpromoting. Advance width 6 to the physical Linux selective-echelon
gate only if all output checks and counted-work checks pass and no
primary paired median falls below that seed's A/A minimum. A smaller
case below its A/A minimum remains an explicit regression signal.

The unchanged requested final gate is four-seed physical isolated
Linux x86-64 one-thread complete-call median and exact bootstrap lower
bound both above 2.00× against selective echelon, plus smaller-case
and two-thread nonregression controls. This is a solver-stage claim,
not an IC online or rho result.
