# Prime-field same-subgroup rho control, frozen before V6 runs

The V5 compact run on the order-67 cell put the pairing cold CPU below the
generic extension-curve negation rho CPU at 64 distinct targets in all three
shared-host rounds: paired pairing/rho ratios `0.8893, 0.5630, 0.9603`.
Both methods independently replayed every scalar. The source generator and
all targets have coordinates in `F_401`, however. A rho implementation can
compute directly over that smaller field. V5 is therefore a representation
control, not a qualified method crossover.

V6 keeps the public targets, selection, pairing, field DLP, extension rho,
seeds, per-attempt step cap, 32 charged retries and prefixes fixed. Add a
same-target negation rho arm using `PrimeCurve(p=401,a=0,b=1)`, with direct
word-field points and the repository's native prime-field walk adder. Check
each point on the prime-field curve, charge group validation and all failed
walks, and retain the prime field multiplication/squaring/inversion ledger.
Run three rounds, rotating arm order across pairing/extension-rho/prime-rho.
The same public-input parsing and validation is a common separately measured
charge; report ratios with and without that charge so it cannot hide a
representation advantage. Emit only after all timers finish. Separately
replay all three scalar arrays with independent field and point arithmetic.

Compare complete cold CPU at `N=1,4,16,64,66`, including pairing setup and
rho validation and restarts. This remains shared-host L0 timing, not an
isolated wall-time or counted-operation result. If the prime-field rho
reference is faster, keep the V5 outcome as the specific weaker-control
observation rather than promoting it to a general ECDLP speed claim.
