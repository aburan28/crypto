# Prime-field two-point SMT protocol

Status: preregistered before the external-solver runs.

This round implements an exact SMT-LIB `QF_BV` frontend for finite
two-point decomposition on a short-Weierstrass curve over a prime field.
It is a solver-stage experiment, not a discrete-logarithm experiment and
not a Weil descent: a literal prime field has no extension coordinates to
descend to Boolean equations.

## Frozen question and input

For a finite target `R` on

```text
E/F_p: y^2 = x^3 + a*x + b,
```

find finite `P1=(x1,y1)` and `P2=(x2,y2)` such that
`P1 + P2 = R`, `0 <= x1,x2 < B`, and `x1 <= x2`.

The cross-repository case is translated without changing field elements
from `aburan28/cryptanalysis`:

- last fixture commit: `a14ca487a5659fc6e6572133f7c89a3326b3e1b5`;
- path: `challenges/ecc/curves/fp-random-b20.json`;
- Git blob: `40b380492c7a250babaa2b5ae86a4cb5592fe598`;
- translated input: `instances/cryptanalysis_fp_random_b20.json`;
- factor-base bound: `B = 4096` (frozen before solving);
- external solver cap: 120,000 ms.

The external solver is cvc5 1.4.1, official Linux x86-64 static archive
`cvc5-Linux-x86_64-static.zip`, SHA-256
`2f8efe58fe27ba7bccbb504533f690b9312d69da14192712460e4a19231f02a1`.
The executable hash is recorded separately after extraction.

## Exact encoding contract

The SMT query declares `x1,y1,x2,y2,lambda` as `w=ceil(log2(p))` bit
vectors and constrains all five to be strictly below `p`. Field addition
and subtraction widen to `w+1`; multiplication widens to `2w`; reduction
by `p` happens before extracting the low `w` bits. Thus host or bit-vector
overflow cannot change a field operation.

Both complete finite-sum branches for odd characteristic are encoded:

1. `x1 < x2`, with
   `lambda*(x2-x1)=y2-y1` and the ordinary affine addition equations;
2. `x1=x2`, `y1=y2!=0`, with
   `(2*y1)*lambda=3*x1^2+a` and the doubling equations.

The inverse-point branch sums to infinity and therefore cannot equal the
finite target. Both input points also satisfy the original curve equation.
A `sat` model is accepted only when native big-integer arithmetic checks
canonical ranges, both curve equations, the branch slope, and an
independent affine group addition back to the target.

## Controls, hypotheses, and stop rules

The correctness controls are fixed as follows:

- unit tests for decimal, hexadecimal, binary, and SMT `(_ bvN w)` model
  values;
- exhaustive small-prime agreement, including planted SAT and complete
  UNSAT cases;
- rejection of a corrupted witness and a malformed/non-prime instance;
- one planted cvc5 integration case before the frozen cross-repository
  case;
- the frozen cross-repository case, once, under the 120-second cap.

The functional hypothesis is that every reported SAT model verifies and
that the frozen 20-bit case completes within the cap. A timeout, `unknown`,
malformed model, nonzero solver exit, or failed independent verification
falsifies that completion hypothesis and is retained as the result. The
bound, target, encoding, and cap are not changed after that outcome.

## Reference and accounting

The native reference enumerates all curve points with `x < B`, stores the
exact signed points, and checks `R-P` membership for every factor point.
It reports factor points and candidate probes. This establishes SAT/UNSAT
for the two-summand relation without trusting the SMT model parser.

The SMT receipt records input/query hashes and byte counts, field and
product widths, solver binary hashes before and after execution, exact
arguments, timeout, exit status, stdout/stderr hashes, elapsed wall time,
parsed verdict, witness, and verification receipt. Wall time is a
practicality diagnostic only; no conversion to group additions is made.

The only comparison table for this round uses solver outcome and verified
correctness. Ratios to Pollard rho, end-to-end `S`, and speedup are `null`:
neither relation collection, linear algebra, nor scalar recovery is in
scope. Consequently this stage-only round does not change the IC
scoreboard or leaderboard.

