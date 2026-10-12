# Frobenius quotient follow-up: independent SAT reproduction

This is an independent re-implementation of the target-only four-point
decomposition formulations recorded in
[`../frobenius_quotient_followup_20260926`](../frobenius_quotient_followup_20260926/README.md).
That note's solver code is not in this repository, so its solver table was
recorded as reported, not reproduced. This directory re-runs the grid with new
code, fresh targets and exhaustive checking of every answer.

It is a **stage diagnostic** on toy fields: `S`, end-to-end cost and speedup
are unset, and no work below `2^61` is claimed. Each cell is a single timing
trial on 4 targets, so no speedup factor is established.

## Files

| file | role |
|---|---|
| `encode.py` | Circuit builder. Each field product of linear forms becomes n² AND gates plus XOR constraints; squaring, sqrt, trace and constant multiplication stay linear. Native XOR, or XOR expanded into clauses (`_cnf`). |
| `run.py` | One fresh single-thread CryptoMiniSat 5.16.0 solver per attempt. Every proposed tuple is checked on the curve; a rejected tuple is blocked and the search continues within the same budget. Resumable with `--resume`. |
| `brute_check.py` | Decides by exhaustive pair matching whether each target is decomposable, and flags any UNSAT or verified answer that disagrees. |
| `test_encode.py` | Checks the S3 formula, the XOR-to-clause expansion, the field-product circuit, and that planted witnesses (with and without phases) satisfy every formulation. |
| `summarize.py` | Builds the tables below from `results/*.jsonl`. |

```sh
pip install pycryptosat==5.16.0
python3 -m unittest test_encode -v
python3 run.py --n 13 --s 4 --budget 300 --phases --out /tmp/n13.jsonl
python3 summarize.py /tmp/n13.jsonl && python3 brute_check.py /tmp/n13.jsonl
```

## What the model is given

The model receives the curve `y² + xy = x³ + 1` over `F_{2^n}`, a
factor-space basis, and the target `R = (a, b)`. It never receives a list of
valid payloads, orbit keys or a decomposition.

The factor space is a seeded random `s`-dimensional subspace of `ker Tr`,
redrawn until it holds at least 4 valid x-coordinates. Its basis is recorded
in each result file. Unknowns are:
- the factor-space coordinates `c_i` of each summand's representative `v_i`;
- with `--phases`, a one-hot Frobenius phase `e_i`, so that `x_i = v_i^(2^k_i)`;
- `t`, the x-coordinate of `±P1 ± P2`;
- `u` (`specialized`) or the slope `ℓ` (`line`).

Constraints:
- `S3(x1, x2, t) = 0` and `S3(x3, x4, u) = 0`.
- Final equation, `specialized`: `S3(t, u, a) = 0`.
- Final equation, `line`: `u = t + ℓ² + ℓ + a` and `t·u = aℓ² + a² + b + a`.
- `subgroup` variants add the quadratic certificate
  `q² + √v·q + 1 = 0` with `Tr q = 0` on every representative.
- All four x-coordinates are distinct.

Two target sets are used per field, each drawn fresh under seed `20260930`:
- **planted:** a signed sum of 4 distinct admissible points, with a random
  phase per summand when phases are unknown;
- **random:** uniform in the odd subgroup `H = 4E`, excluding the identity.

## Results

Each cell reads *verified / attempts*, followed by the total wall seconds
across attempts (build, load, solve and verification included). The unsat
column counts attempts where the solver proved no solution exists. Every
unsat and verified answer below agrees with `brute_check.py`, with 0
mismatches in every file.

### Representatives only (no phase): the model is much smaller than the note's

At n = 13 the factor space holds 5 valid representatives. Only 4 of the 8
targets used are decomposable at all: exactly the 4 planted ones.

| n, budget | targets | spec_xor | spec_cnf | subgroup_xor | line_xor | line_subgroup | line_subgroup_cnf |
|---|---|---|---|---|---|---|---|
| 13, 2 s | planted | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 1/4 |
| 13, 2 s | random | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 |
| 13, 300 s | planted | 4/4; 37 s | 4/4; 74 s | 4/4; 16 s | 4/4; 56 s | 4/4; 60 s | 4/4; 26 s |
| 13, 300 s | random (all unsat) | 0/4; 465 s | 0/4; 1009 s | 0/4; 91 s | 0/4; 513 s | 0/4; 118 s | 0/4; 213 s |
| 19, 2 s | planted and random | 0/8 in every column |  |  |  |  |  |
| 19, 300 s | planted, 2 of 4 targets | 0/2 | 0/2 | 0/2 | 0/2 | 0/2 | 0/2 |

The n = 19, 300 s grid was stopped by the harness's background-job time limit
after 12 of 48 attempts. Every completed attempt timed out, and the run can
be resumed.

### Frobenius phase unknown (the note's setting)

At n = 13 the 5 representatives give 65 admissible x-coordinates, and all 8
targets are decomposable. This matches the note's saturation claim for
degree 13.

| n, budget | targets | spec_xor | spec_cnf | subgroup_xor | line_xor | line_subgroup | line_subgroup_cnf |
|---|---|---|---|---|---|---|---|
| 13, 2 s | planted | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 |
| 13, 2 s | random | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 |
| 13, 300 s | planted | 4/4; 176 s | 4/4; 674 s | 4/4; 164 s | 4/4; 299 s | 4/4; 160 s | 4/4; 226 s |
| 13, 300 s | random | 4/4; 285 s | 3/4; 611 s | 4/4; 131 s | 4/4; 175 s | 4/4; 138 s | 4/4; 529 s |

At n = 19 the factor space holds 22 valid representatives, giving 418
admissible x-coordinates. The two fastest n = 13 models were then given a
longer budget on 4 planted targets. All 4 targets are decomposable (checked
exhaustively).

| n, budget | targets | subgroup_xor | line_subgroup |
|---|---|---|---|
| 19, 1200 s | planted | 0/4; 4921 s | 0/4; 4910 s |

Each n = 19 model has about 7.5k variables and 5.4k AND gates, against 3.6k
and 2.5k at n = 13.

Tuples proposed and rejected by the curve check (300 s, n = 13, phases, summed
over all 8 targets):

| spec_xor | spec_cnf | subgroup_xor | line_xor | line_subgroup | line_subgroup_cnf |
|---|---|---|---|---|---|
| 20 | 12 | 0 | 23 | 0 | 0 |

## Reading

- **2 s says little.** With unknown phases, 0 of 48 attempts at n = 13 finish
  in 2 s, yet 47 of 48 finish in 300 s on the same targets. Solve times range
  from 12 s to 303 s, with a median of about 50 s. The note's sporadic
  4-of-56 at 2 s is consistent with this. The earlier conclusion that
  "success remains sporadic" describes the budget, not the formulations.
- **The subgroup certificate helps and never hurts here.** Formulations
  carrying it proposed no rejected tuples. They were the fastest native-XOR
  models on both target sets, and without phases they proved UNSAT 4–5 times
  faster (91–118 s against 465–513 s for 4 targets). Single trials, 4
  targets each.
- **Native XOR against expanded clauses.** At 300 s, expanding XOR into
  clauses was slower in every paired comparison: `spec_cnf` against
  `spec_xor`, and `line_subgroup_cnf` against `line_subgroup`, by 1.4–3.8
  times. At 2 s the note saw no difference, and at that budget nothing
  finishes. This is not a robust factor.
- **Line against specialized.** There is no consistent winner. With the
  certificate the two are within noise; without it `line_xor` was faster on
  random targets and slower on planted ones.

- **Degree 19 is out of reach at these budgets.** With phases unknown, n = 13
  planted targets finish in a median of about 50 s. At n = 19, none of 8
  planted attempts finishes in 1200 s. That is at least a 24× increase for
  6 more field bits, from one lower-bounded sample per cell. Without phases,
  n = 19 planted attempts also time out at 300 s.

## Degree 29: not attempted

The note identifies n = 29 as the first informative sparse control, where
about 8% of target/label cases are compatible. No SAT solve was attempted
there, because the n = 19 planted instances already exceed 1200 s. A run at
n = 29 at comparable budgets would only record timeouts, which are not
evidence about the formulation. Revisit when an n = 19 planted instance
solves, either at a longer budget (resume with `--resume --budget ...`) or
with a better encoding.

## Not measured

Relation-matrix rank, Groebner solving degree, repeated timing trials, and
end-to-end discrete-log cost were not measured. The representative-only
n = 19, 300 s grid is partial (12 of 48 attempts) and can be resumed.
