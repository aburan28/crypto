# Radical unordered-pair image: matched solver-stage experiment

This directory freezes an exact polynomial-presentation experiment for elliptic-
curve index-calculus decomposition systems.  It is a solver-stage engineering
diagnostic.  It does not measure a complete ECDLP, does not compare total cost
with Pollard rho and does not handle an unknown scalar.

## Result

For a factor-base coordinate set `W` of size `b`, replace an ordered pair
`(x,y)` by `(s,t)=(x+y,xy)` and impose the **radical vanishing ideal** of the
unordered-pair image.  Its degree-compatible regularity is exactly `b`; ordered
membership `<h(x),h(y)>` has regularity `2b-1`.

In sparse four-coordinate Semaev trees, this preprocessing lowers msolve's
maximum F4 degree on every one of the 34 matched targets.  All 276 timed
processes complete.  Every direct/symmetric root set, root recovery and
nondegenerate elliptic-curve group relation is checked independently.

The complete table is in [RESULTS.md](RESULTS.md).  Prime main and holdout cells
with `b >= 5` and the binary `k=2,3` cells improve the cold matched batch after
charging one-time radical preprocessing.  The tiny prime `b=3` batch regresses
after setup and remains in the table.  A prior binary `k=4` feasibility attempt
is retained as censored and supports no comparative claim.

## Boundary and classification

The exact algebraic floor is `d_reg >= b`: there are `b(b+1)/2` unordered
pairs and exactly that many monomials `s^i t^j` with `i+j <= b-1`.  The radical
presentation reaches ratio `d_reg / b = 1`.  This is a presentation result, not
an attack-complexity floor.

The direct sparse pair tree is the matched solver-stage reference.  Total common
operations, `S`, matched-rho ratio and generic-floor ratio are all **null**.
Classification: **engineering solver-stage diagnostic**.  No exponent or
end-to-end speedup is claimed.

The frozen [contract](contract.json) states the acceptance and falsification
conditions.  The repository WDSat regression is inapplicable because it executes
the Rust Boolean solver, while this experiment covers Sage/msolve prime-field and
radical-image systems.  Run 001 is an equivalent matched suite under the parent
[index-calculus accounting contract](../index_calculus_baseline_20260914/ec_index_calculus_contract.json).

## Reproduce

Requirements:

- SageMath 10.9;
- msolve 0.9.5 on `PATH`;
- one CPU thread; and
- enough space for approximately 10 MB of accepted raw evidence.

From the repository root:

```sh
DOT_SAGE=/tmp/radical-pair-sage \
  sage -python research/radical_pair_image_20260921/test_radical_pair.py

DOT_SAGE=/tmp/radical-pair-sage \
  sage -python research/radical_pair_image_20260921/run.py \
  --output research/radical_pair_image_20260921/results/new-run \
  --repetitions 3 --timeout 120

python3 research/radical_pair_image_20260921/summarize.py \
  research/radical_pair_image_20260921/results/new-run
python3 research/radical_pair_image_20260921/audit.py \
  research/radical_pair_image_20260921/results/new-run

python3 research/radical_pair_image_20260921/thread_scaling.py run \
  --run research/radical_pair_image_20260921/results/run-003 \
  --output research/radical_pair_image_20260921/results/new-thread-run
python3 research/radical_pair_image_20260921/thread_scaling.py verify \
  --output research/radical_pair_image_20260921/results/new-thread-run
```

The runner refuses to overwrite an output directory.  It freezes msolve inputs,
raw stdout/stderr, source and input hashes, target certificates, phase setup
times and all repetitions.  Execution order rotates by cell, target and
repetition.  `msolve -g 1 -v 2` retains the leading ideal and detailed F4
protocol without writing the much larger expanded reduced basis.

## Evidence

- Accepted run: [`results/run-003/`](results/run-003/)
- Summary: [`results/run-003/summary.json`](results/run-003/summary.json)
- Audit: [`results/run-003/final-audit.json`](results/run-003/final-audit.json)
- Exact thread scaling: [`results/thread-scaling-001/summary.json`](results/thread-scaling-001/summary.json)
- Superseded binary-certificate run: [`results/run-002/SUPERSEDED.md`](results/run-002/SUPERSEDED.md)
- Superseded absolute-path run: [`results/run-001/SUPERSEDED.md`](results/run-001/SUPERSEDED.md)
- Rejected output-mode trial: [`results/run-000-rejected/REJECTED.md`](results/run-000-rejected/REJECTED.md)
- Prior bounded binary `k=4` stop: [`results/binary-k4-censored.json`](results/binary-k4-censored.json)
- Canonical scoreboard: [`../../docs/index-calculus-scoreboard.html#radical-pair-image-20260921`](../../docs/index-calculus-scoreboard.html#radical-pair-image-20260921)

The broad use of elementary symmetric variables and elimination ideals is prior
art.  This package claims only the stated lemma, exact finite controls and
measured solver-stage behavior.
