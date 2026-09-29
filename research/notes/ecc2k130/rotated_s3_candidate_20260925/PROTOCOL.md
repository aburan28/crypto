# Frozen protocol: recursive-S3 affine-candidate diagnostic

Status: preregistered before this experiment's first candidate outcome. This
is a small-rung, solver-free algebraic-semantic diagnostic, not an S6/S7
exporter, a SAT/Gröbner benchmark, relation yield, or ECDLP speed claim.
It depends on the complete 16-point corpus in [PR #767](https://github.com/aburan28/crypto/pull/767)
and the independent signed-point semantic gate in
[PR #770](https://github.com/aburan28/crypto/pull/770). `FROZEN.json` pins
their source/input/archive hashes, this protocol, and both implementations.
Do not replace a point target or a failed arm after inspecting output.

## Hypothesis and finite reference

Use exactly the n13-m5 and n19-m6 beta=3, d=2 rotated factor slots and
eight subgroup Q points per arm in #770's `INPUTS.json`. Each Q has four
separate exact full-point targets Q+T, T in [O,(0,1),(1,0),(1,1)]. The field,
curve, rational lifts, and mask coordinates are exactly those of #767/#770.
The complete rational reference is every labelled signed factor-point tuple,
with all prefixes added by the group law. Reproduce all 64 archived Q+T
counts and all exact target mask sets before any candidate classification.

For every rational factor-x mask tuple, enumerate every affine recursive
S3 prefix path c2,...,c(m-1), where ck is the x-coordinate of the sum of
the first k factors. At each step solve **all** roots over GF(2^n) of

    A c^2 + B c + C = 0,
    A=(a+b)^2, B=ab, C=(ab)^2+1.

When A=0 and B!=0 solve the linear equation; A=B=0 has no root. When
A!=0 and B=0 use the unique square root. Otherwise normalize to u^2+u=h,
reject trace(h)=1, and enumerate u and u+1 using odd-degree half trace.
The terminal equation uses c(m-1), the last factor x and the fixed affine
target x. The x-only equation cannot distinguish a target from its negative;
the exact point oracle makes that distinction. Do not assign a fabricated x
to O. A chain with any O prefix is classified as an exceptional rational
witness and excluded from this affine-only candidate count.

The primary falsification question is whether an affine-only chain introduces
mask/path candidates without an exact Q+T signed-point witness, or misses any
all-affine signed witness. Keep separate (i) candidate paths with a nonrational
intermediate x, (ii) paths whose intermediates rationally lift but have no
consistent signed witness for the fixed full target, (iii) candidate masks
absent from the full oracle, and (iv) true masks supported only through O
prefixes. Also count exact point tuples with O prefixes and each identity,
inverse-to-O, doubling, and ordinary addition branch. A negative affine
candidate census is **not** a full PDP UNSAT unless exceptional branches are
also accounted for.

Apply the explicit rational-PDP trace-parity prefilter from #770 as a separate
counter: all mask tuples, accepted masks, rejected masks, affine candidate
paths before the filter, and paths after the filter. The all-four-T projected
problem retains both parity classes. Do not infer a 2x speedup or treat the
trace relation as automatically present in a raw Semaev ideal.

## Independent replay, caps, and decision

The producer uses the #762 polynomial-reduction/Euclid field implementation,
half trace, and direct recursive root enumeration. The verifier uses #767's
independent bit-serial/Fermat field and point law, solves u^2+u=h through a
GF(2)-linear pivot basis rather than half trace, and independently enumerates
all rational point tuples and all algebraic paths. It compares exact per-target
candidate path/mask lists, classifications, archived full-point counts, and
field-equation evaluations. Both programs check every returned root in S3.

Run n13 then n19 sequentially in separate cold children. Each producer and
verifier child has a 600-second wall and 512-MiB peak-RSS acceptance cap;
the runner retains stdout, stderr, exit, UTC interval, hashes, partial files,
and failure/censor receipts. A cap produces a censored arm, not a refutation;
do not substitute a nearby field or shrink the target set. Record field
multiplications/squares/inversions, root calls, path expansions, point
additions, input/output bytes, CPU, wall, and RSS. These are separate native
stage counters, not a common operation-equivalent S or matched-rho ratio.

The gate passes only if both arms finish within caps, independent arithmetic
and point counts match, every all-affine rational witness appears among
candidate paths, and every candidate classification is independently replayed.
A zero-spurious outcome would support this toy affine-chain semantics only;
the next exporter must still encode O/inverse branches, exact point lifting,
the n131 challenge field/model importer, and independently checked solutions.
If n19 exceeds its cap, preserve its censored receipt and accept no n19
candidate verdict. Update the canonical scoreboard with the measured stage
decision, leaving full-DLP cost and speedup unset.

From the repository root, run hash-only preflight before the first outcome:

```sh
python3 research/notes/ecc2k130/rotated_s3_candidate_20260925/ci_replay.py
```

After preregistration/review, run:

```sh
python3 research/notes/ecc2k130/rotated_s3_candidate_20260925/run.py --out /private/tmp/rotated-s3-candidate-run-20260925
```
