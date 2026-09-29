# Frozen protocol: projected rotated m5/m6 PDP corpus admission

Status: preregistered before this corpus oracle or any solver run. This is a
point-target corpus admission, not a solver benchmark, relation-collection
result, rank result, or ECC2K-130 discrete-log claim. It depends on the
normal-basis construction and independently replayed n13 support in
[PR #762](https://github.com/aburan28/crypto/pull/762), and the separate
solver-interface admission in [PR #763](https://github.com/aburan28/crypto/pull/763).
No solver is timed here.

The solver-admission #763 prerequisite is merged at
`ec990d85b61eacccefa9dee97f2a835e7a20b785`. The parent #762 evidence is merged on main at
`2de218f583d2dddeafdd0180dd02329dc57d53e4`; its immutable
`evidence/raw.tar.gz` has SHA-256
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`.
The n13 model and beta are inherited from its frozen source ancestor
`a6be357281fd03dab15ce690d46e4e571209febf`.

## Frozen arms and reason for the second rung

Use `E: y^2+xy=x^3+1`, all rational lifts of each x, and no point at
infinity in a factor set. The two arms are exactly:

| arm | field polynomial | q | beta | m | d | lambda |
|:--|:--|--:|--:|--:|--:|--:|
| n13-m5 | `x^13+x^4+x^3+x+1` (`0x201b`) | 2003 | 3 | 5 | 2 | 89 |
| n19-m6 | `x^19+x^5+x^2+x+1` (`0x80027`) | 130873 | 3 | 6 | 2 | 41811 |

For each arm form `V_i=span_F2{beta^(2^i), beta^(2^(m+i))}` for
`0<=i<m` and `F_i={P in E(F_(2^n)):x(P) in V_i}`. Verify Rabin
irreducibility, normal-basis rank n and trace 1, rank of all m slices
equal to 2, rank of their sum equal to 2m, and `tau(F_i)=F_(i+1)`.
The n19 polynomial is the smallest-mask sparse irreducible used by the
repository's Koblitz constructor. These polynomial, beta and arm choices
were made by field-only preflight before seeing n19 support outcomes.

Independently verify `#E(F_(2^13))=8012=4*2003` and
`#E(F_(2^19))=523492=4*130873` from `t_0=2,t_1=-1,
t_n=-t_(n-1)-2*t_(n-2)`, with both q prime. Select `G` as the first
affine point in increasing `(x,y)` whose `H=[4]G` is nonzero. Verify
`[q]H=O`, `tau(H)=[lambda]H`, `lambda^2+lambda+2=0 mod q`, and
`lambda^n=1 mod q`. The n13 H must match PR #762's `[4793,2429]`.

PR #762 found projected misses at n13 m5 d2 and saturation at n13 m6 d2.
The n19 m6 d2 arm is chosen before its support measurement because each
`F_i` has at most `1+2(4-1)=7` points, hence at most `7^6=117649`
labelled tuples. Since `q=130873`, at least `13224` projected subgroup
targets must be unsupported even without accounting for collisions. This
is a necessary counting guarantee, not a predicted hit rate or solver cost.

## Complete oracle and deterministic target freeze

For each arm, enumerate **every labelled tuple** in `F_0 x ... x F_(m-1)`
by exact group law. Save the full point-sum histogram, its `[4]`-projected
histogram, multiplicities, one witness for every supported sum, factor
points, and operation/timing receipts. The complete oracle decides the
eight selected target labels; a sample of tuples is not a negative proof.

The selection domain is `ECC2K130-ROTATED-PDP-CORPUS-20260925-v1`.
SHA-256 digests are read big-endian. For planted targets, scan counters
`0,1,...` with seed `{domain}/{arm}/planted/{counter}` and reduce modulo
the exact labelled tuple count. Decode in `itertools.product` order.
For each resulting full sum `S`, set `R=[4]S` and
`Q=[4^(-1) mod q]R`; skip `R=O` and duplicate `R`. Keep the first four
distinct `R`. Record `S`, the tuple, `R`, `Q`, and the unique 4-torsion
`T=S-Q`; independently verify `S=Q+T` and `[4]Q=R`.

For negatives, scan seed `{domain}/{arm}/negative/{counter}`, set
`k=1+(digest mod (q-1))`, `Q=[k]H`, `R=[4]Q`, and skip duplicate `k`
or any `R` in the **complete** projected histogram. Keep the first four
distinct absent `R`. Preserve every attempted counter/k and the exact
support membership decision so no favourable target can be hand-picked.
For each accepted positive and negative Q, record all four `Q+T` full
point multiplicities and their witnesses or exact absence. Stop with a
preserved failure receipt if either class has fewer than four distinct
targets after 100000 counters; never shrink or replace an arm.

Freeze the resulting sixteen point targets, source and input SHA-256s,
factor lists, full histogram, and accepted/attempted selection stream in
the PR before passing the targets to any solver. A later solver run must
read this committed corpus byte-for-byte and use a new preregistered
solver protocol. No target is selected by SAT, FES, crossbred, F4/F5,
Msolve, WDSat, or solver wall time.

## Verification, caps and decision

An independent bit-serial field implementation with Fermat inverses
rebuilds all factor points and **directly enumerates every labelled
tuple**, rather than trusting the producer's accumulated histogram.
It compares both full and projected multiplicities, verifies every
planted witness and each `Q+T` identity, independently reproduces the
SHA selection stream, and confirms every negative through the complete
projected histogram. The parent #762 source/input hashes and its n13
resulting support archive SHA are pinned in the evidence receipt once
that PR is merged; this corpus does not alter those inputs.

Charge shared field/basis/generator setup, factor lifting, all group
additions and scalar multiplications, complete histograms, all SHA
candidates, raw serialization, and independent verification. Report
operation counts as primary, wall/CPU/peak RSS as secondary. Each
producer arm has a 300-second wall acceptance cap and 512-MiB peak-RSS
acceptance cap; each verifier arm has a 600-second and 512-MiB cap.
`SIGALRM` stops a wall overrun; peak RSS is a post-measurement acceptance
gate. Preserve stdout, stderr, exit status, partial files and failure
receipt on a cap or validation failure. Do not substitute a nearby rung.

Admit the fixed corpus only if both arms have exactly four planted and
four exact projected-negative targets, all source/input hashes and
independent checks pass, and each arm meets its caps. Otherwise report
the failed gate and its preserved evidence. No wall-time speedup or
ECC2K-130 extrapolation follows from a successful corpus admission.

## Exact source and input freeze

Before any corpus support outcome, the source, input manifest, merged parent
source and archive were pinned in `FROZEN.json`. The corpus `FROZEN.json` SHA-256 `7e0ed24783547ac8a80a83cd62a27976f72798730483121db0c2e55341e45e6b`.
`run.py` and CI check this independent protocol anchor before any target is
measured or replayed. Use Python >=3.12 from the repository root:

```sh
python3 research/notes/ecc2k130/rotated_pdp_corpus_20260925/ci_replay.py
python3 research/notes/ecc2k130/rotated_pdp_corpus_20260925/run.py --out /private/tmp/rotated-pdp-corpus-run-20260925
```

The hash-only `ci_replay.py` invocation reads no new corpus support outcome.
When `evidence/raw.tar.gz` is committed, the same command with `--evidence
evidence-path` independently re-enumerates both complete point histograms.
