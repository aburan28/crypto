# Fixed eight-bit matrix-F5 and isogeny-neighbor measurements

This is an executable **toy algebra diagnostic**, bounded to GF(256) and at most
six Boolean variables. It measures the repository's Boolean matrix-F5 **row
criterion**, not a full incremental, signature-compatible F5 implementation.
It neither accepts cryptographic targets nor recovers discrete logarithms.

```sh
python3 -m unittest discover -s research/toy_f5_neighbors_20260924 -v
python3 research/toy_f5_neighbors_20260924/experiment.py --output /tmp/toy-f5-run
python3 research/toy_f5_neighbors_20260924/verify_replay.py /tmp/toy-f5-run
```

Python 3.10+ and the standard library suffice. The frozen run uses approximately
33 seconds here; this is a practicality observation, not a performance claim.
The raw artifact is gzip-compressed JSON. Each system retains its generators,
roots, per-degree matrix traces, row-space hashes, verified lifting counts,
transported inputs and three repetition hashes. Repetitions are deterministic
counter checks, not independent statistical samples.

## Fixed protocol and boundaries

`contract.json` was written before the measurement run. The reference and
candidate use the same subgroup points, selector labels, support size and
transported targets. Frozen seed 101 and holdout 503 each select four nonzero
x-coordinates (eight signed points). Two or three summands require four or six
Boolean selector bits. Both ordered and canonical encodings are measured on all
32 subgroup targets, including infinity, on the source and all four neighbors.
That is 1,280 systems and 7,680 solver runs with three repetitions.

The canonical variant adds an equation requiring nondecreasing selector indices.
It keeps one representative per x-multiset; all signed lifts are still checked.
It is a symmetry-breaking constraint, not an invariant-coordinate elimination
or a proof that quotient equations always reduce degree. Repeated summands and
opposite signs are allowed. Support excludes the order-two point.

The predeclared test requires every neighbor to reduce aggregate F5 matrix XORs
by at least 10% on both splits and encodings, without any increased paired
completion degree. It fails. See `RESULTS.md`. Native-neighbor factor bases and
isomorphic-coordinate controls are not measured; these are exactly transported
workloads, not independent curve samples or evidence that conductor causes a gain.
The GF(128) enumeration pilot remains a separate, unchanged experiment. GF(256)
is used here because the source has four rational cyclic order-three kernels.

The boundary for the present table is the corresponding source matrix stage.
There is no measured full-DLP reference or calibrated universal lower bound for
this stage. Full-DLP operations, S, rho/floor ratios and speedup are **null**.
The WDSat regression protocol is inapplicable to this fixed Boolean reference;
this matched frozen/holdout suite is its bounded equivalent for stage diagnostics.
No exponent, asymptotic or end-to-end improvement is claimed.

## What degree and cost mean

The Boolean ring is F2[x1,...,xn]/(xi^2+xi). At budget D the matrix contains every
s*f_i with deg(s)<=D-deg(f_i). The code checks equality of F4 and F5 reduced row
spaces at **every** budget. Completion means the span contains every generator
and is closed under multiplication by all variables. This certifies the whole
Boolean ideal independently of the exhaustive truth-table oracle. A second
certificate checks rank = 2^n - number of roots and verifies every root.

The reported number is the first **Macaulay ideal-completion budget**. It is not
intrinsic degree of regularity, a full-F5 solving degree, or the largest reduced
polynomial degree (budgets can exceed n). F4 and the F5 filter must have identical
completion degrees for the same input generators. Neighbor differences describe
how their generators present the ideal, not a change in its transported zero set.

The criterion is the one documented in `src/cryptanalysis/matrix_f5_f2.rs`:

    W_i(D) = V_(i-1)(D-d_i) + span{s*(f_i+1): deg(s)<=D-2*d_i}
    V_j(e) = span{s*f_k: k<=j, deg(s)<=e-d_k}

A row t*f_i is pruned only when t is a leading monomial of W_i(D). Prefix
elimination is reused per distinct lower budget. Frobenius syzygies use f_i+1;
using f_i itself is unsound in the Boolean ring, and a regression test preserves
that counterexample. The reference uses degree then integer-mask monomial order.
The Rust implementation was inspected, but **not executed** (Rust is unavailable
in this environment). Its source SHA256 at implementation was
`e10c2290b40750fbf65e9ba8e65193d1f253939d224b0e3a3a584fc0bd9fa381`.
This is an independent bounded implementation, not a native Rust benchmark.

Each matrix row fits in one 64-bit word. The table charges elimination **plus
criterion** XORs over the whole degree sweep; it excludes the separately logged
completion-certificate XORs. Symbolic products, packing, field-table setup,
point enumeration, encoding, extraction and verification are additional work.
They have no calibrated conversion to this unit. Consequently the table is a
matrix-stage ledger, not total solver cost. Encoding records field-equation
evaluations and Mobius-transform XORs; verification records exhaustive signed
lift checks. Isogeny certificates record their checked point/pair counts.

## Actual equations and exact maps

`curves.py` fixes modulus 0x11B and the family y^2+xy=x^3+b. It evaluates Semaev
S3 and the quadratic resultant defining S4; it does not construct equations from
the planted solutions. Output field bits are interpolated to Boolean ANF on the
<=64 selector assignments. This exhaustive interpolation is explicit setup,
not a scalable Weil-descent encoding. Infinity targets use the appropriate
lower summation polynomial. Every recovered x-tuple is lifted and verified by
independent point addition; extraneous resultant roots would be rejected and
counted (none occurred in this suite).

For a kernel {O,Q,-Q}, q=x(Q), the normalized degree-three map has

    X = x + q/(x+q) + q^2/(x+q)^2
    Y = y(P) + y(P+Q) + y(P-Q)
    b' = b + q + q^2.

The degree-three numerator x^3+(q^2+q)x evaluates to q^2!=0 at its denominator's
root. We check all rational images, the exact kernel, the rational x formula,
1,024 subgroup homomorphism pairs, and both dual compositions on 288 points per
side. The normalized dual here needs negation to compose to [3]. Restricting to
the order-32 subgroup makes the degree-three maps injective; this tiny 2-power
group is not a model of a cryptographic prime-order subgroup.

The source descends from F2, where trace -1 gives fundamental discriminant -7
and maximal endomorphism order. Over GF(256), trace -31 gives -63=-7*3^2.
Since 3 is inert in Q(sqrt(-7)), the four degree-three edges from the maximal
order descend to conductor 3. Rational three-torsion drops from 9 to 3 on every
codomain. These algebraic facts supplement the finite map checks; point checks
alone would not prove an arbitrary rational formula is an isogeny.

References: [Semaev summation polynomials](https://eprint.iacr.org/2004/031),
[Sutherland, Isogeny volcanoes](https://arxiv.org/abs/1208.5370), and the F5
criterion's derivation and references in the repository source above.
