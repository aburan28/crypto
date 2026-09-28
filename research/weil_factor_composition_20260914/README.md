# Weil-friendly composition of factor bases

This experiment preserves the factor base as a union of Frobenius translates,
but solves its components in their own small coordinates. It also investigates
which seed spaces keep **all mixed products** small. A small self-product alone
is insufficient. This is an engineering and structural experiment; calibrated
total cost, `S`, rho/floor ratios and attack exponents remain unmeasured.

The [contract](contract.json) was frozen before measurement. Its reference is
`c5d86171608af64353da1a2ade61f20538ee1527`, the F4 suffix-XOR implementation.
The native ambient F4, enumeration and pair-table routines remain the controls.
The external frozen WDSat S4 executable does not call this native S3 encoding;
the equivalent native regression reruns every previous frozen/fresh fixture.

## Construction and exact solver

Let the curve be defined over `GF(q)`, where `q=2^k`, inside `GF(2^n)`. For an
independent binary seed basis of dimension `ell`, form

```
V_i = Frob_q^i(V),  U = union_i V_i,  F = {P : x(P) in U}.
```

The old Frobenius-union constructor already produces this exact factor base.
For a nonlinear union its algebraic adapter uses the ambient `n`-dimensional
basis per summand, then checks membership while lifting. The new adapter retains
the seed coordinates and searches all `c(c+1)/2` unordered component pairs. It
uses `2*ell` Boolean variables per S3 system, with exact field equations and
point lifting. Every failed component and overlapping solution is retained in
the work accounting; the union is not silently replaced by one component.

For target abscissa `r`, the characteristic-two S3 equation is

```
(xy)^2 + r*xy + r^2*(x+y)^2 + b = 0.
```

For bases `(u_i)` and `(v_j)`, precompute `u_i*v_j`, their squares, and basis
squares. The target-dependent quadratic coefficients are
`(u_i*v_j)^2 + r*(u_i*v_j)`; the linear coefficients are `r^2*u_i^2` and
`r^2*v_j^2`. Inverse Frobenius sends the component pair `(i,j)` to `(0,j-i)`
and rotates the target. Because curve coefficients are fixed by `Frob_q`, only
`c` relative-offset tensors are needed. This reduces stored preprocessing, not
the number of component queries. Coordinates are mapped back with the original
ordered component bases, so canonicalizing a space never changes its coordinate
map accidentally.

The optional linear projection row-reduces the quadratic coefficient columns
first. Rows with zero quadratic part are **all linear consequences of that
coordinate row space**. A second elimination solves those linear equations,
rejects contradictions, and substitutes the pivot variables as affine functions
of the remaining variables. The existing Boolean F4/splitting solver handles the
remaining system. Affine reconstruction restores the original coordinates.
This is row-space preprocessing, not a claim to find every linear polynomial
in the full Gröbner ideal before solving it.

Both variants enumerate through non-lifting roots and share one node budget
across the entire component cover. Budget exhaustion is explicit. The initial
test caught an incorrect one-root limit per chart; that was fixed before any
accepted measurements. The identity target has no S3 abscissa and is handled
directly as a group sum. Full-DLP direct-scalar shortcuts are disabled.

The opt-in API is `WeilChartPlan::new` and `KoblitzIcOptions::weil_charts`.
It validates the field, curve, independent seed, Frobenius closure, exact union,
and two-summand arity. Defaults stay on the previous path. This experimental
packed-row implementation supports `n<=63`, `1<=ell<=7`, and `m=2`; it is not a
131-bit or higher-arity solver. Preprocessing is caller-owned and charged cold.
No new Redis dependency is required: immutable plans are reused in memory across
targets. This round does not implement serialization of chart plans into Redis.

For one complete example using the scaled GF(8) components on the degree-15
curve, run `cargo run --release --example weil_factor_composition -- dlp 8 charts_linear 503 7`.
The example builds the plan, passes it through the public index-calculus options,
checks the recovered scalar and emits the full relation-attempt verification log.

## Why scaled subfields are a useful composition

Take a proper subfield `K=GF(2^d)` and a nonzero scalar `alpha`, and use
`V=alpha*K`. Frobenius preserves `K`, so every component is a scalar copy of
the same subfield:

```
V_i = alpha_i*K,
span_F2(V_i*V_j) = alpha_i*alpha_j*K,
dim_F2 span(V_i*V_j) = d.
```

This holds for every relative offset, including mixed components. More precisely,
the quadratic coefficient space is the image of the product span `W` under the
binary linear map `L_r(p)=p^2+r*p`. For nonzero `r`, its kernel is `{0,r}`.
Consequently its dimension is exactly `dim(W)-1` when `r in W`, and `dim(W)`
otherwise. At `r=0` squaring is invertible and the dimension is `dim(W)`.
The linear columns and constant can still increase the rank of the whole system.
This is the source of cheap linear consequences for Weil restriction.

There is another way to see the opportunity. Write `x=alpha_i*u`,
`y=alpha_j*v`, with `u,v in K`, and introduce `w=u*v in K`. The field equation
becomes linear over F2 in the **3d coordinates of `(u,v,w)`**:

```
(alpha_i*alpha_j)^2*w^2 + r*alpha_i*alpha_j*w
    + r^2*alpha_i^2*u^2 + r^2*alpha_j^2*v^2 + b = 0.
```

The remaining requirement is `w=u*v`. If the linear system is inconsistent,
the chart is refuted immediately. If it has a unique solution, one multiplication
checks it. Otherwise only its remaining degrees of freedom require nonlinear
solving. The implemented row projection eliminates the product coordinates
implicitly; this auxiliary-variable explanation is a derivation, not a separate
measured solver. A small product space thus gives a precise mechanism beyond
merely reducing the number of variables.

There are two important controls. First, setting `alpha=1` may put the factor
base in a small point subgroup: cofactor projection can erase every useful
column. The measured column count and coverage, rather than a subspace dimension,
decide whether the base is useful. Second, the definition field and x-coordinate
subfield need not coincide. The frozen `GF(8)` curve over `GF(2^15)` with
`x in GF(32)` tests a complementary-subfield construction. This space is
Frobenius-stable over F2 but is not GF(8)-linear. It therefore explores a shape
outside the existing GF(8)-linear invariant-space search. Neither construction
guarantees adequate relation yield.

The composition should retain a **common** small subfield. Mixing two different
subfields has a different cost: the binary span of `GF(2^d)*GF(2^e)` is their
compositum `GF(2^lcm(d,e))`. For example, mixing GF(8) and GF(32) as summand
components fills all 15 field dimensions in their mixed product. Using GF(32)
abscissae on a GF(8)-defined curve is distinct from mixing both abscissa spaces.
This explains why scalar copies of one subfield are the first composition to
test, with the point-subgroup admission check retained.

Frobenius-invariant factor bases and subfield-curve index calculus are established
ideas; this round does not claim their invention. See Galbraith, Granger, Merz
and Petit, [On Index Calculus Algorithms for Subfield Curves](https://eprint.iacr.org/2020/1315).
The component-coordinate implementation, complete-cover comparisons and the
specific measured candidates are the contribution of this experiment.

## Prime-degree obstruction and counting boundary

For nonzero binary subspaces `A,B` of a prime-degree extension, linear Kneser
gives

```
dim span(A*B) >= min(n, dim A + dim B - 1).
```

The stabilizer of a proper product span cannot be a larger intermediate field
when `n` is prime. Thus the scaled proper-subfield construction has no nontrivial
instance at degree 131. This is a product-space bound, not an impossibility
theorem for fast index calculus. The field-extension form of Kneser's theorem
and its hypotheses are given by Bachoc, Serra and Zemor in
[Revisiting Kneser's Theorem for Field Extensions](https://arxiv.org/abs/1510.01354).

The structural script verifies the chosen degree-131 polynomial's irreducibility
before doing any field arithmetic. It measures every relative Frobenius product
for power spans, Frobenius spans and independent random spaces at dimensions
8, 16, 32 and 44. Symmetry of the relative ranks is checked independently. These
are exact field-linear-algebra diagnostics, with no point counts or DLP timing
extrapolation. Curve columns and coverage for these structural-only rows are
explicitly null; the small actual-curve census is in the stage archive.

The power-span failure already has a simple explanation at the first offsets.
For `V=span(1,z,...,z^(ell-1))`, its product with its square-Frobenius image
contains every power from `z^0` through `z^(3*(ell-1))`. At `ell=44`, these are
130 independent powers in a degree-131 field, versus only 87 for `V*V`.
At offset two, products contain every power through `z^215`, including an entire
field basis, so the rank is 131. This argument is independent of the chosen
irreducible polynomial. The measured remaining offsets check how broadly this
loss persists; they do not prove that every possible seed space must behave so.

With `B` signed factor-base points and a subgroup of order `N`, unordered pairs
give at most `B(B+1)/2` distinct sums. Consequently uniform nonzero subgroup
targets have coverage at most `min(1,B(B+1)/(2*(N-1)))`. Counting the same-component
pairs alone can reduce per-target work but loses coverage; its exact coverage is
recorded as an ablation and never represented as a complete oracle for the union.
Even a cheaper oracle must pay for the additional relation trials and matrix work.

## Protocol and reproduction

Ten cases include ordinary polynomial spans, independent random spans, scaled
proper subfields and unscaled controls. Every stage has eight distinct, uniformly
sampled nonzero subgroup targets; the target RNG is separate from seed-space
generation. Seed 17 is development and 937 is the fresh holdout. Complete S3
root sets are checked independently against exhaustive field evaluation. Each
first point witness is checked by group addition, and refutations are checked
against a separately built complete pair-sum truth set withheld from the solvers.

All ten cases enter the full-DLP comparison at four seed/secret combinations and
three repetitions, with ambient F4, projected charts, enumeration and pair tables.
Every relation-attempt input of a completed run is cross-checked against the
independent pair oracle. Solves must recover the scalar and verify `[k]P=Q`.
Zero-column, incomplete, timeout and low-yield controls are retained. Matched
signed-Frobenius rho runs cover each distinct curve and seed/secret combination.

The runner executes sequentially on one CPU, with 2 GiB of address space, one
thread and caches off. Variant order reverses on alternating repetitions. Native
regression uses separate reference/candidate binaries; unchanged baseline and
new chart variants are also compared in the same experiment binary, with identical
arithmetic and compiler settings. Source and binary hashes are preserved.
This was a shared host, not an isolated timing machine. Correctness-test
compilation ran on a separate CPU during part of the measurement; the benchmark
itself remained pinned to its one CPU. Bootstrap intervals are descriptive timing
diagnostics and repeated seeds are not independent algorithm instances.

Use Rust 1.90.0, release mode, no `RUSTFLAGS`, and copy the saved
`validation/dependencies.lock.txt` to `Cargo.lock`. Build the reference's
`f4_linear_algebra_bench` example in a separate checkout at the contract revision.
Build the candidate's two examples and copy binaries into separate directories.
Use separate Cargo target directories, or run `cargo clean --release -p crypto`
before switching worktrees. Then:

```sh
cargo build --release --locked --features redis-cache \
  --example weil_factor_composition --example f4_linear_algebra_bench
cargo test --release --locked --features redis-cache --lib weil_charts::
python3 research/weil_factor_composition_20260914/run.py \
  --reference-native /path/to/reference/f4_linear_algebra_bench \
  --candidate-native /path/to/candidate/f4_linear_algebra_bench \
  --experiment /path/to/candidate/weil_factor_composition \
  --output research/weil_factor_composition_20260914/results/new-run
python3 research/weil_factor_composition_20260914/structure.py \
  research/weil_factor_composition_20260914/results/new-run/structure.json
python3 research/weil_factor_composition_20260914/compare.py \
  research/weil_factor_composition_20260914/results/new-run
```

The raw JSONL archive stores every process's arguments, status, stdout and stderr.
Cold stage cost includes setup, target generation, field setup, eight first-witness
queries and their verification. Full root enumeration and truth-table construction
are separately timed validation. Cold DLP cost includes curve/base/plan setup,
target generation, all collection and matrix work, and final scalar verification.
Internal phase times are subsets of that cost and must not be added a second time.
Independent oracle cross-checks are validation, not algorithm inputs. Process wall
time also includes validation and output; it is preserved separately.

Tensor bytes count only the immutable coefficient payload, `16*c*ell*(ell+1)`
bytes in this single-word implementation. They exclude allocator metadata, field
tables, bases, coordinates, solver matrices and cached results. Operation counters
cover named kernels only; they are not total calibrated attack operations.

See [RESULTS.md](RESULTS.md) for every variant, retained failures and paired timing
intervals, and the [canonical scoreboard](../../docs/index-calculus-scoreboard.html#weil-factor-composition-20260914).
