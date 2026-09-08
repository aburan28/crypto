# Bit-sliced Weil descent: review of EC-Index-Calculus-Benchmarks

**Upstream:** [`mtrimoska/EC-Index-Calculus-Benchmarks`](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks)
(~2k lines of C + a standalone Sage reimplementation, 2018)
**Audit scripts:** `research/ec-index-calculus-review/`
**Write-up:** <https://claude.ai/code/artifact/25b658b7-4d69-4c3c-82ad-988a50dae691>
**Provenance:** Review of the prior art closest to
`src/cryptanalysis/semaev_sat.rs`, which was stuck at `n = 4`–`5` while
that repository reaches `n = 19` on the same class of system.  The
question was what it does differently; the answer is in
["What we took from it"](#what-we-took-from-it), and our pipeline now
reaches `n = 19` too.

## What it is

An instance generator for a SAT-based index calculus attack on binary
elliptic curves: Weil descent of the **symmetrised fourth summation
polynomial** over `F_{2^n}`, for the Koblitz curve
`E: y² + xy = x³ + x² + 1`.  Emits every instance in three formats —
Magma/Gröbner, XOR-CNF, and a custom ANF.

Verified locally: the toolchain builds with `gcc *.c -o weill`, and
re-running `script_benchmarks.sh` regenerates the shipped XOR-CNF, ANF
and Magma artifacts **byte-identically**.

## The core idea: one polynomial, one bit-vector

An entire multivariate Boolean polynomial lives in a single flat
bit-vector.  Bit `d·T + k` means *"the coefficient of `e^d` contains
monomial number `k`"*, where `T` is the number of monomials in the
coefficient ring.  `params.c` supplies an O(1) rank/unrank pair
(`indices_to_num` / `num_to_indices`) between monomials and integers,
ranking a monomial as a multiset of variable *groups* plus a
mixed-radix digit within each group.

Three consequences, in decreasing order of cleverness:

1. **Frobenius is a bit-shift.**  `power()` computes *even* powers with
   zero multiplications — shifts and masks only.  This needs two facts
   simultaneously: char 2 gives `(Σ c_d e^d)² = Σ c_d² e^{2d}` with no
   cross terms, and the Boolean field equations give `a² = a`, hence
   `c_d² = c_d`.  So squaring leaves every coefficient's bit pattern
   untouched and merely relocates it.  Odd powers cost exactly one
   multiply.  Of the twelve terms of the descended `S_4` — including
   `e₁⁴`, `e₂⁴`, `e₃⁴`, `e₃³` — only four need a multiply at all.

2. **Constant multiplication and reduction are shift-XOR.**
   `multiply_with_const` multiplies by the known `x_R ∈ F_{2^n}` by
   walking that constant's set bits with a *delta-encoded* shift: one
   running accumulator, shifted by the gap to the next set bit rather
   than re-shifted from scratch.  `modulo` reduces mod the field's
   irreducible the same way.  Reduction is **lazy** — degrees grow to
   the exact bound `_n = 4n + 8l − 11` and `modulo` runs once, after
   all twelve terms have accumulated.

3. **One engine, two instantiations.**  `set_params` reconfigures the
   same `multiply`/`power`/`modulo` for both halves of the descent:
   stage A in the x-monomial ring (degree ≤ 3, `_n = n`) for the
   `X ↔ e` correspondence; stage B in the e-monomial ring
   (degree ≤ 2, `_n = 4n + 8l − 11`) for the Semaev part.

## Secondary techniques

- **Free monomial elimination.**  `compute_offsets_X/e` OR-reduce the
  bit-vector against itself shifted by `T`, collapsing every degree
  block onto block 0, which yields the set of monomials occurring
  *anywhere*.  Monomials that never occur get no SAT variable; the
  survivors are renumbered contiguously.

- **ANF instead of Tseitin.**  The `.anf` format is DIMACS-like with
  inline degree markers on XOR lines:
  `x T 25 .2 1 7 .2 1 13 .2 7 13 0` means
  `25 ⊕ (1∧7) ⊕ (1∧13) ⊕ (7∧13) = 0`.  No auxiliary variable per
  monomial.

- **Staged codegen through the C compiler.**  `create_semaev` writes
  `terms.h`; `compute_offsets_*` write `offset.h` / `e_offset_weil.h`;
  `script_benchmarks.sh` re-runs `gcc *.c` *between* phases so each
  stage's results become compile-time constants for the next.
  Non-reentrant, six compiler invocations per instance.

- **Planted-solution benchmarks.**  SAT instances are built by
  constructing a decomposition (sum three factor-base points), so a
  solution provably exists and is recorded in the `INFO` file.  The
  generator asserts `f3 == 0` before writing, and the Sage script is an
  independent second implementation of the same modelling.

## Measured: what the encoding choice costs

Instance `n19l6-1-S`, one system, three formats:

| Encoding | Variables | Clauses | Bytes |
|---|---:|---:|---:|
| ANF (`.anf`) | 51 | 52 | 37,946 |
| XOR-CNF (`X….dimacs`) | 767 | 2,416 | 45,476 |
| Plain CNF (`….dimacs`) | 4,986 | 19,444 | 324,032 |

The ANF→plain-CNF byte ratio is stable across every family in the
repository: **8.4×** at `n15l5`, **8.5×** at `n17l6`, **8.5×** at
`n19l6` — a property of the encoding, not of one instance.

## Audit: one shipped instance is mislabelled

The `U` instances are generated from a *random* `x_R`, which is
unlikely but not guaranteed to be undecomposable over the factor base.
The label is an expected outcome, not a certificate.

`research/ec-index-calculus-review/audit_labels.py` decides all 30
exhaustively (every sorted triple in the `l`-dimensional subspace;
complete because `f₃` is symmetric in `X₁, X₂, X₃`).

- All 30 planted SAT solutions verify.
- 29 of 30 `U` instances are genuinely unsatisfiable.
- **`INFOn19l6-19-U` is satisfiable**, with
  `X₁ = a⁴ + a²`, `X₂ = a⁵ + a`, `X₃ = a⁵ + a³ + 1`.

`verify_witness.py` confirms this against the shipped
`Xn19l6-19-U.dimacs` itself rather than only against our own algebra:
propagating the assignment leaves 0 variables unassigned and 0 clauses
violated, with one positive control (a planted solution on a genuine
`S` instance) and two negative controls.

Anyone using these 60 files as a regression corpus should expect one
SAT answer where the filename says `U`.

## What we took from it

All of the below is implemented; see
[`RESEARCH_SAT_SEMAEV.md`](./RESEARCH_SAT_SEMAEV.md) for the details
and the measured numbers.  The pipeline went from stalling at
`n = 4`–`5` to deciding the `n = 19, l = 6` instance end to end.

1. **Factor-base restriction** — `weil_descend_s3_subspace` confines
   each unknown to the `l`-dimensional subspace, turning `n` equations
   in `2n` unknowns into `n` equations in `m·l ≈ n`.
2. **Symmetrised `S₄` in char 2** — `binary_semaev_s4`, with squaring
   implemented as coefficient relocation rather than multiplication.
   Specialised to `b = 1`; general `b` still needs re-deriving.
3. **XOR-native Gauss-Jordan** in `sat.rs`, hooked at the `propagate()`
   fixpoint.  `analyze()` and `backjump()` were untouched, because a
   reduced row yields a reason clause directly.
4. **Native parity constraints instead of Tseitin** —
   `XorEncoding::Native` vs `XorEncoding::Cnf`, both kept so the
   difference stays measurable.  It is 10.6× the clauses at `n = 19`.
5. **Bit-sliced coefficients: deliberately deferred.**  It buys descent
   time, and descent is 11 ms against a four-minute solve.  The
   algorithmically important half — Frobenius as relocation, and lazy
   reduction — is implemented; the bit-packing is not.
6. **Reference corpus** — `semaev_corpus.rs`, parameters only (see
   below), with `n19l6-19-U` recorded as satisfiable.

Also fixed: `collect_xor_lits` no longer smuggles a parity flag in-band
as literal `0`, which is the DIMACS clause terminator.

### An unexpected cross-check

Our encoder, written from the algebra rather than from upstream's code,
produces **exactly** upstream's instance size for the `n = 19, l = 6`
family: 767 variables, 2 364 ordinary clauses, 52 parity rows, against
their `p cnf 767 2416` with 52 `x`-lines.  Two independent
implementations agreeing to the clause is a much stronger check on the
modelling than any single instance solving.

### On vendoring

The corpus module carries instance *parameters*, not upstream's
generated DIMACS/ANF files.  Upstream is GPL-3.0 and this crate
declares no licence, so copying its output would be a licensing
decision rather than an engineering one; regenerating from parameters
sidesteps that and exercises more of our own pipeline besides.

## Landmines in the upstream code

| Where | What | Bites when |
|---|---|---|
| `multiply()` | Accumulates with `vect_bin_set_1` (OR, not XOR), so duplicate monomials never cancel — not a correct general `GF(2)` multiply | Safe at every current call site (factors have disjoint variable groups, or provably multiplicity-free products); the invariant is undocumented |
| `semaev_out.c` | Loop variable `i` clobbered inside the `k` loops by `num_to_indices` | Harmless only because `eᵢ` is homogeneous of degree `i`, so the clobbered value equals the original — extending to `S_5` breaks it silently |
| `degree()` | Scans from `_vect_bin_size`, one past the last valid bit; indexes `t[-1]` when `_n·T ≡ 0 (mod 64)` | Not triggered by any shipped parameter set |
| `CMakeLists.txt` | `add_subdirectory(src)`, and there is no `src/` | Immediately; only the `gcc *.c` path builds |
| `…to_unsym_out()` | `first = 1` reset inside the inner loop, so `+` separators are never emitted and the Sage output is malformed | Never — dead path the pipeline does not assemble |
| `init(v)` | Every vector is a fixed `__ARRAY_SIZE__` stack array (472 KB at the shipped constant), several live per function | Only the live prefix is zeroed, so cheaper than it looks, but the frames are real |

## Reproducing

```bash
git clone https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks.git
python3 research/ec-index-calculus-review/audit_labels.py <clone>
python3 research/ec-index-calculus-review/verify_witness.py \
    <clone>/benchmarks/Xn19l6-19-U.dimacs 6 20,34,41
```

## References

- M. Trimoska, S. Ionica, G. Dequen, *A SAT-based approach for index
  calculus on binary elliptic curves* — the work this artifact
  accompanies.
- P. Gaudry, *Index calculus for abelian varieties of small dimension
  and the elliptic curve discrete logarithm problem*, 2009.
- J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using symmetries in
  the index calculus for elliptic curves discrete logarithm*, 2014 —
  the symmetrisation this repository applies to `S_4`.
- M. Soos, K. Nohl, C. Castelluccia, *Extending SAT solvers to
  cryptographic problems*, SAT 2009 — XOR-native reasoning.
