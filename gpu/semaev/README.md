# gpu/semaev — the pairs-and-solve decomposition sweep

CUDA kernels for the decomposition oracle of
[`RESEARCH_SEMAEV_DECOMPOSITION.md`](../../RESEARCH_SEMAEV_DECOMPOSITION.md),
which records exactly this as deferred:

> **Parallelism.** The pair loop is embarrassingly parallel over `X₁` and
> nothing in it shares state. Four cores would be ~2 bits of `l`.
> Recorded, not done.

A GPU is ~10⁴ lanes rather than four, so the same argument gives ~13 bits
of `l`. That matters because `l = 16…18` is where the only question that
would change the note's conclusion — *is there a sub-`2^{2l}` oracle?* —
can actually be tested.

**No kernel here has been run on a GPU.** The environment this was
written in had no NVIDIA device and no `nvcc`. The arithmetic, the
subspace polynomial, the quartic and the sweep are verified on the CPU
against an independent Python oracle; see
[How this is tested](#how-this-is-tested).

## The algorithm, unchanged

`S₄` is a quartic in its last argument. Fix `X₁` and `X₂` and it becomes
a degree-4 polynomial `q(X₃)`; its four roots are the only candidates.
The factor base is an `F₂`-subspace `V`, and the polynomial vanishing
exactly on `V`,

```text
  L_V(t) = Π_{v ∈ V} (t + v) = Σ_{i=0}^{l} aᵢ · t^{2^i},
```

is **linearized**: degree `2^l`, but only `l + 1` non-zero coefficients.
So `L_V mod q` costs `l` squarings of a degree-3 polynomial, and
`gcd(q, L_V mod q)` has as its roots exactly `q`'s roots in `V`. No
search over the factor base at all — `O(l)` field operations per pair
where evaluating over the factor base would be `2^l`.

## Layout

| File | What it is |
|---|---|
| `sref.py` | Pure-Python oracle: the field, `L_V`, the quartic, and **exhaustive triple enumeration** as ground truth |
| `gf2n.cuh` | One-word `GF(2^n)` for `n ≤ 63`, branch-free |
| `decomp.cuh` | `L_V`, the quartic, `rem`/`sqr_mod`/`gcd`, and `decompose_row` |
| `test_cpu.cpp` | Verification harness — compiles the `.cuh` headers with g++ |
| `bench2.cu` | Device driver: selftest, throughput, occupancy |

## Build and test

```bash
make test            # needs only Python 3 and a C++17 compiler
make test L=4        # a different subspace dimension
make ladder          # l = 3, 4, 5 in sequence
make bench2          # CUDA binary (needs nvcc)
```

## How this is tested

| Check | What it covers |
|---|---|
| Field arithmetic | `mul`/`sqr`/`inv` against the oracle, then `sqr == a·a` and `a·a⁻¹ == 1` swept over the whole field (they are different code paths) |
| Subspace polynomial | `L_V` must vanish on **all** of `V` and on nothing else — swept completely at these sizes |
| Quartic coefficients | all five, for four `(X₁, X₂, x_R)` triples, against the oracle's independent derivation |
| **The sweep** | every target's verdict against **exhaustive triple enumeration**, plus every witness re-checked against `f₃` and against the factor base |

The last row is the decisive one, and it is the same gate the Rust module
uses (`decomposition_agrees_with_exhaustive_search`). A wrong gcd or a
wrong subspace polynomial does not merely run slower — it answers a
*different question*, and only comparison against exhaustive triples
catches that. Both verdicts must occur in every run; the harness fails if
they do not, so it cannot pass by always saying the same thing.

### What the CPU harness cannot cover

The launch structure. `./bench2 selftest` on real hardware compares both
kernels against the verified host sweep on every target, checks each
witness, and fails if only one verdict occurred.

## Two kernels, because the load is triangular

Thread 0 walks `2^l` pairs and the last thread walks one, so a static
row-per-thread mapping wastes about half the lanes.

- `sweep_kernel` — one thread per `X₁`, grid-stride, so a smaller grid
  re-balances.
- `sweep_kernel_flat` — one thread per **pair**, which removes the
  imbalance at the cost of redoing the quartic setup per pair.

Which wins is a hardware question. Both are checked for correctness;
neither is assumed faster.

## The one place a GPU is at a disadvantage

The CPU path uses `pclmulqdq` — a single carry-less multiply — and folds
the 128-bit product through a byte-indexed table. The note measures that
as roughly `6×` over a shift-and-xor loop.

NVIDIA hardware has no carry-less multiply, so `gf_mul` is the
interleaved shift-reduce loop: `n` iterations, branch-free (the
conditionals are arithmetic masks, so a warp never diverges on operand
values), but still `n` times the work of one `pclmulqdq`.

That is the honest counterweight to the lane count: a GPU brings ~10⁴
lanes and gives back ~6× per lane on the field's hottest operation. The
fix is **bit-slicing** the field across threads rather than packing it
into a word — a different data layout, not a tuning of this one, and not
implemented.

## One deliberate difference from the CPU path

`semaev_decomp::decompose` spends one field inversion per *row*, by
batch-inverting a whole row's leading coefficients (Montgomery's trick),
and uses pseudo-remainders so the gcd needs none. A batch inversion is a
sequential prefix product, so it does not belong inside a thread that
owns one `X₁`. Here each thread makes its own quartic monic with one
inversion per *pair*: more inversions in total, no row-wide dependency.
`./bench2 selftest` checks the two agree.

## What this does not claim

The note it comes from is explicit, and nothing here changes it:

> So: the oracle got roughly a thousand times faster and the attack's
> exponent did not move. That is not a disappointment, it is the shape of
> the problem.

Relation collection costs `Θ(2^n)` with **any** oracle polynomial in the
factor base, because a larger factor base needs fewer tries and makes
each try proportionally more expensive, and the two cancel exactly. A GPU
is another constant on that column, and under `AGENTS.md` it is an
**engineering** row: wall-clock is a practicality note, never the metric.

What it buys is the same thing the `926×` in the note bought — a usable
experiment. `l = 12` went from a day and a half per target to 15 seconds
and made that rung testable; this is the step that would make `l = 16…18`
testable, and that is the range where a sub-`2^{2l}` oracle could be
looked for. [`examples/fixed_x1_oracle.rs`](../../examples/fixed_x1_oracle.rs)
already measured one candidate route and found its deciding degree grows
as `l/2`, so the search is for a different route, not a faster version of
that one.
