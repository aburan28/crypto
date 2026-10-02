# Direct polynomial-basis sigma map: frozen static screen

Source parent: `efac7fc1f1a0eac8271313b3e9bea8194fe78f01`.

Class: compatible arithmetic engineering.  This screen does not run a walk,
search for distinguished points, solve a collision, or recover a key.  It asks
whether the two polynomial/normal basis conversions and the normal-basis
Frobenius network in the sigma walk can be replaced by one exact linear map.

## Frozen target

For every reduced 131-bit polynomial-basis input `p` and every `j` in 3..10,
the candidate must compute

```
L_j(p) = toPolynomial131(
             fromPolynomial131(p)
           + sigma131(fromPolynomial131(p), j))
```

with bit-for-bit equality.  The intended call site applies the same selected
map independently to X and Y.  The Hamming-weight selection remains in the
normal basis and is outside this experiment.

## Native construction and oracle

The experiment uses one C++17 program; Python and generated Python artifacts
are excluded.  It obtains each 131-by-131 matrix from the repository's actual
`fromPolynomial131`, `sigma131`, and `toPolynomial131` implementations.  The
program must:

1. check all 131 input basis vectors for all eight powers against the composed
   repository oracle;
2. check every emitted candidate on those same 1,048 map/basis cases;
3. check deterministic dense and edge vectors for each power; and
4. report the matrix rank and reject non-canonical high output bits.

Because every implementation under test is GF(2)-linear, equality on all
basis vectors proves equality on every 131-bit input.  Dense vectors are a
separate implementation guard, not additional mathematical coverage.

## Circuit families and accounting

At least one complete direct circuit family must be emitted and compiled.
The initial families are:

- a constant-table family over fixed input chunks, with table bytes, loads,
  XORs, masks, and dynamic-index costs reported separately;
- a word-level masked-shift family assembled from matrix diagonals, with every
  32-bit shift, mask, OR, XOR, load, temporary, and output word counted; and
- an in-place masked-CNOT/elimination family if the native synthesizer can
  factor the rank-130 map without increasing the first family's ALU count.

Counts are deterministic source operations, not throughput predictions.
Apple Clang 17 `-O3` host assembly instruction and text-byte counts are a
second static screen.  They are reported for identical noinline wrappers and
fixed `j`; host assembly is not described as Blackwell SASS.

The baseline is measured twice: the shipping 9-word
`fromPolynomialProduct131` path and the existing reduced-input inverse path,
each followed by `sigmaWalkNetwork131`, addition, and `toPolynomial131`.
Memory operations and ALU operations remain separate because prior GPU work
shows that removing read-only loads need not track throughput.

## Decision gate

A candidate is eligible for a later GPU A/B only if it passes the full native
oracle and one of these predeclared static gates:

- at least 20% fewer compiled host integer instructions than the reduced-input
  composed path, without more than 64 KiB of total constant data for all eight
  maps and without data-dependent global-memory reads; or
- at least 30% fewer deterministic word-level ALU operations than that path,
  with no more read-only loads than the current 56 shared-mask loads per
  coordinate and no increase in live 32-bit temporaries above 16.

Otherwise the result is retained as negative static evidence and no GPU is
launched.  Passing the static gate authorizes only a separately frozen,
matched GPU protocol; it is not itself performance evidence.
