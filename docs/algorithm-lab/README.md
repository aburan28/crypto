# Interactive algorithm lab

Open [algorithm-lab.html](../algorithm-lab.html) directly, or visit the page under
`/crypto/scoreboard/algorithm-lab.html` after Pages publication. The page loads only
local JavaScript and CSS and needs no server-side computation.

These are exact, bounded teaching examples, not performance measurements or a
new cryptanalytic method. Benchmark figures remain in the performance dashboard
and canonical cost ledger. No benchmark is rerun or runtime solver source changed by the lab.

## Factor bases and relation matrices

The curve is `y² = x³ + 2x + 2` over `F_17`. Its group has order 19 and generator
`G = (5, 1)`. The demo takes affine points with `y ≤ 8`, sorted by `x`, and retains
the first 2–9 as the factor base. It does not claim this small coordinate cut is
the subspace/orbit construction used by the repository.

Unordered pairs with repetition are enumerated and re-added in the group. Known
probes `[a]G` match those sums, giving `sum(c_i * ell_i) = a (mod 19)`. The matrix
engine receives only these coefficients and right-hand sides. Mixed pairs are
considered before doubles to produce instructive row operations; a maximal
independent subset is selected for the default walkthrough. This enumeration and
selection are teaching overhead, not a timed index-calculus result.

At size 4 the base is `(0,6), (3,1), (5,1), (6,3)`. Its logarithms are `7,4,1,2`.
The first pair sums to `(13,10) = [11]G`. Each recovered log is checked by point
multiplication. When the selected target decomposes, its scalar is reconstructed
from the recovered logs and verified again.

The duplicate-row mode removes one independent equation and duplicates the
first. The incorrect-relation mode changes one right-hand side by 1: it fails the
group check and also produces an inconsistent augmented matrix. Real relation
collectors reject that input before matrix solving.

The matrix demonstration uses dense Gauss–Jordan elimination over `F_19`; it does
not execute the repository's sparse filtering/block-Wiedemann implementation.
See [the repository workflow](https://github.com/aburan28/crypto/blob/main/docs/ic/README.md#linear-algebra-relation-filtering-and-block-wiedemann).

## SAT

The four variables are `x, y, z, w`. Each example defines `z = x AND y`, encoded as
`(not x OR not y OR z) AND (x OR not z) AND (y OR not z)`.

- **AND + XOR:** `x XOR y XOR w = 1` and `z XOR w = 0`.
- **Search:** `(x OR y) AND (not x OR y) AND (x OR not y)` and `z XOR w = 0`.
  Trying `x = 0` first causes a conflict and backtrack; `1,1,1,1` satisfies it.
- **Contradiction:** both `x XOR y = 0` and `x XOR y = 1`, plus `z XOR w = 0`.
  Every one of the 16 assignments fails.

The bounded DPLL interpreter performs CNF unit propagation, single-row XOR
propagation, and chronological backtracking. It has no clause learning,
restarts, watched literals, or Gaussian combination of parity rows. The real
Rust solver has CDCL and native parity reasoning. The Boolean teaching example
is separate from the prime-field curve example; it is not that curve's Semaev
encoding. In the real binary-field pipeline, SAT produces decompositions before
relation linear algebra.

CNF mode forbids each failing parity assignment with one clause. The native and
CNF modes preserve exactly the same original truth table. Downloads are ordinary
DIMACS CNF with `1=x, 2=y, 3=z, 4=w`, suitable for an external SAT solver. Trace
counters count only this interpreter's actions and are not runtime speedups.

## Validation

```sh
node scripts/site/test_algorithm_lab.js
python3 -m unittest discover -s scripts/site -p 'test_*.py'
# Optional real-browser checks (requires Playwright and installed Chromium):
python3 scripts/site/test_algorithm_lab_browser.py
```

The engine checks cover all 361 point additions against a fixed scalar-multiple
table, all 8 base sizes and 144 nonidentity target selections, every row-operation
invariant, rank deficiency and inconsistent relations, all 16 assignments for
each SAT encoding, both branch orders, and the exported DIMACS meaning. The
point table and default matrix are also checked independently with SageMath.

References: [SageMath finite-field elliptic curves](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_finite_field.html),
[CryptoMiniSat](https://github.com/msoos/cryptominisat), and the repository's
[Semaev encoding](https://github.com/aburan28/crypto/blob/main/src/cryptanalysis/semaev_sat.rs)
and [SAT solver](https://github.com/aburan28/crypto/blob/main/src/cryptanalysis/sat.rs).
