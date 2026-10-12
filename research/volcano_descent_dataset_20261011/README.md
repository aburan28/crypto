# Tall isogeny volcanoes: descent algorithm, dataset and per-layer ECDLP measurements

**Date.** 2026-10-11. **Scope.** Prime-field ordinary curves, toy sizes (fields of 40 to
68 bits). Binary fields: one ECC2K-130 step only (see `ecc2k130_level1.py`). Nothing here
concerns a deployed curve, and no ECDLP speedup is claimed. Results are in
[`RESULT.md`](RESULT.md).

**Conductor.** `conductor check` could not reach a control plane (connection refused on
localhost:8080), so this work ran on the hook's warn-through path. Scope is confined to this
directory.

## What was built

| File | Purpose |
|---|---|
| `volcano_descent.py` | The algorithm. CM parameter search, exact BFS from the crater, blind descent to the floor, deterministic curve and group-structure routines. |
| `build_dataset.py` | Deterministic campaign of 102 tall volcanoes, fail-closed checks, writes `data/*.csv.gz`. |
| `ecdlp_rho.py` | Counted-operation Pollard rho with Pohlig-Hellman. Modes `plain`, `neg`, `aut` (full automorphism orbit). |
| `run_ecdlp.py` | Experiments E1 (generic layers) and E2 (automorphisms, D_K = -3, -4). |
| `analyze.py` | Tables for `RESULT.md` and a leakage-aware ML baseline. |
| `ecc2k130_level1.py` | One real descent step on ECC2K-130 over F_(2^131) with ell = 263. |
| `test_volcano_descent.py` | 11 unit tests, about 10 seconds. `python3 -I test_volcano_descent.py` |

Dependencies: Python 3.13, `cypari` (PARI 2.15.4), and for `analyze.py` numpy, pandas,
scikit-learn, tabulate.

## The algorithm

Let E / F_p be ordinary with trace t and `t^2 - 4p = D_K f^2`. For a prime ell, the rational
ell-isogeny graph component of E is a volcano of depth `k = v_ell(f)`. Level i holds the
curves whose endomorphism ring has conductor `c * ell^i`, where `c = f / ell^k`. Neighbours
of a curve are the F_p-roots of the modular polynomial Phi_ell(X, j). A non-floor curve has
ell + 1 rational isogenies (one up, ell down; at the crater the "up" ones are horizontal).
A floor curve has exactly one.

* **Exact exploration** (`explore_from_crater`). Build the curve by CM so its crater is
  known: the roots of the Hilbert class polynomial of `D_K c^2` mod p. Run a multi-source BFS
  over Phi_ell roots. BFS distance to the crater set equals the level, so every vertex gets a
  ground-truth label.
* **Blind descent** (`descend_to_floor`). Needs only a curve, its trace and ell. Classify the
  first step by one non-backtracking random walk per neighbour. A descending walk ends at a
  degree-1 vertex in exactly `k - i` steps, any other walk takes at least `k - i + 2`. Then
  step to any neighbour other than the one just left until none remains. Kernels are counted
  with multiplicity, which matters at `j = 0, 1728` where the extra automorphisms merge
  kernels onto one j-invariant. The first version of this code treated the crater at `j = 0`
  as a floor; the regression test `test_descent_from_automorphism_craters` covers that.

## Fail-closed checks (every volcano, before any row is written)

1. Level sizes equal the class numbers `h(D_K (c ell^i)^2)` (PARI `qfbclassno`).
2. Every non-floor vertex has ell + 1 rational kernels, every floor vertex one.
3. The measured group structure of every stored curve equals `O_d/(pi - 1)` with
   `n1 = gcd(A, f/d)`, where `pi - 1 = A + f*omega` and d is the curve's conductor (Lenstra).
   The group is measured from point orders against the known factorisation of N, not read
   from the label.
4. Point counts (PARI SEA) equal `p + 1 - t` on a sample.
5. Blind descent from 24 random vertices per volcano reaches a floor curve at the correct
   level in the correct number of steps.

## Data (`data/`, gzip CSV)

* `volcanoes.csv.gz` one row per volcano: `ell, D_K, c, k, m, f, p, t, N`, factorisation of N,
  `largest_prime_bits`, per-level sizes, per-check results.
* `vertices.csv.gz` one row per curve (all curves for volcanoes up to 4,200, otherwise 300
  sampled per level; column `sampled`).
* `edges.csv.gz` every rational ell-isogeny between curves of fully stored volcanoes, with
  multiplicity and type `up`, `down` or `horizontal`.
* `descents.csv.gz` the blind-descent trials.
* `ecdlp_E1.csv.gz`, `ecdlp_E2.csv.gz` one row per solved discrete log.

Vertex columns, by role (this matters for machine learning):

* **Labels (ground truth, not observable from the curve):** `level`, `height = k - level`,
  `end_conductor`, `end_disc`, `end_class_number`, `comp` (ell-isogeny component).
* **Curve data:** `j, a, b` for `y^2 = x^3 + a x + b`; volcano-level `p, t, N` in `volcanoes`.
* **Cheap observables:** `n1, n2` (group structure), `aut_order`, `sig_2 ... sig_13`
  (number of rational q-isogenies, with multiplicity, for q in 2,3,5,7,11,13).
* **Graph-derived (leak the level near the floor):** `sig_ell` for the volcano's own ell (1 on the
  floor, ell + 1 elsewhere), `n_down`, `n_horizontal`, `ell_rational_mult`, `up_j`.
* `m` in `volcanoes` is the requested lower bound on `v_ell(A)`. The actual group structure is in
  `n1`. When `m > 0` the group structure alone reveals the level: `n1 = ell^min(v_ell(A), k - level)`.

## Reproduce

```
python3 -I test_volcano_descent.py
python3 -I build_dataset.py --out data              # about 4 minutes on 4 cores
python3 -I run_ecdlp.py --experiment E1 --trials 16 --per-level 2 --max-volcanoes 6
python3 -I run_ecdlp.py --experiment E2 --trials 20 --per-level 2
python3 -I analyze.py --write RESULT_TABLES.md
python3 -I ecc2k130_level1.py                       # about 1 minute
```

All randomness is seeded from a SHA-256 of the campaign seed and configuration, and PARI's
generator is seeded per volcano. Rebuilding gives identical `volcanoes.csv.gz`.

## Known limits

* Curves are CM-constructed. The crater is known by construction, so the dataset describes
  isogeny classes with a chosen large conductor, not random curves.
* Prime fields only. In characteristic 2 PARI's `polmodular` accepts only a prime-subfield
  evaluation point, so a binary-field volcano can be entered but not walked past the crater.
* The rho solver is pure Python and its counted unit (one affine group operation) is not the
  ecbench counted unit. Do not mix the two.
* Volcano depth is limited by `p >= |D_K| f^2 / 4`. The deepest here is 15 levels (ell = 2).
