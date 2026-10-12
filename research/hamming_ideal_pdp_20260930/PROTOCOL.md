# Hamming ideals and an ISD-like Gröbner oracle for point decomposition, 2026-09-30

Pre-registered protocol. Written before any cell was run. Deviations are
recorded in `RESULT.md`, never here.

## What is being tested

La Scala, Marchesin and Tiwari, *Hamming ideals and Gröbner bases for
ISD-like syndrome decoding* (2026), model the Hamming-weight constraint
`wt(v) = t` on `v ∈ F_2^n` as an ideal: by Lucas' theorem the binary digits
of `wt(v)` are the Boolean elementary symmetric functions `e_{2^k}(v)`, and
three lifted presentations with auxiliary variables (C-Hamming, FC-Hamming,
QFC-Hamming) keep the generators at bounded degree. They then decode with an
ISD-like strategy (`GBDecode`) in which only part of an information set is
fixed and the rest is solved by a truncated Gröbner computation inside a
branching search (`MultiSolve`, oracle `OracleT`). Their finding on Classic
McEliece Category 1 is negative: the Prange point (fix everything, linear
algebra only) is the cheapest configuration, because every coordinate left
free costs Gröbner calls that outweigh the combinatorial gain.

This repository's point decomposition problem (PDP) on the Koblitz family
`K_0 : y² + xy = x³ + 1` over `F_{2^n}` asks whether a target `R` is a sum of
`m` points whose abscissae lie in a factor base `F`. The decomposition note
[`RESEARCH_ECC2K130_DECOMPOSITION.md`](../notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)
§6 states two facts about prime `n`: a Frobenius-stable factor base of
subspace shape does not exist below dimension `n − 1`, and a Frobenius-stable
*set* (a union of orbits) "has no low-degree membership polynomial", so it
must be materialised (`2^l` stored points) and the algebraic oracle is lost.
The solver module `src/cryptanalysis/koblitz_groebner.rs` states the same
limitation from the other side: restricting summands to a subspace is "the
property the linearised-polynomial factor base has and a weight-bounded one
does not".

The paper contradicts the second half of both statements. In a **normal
basis** of `F_{2^n}`, Frobenius is a cyclic shift of coordinates, so the set
`F_w = { P : wt_NB(x(P)) ≤ w }` is Frobenius-stable at every `n`, needs no
storage (membership is a popcount), and by the paper has a membership ideal
with `O(n)` auxiliary variables and generators of degree at most
`⌊log₂ n⌋ + 1` (FC) or `2` (C, QFC). Whether that ideal makes a *useful*
oracle is the experiment.

## Hypotheses

- **H1 (the paper's construction transfers).** With the weight constraint
  encoded as a C-, FC- or QFC-Hamming ideal and the summation polynomial
  `S₃` Weil-descended in normal-basis coordinates, the `MultiSolve` oracle
  answers "is `R = P₁ + P₂`, `P_i ∈ F_w`" correctly on every target, and the
  number of `GroebnerSafe` calls per target grows more slowly than the
  exhaustive candidate count `|F_w|` as `n` grows.
- **H0 (the paper's negative result transfers).** Tame calls occur only
  when nearly all coordinates of one summand have been assigned, so the
  oracle is the exhaustive oracle with a Gröbner computation bolted onto
  each candidate, and the calls-per-target exponent in `|F_w|` is `≥ 1`.

The expected outcome is H0. The measurement decides.

## Frozen inputs

| item | value |
|:--|:--|
| curve | `K_0 : y² + xy = x³ + 1` over `F_{2^n}`, `a = 0, b = 1`, cofactor 4 |
| field degrees | `n ∈ {7, 11, 13, 17, 19}`; 7 is a smoke size, the four primes 11–19 are the fit |
| reduction polynomial | lowest-integer irreducible of degree `n` with constant term 1, the `find_irreducible` convention of `scripts/ecc2k130_point_decomposition.py` |
| normal basis | `{α^{2^i}}` for the smallest integer-encoded `α` whose conjugates are independent; recorded in every manifest |
| summands | `m = 2` (one `S₃` link; the quadratic case) |
| weight base | `F_w = {P : wt_NB(x(P)) ≤ w}` with `(n, w) = (7,2), (11,2), (13,2), (17,3), (19,3)`, chosen so the yield `C(|F_w|,2)/#E` lies near 1 |
| subspace base | random `F_2`-subspace `V` of dimension `l = ⌈log₂ #{x : wt(x) ≤ w}⌉`, seed-derived; `F_V = {P : x(P) ∈ V}`; the repository's oracle shape, run through the same solver |
| targets | 24 per `n`, seeds `0..23`; even seeds planted (`R = P₁ + P₂`, `P_i` uniform in `F_w`), odd seeds uniform in the odd-order subgroup; every target labelled by the exhaustive oracle for *each* base |
| encodings | `SUB` (subspace, polynomial-basis coordinates), `MONO` (control: all degree-`w+1` monomials in each summand's coordinates), `C`, `FC`, `QFC` (the paper's Theorems 3.5, 4.5 and Section 5) |
| digit constraints | `wt ≤ w` is the digit condition on the root ESF variables: `y_{2^k} = 0` for `2^k > w`, plus the algebraic normal form of `[Σ_{2^k ≤ w} t_k 2^k > w]` when `w + 1` is not a power of two (`w = 2`: `y₁y₂ = 0`) |
| monomial order | degrevlex, summand coordinates first (summand 1 then 2, decreasing index), auxiliary variables after, in tree order — the paper's `XRev-BitFold` |
| solver | `MultiSolve` with `OracleT`; `GroebnerSafe` = Boolean F4 truncated at pair degree `d = 6` with a budget `τ = 2³⁰` row-word XORs per call in place of the paper's 20-minute timeout; tame iff the basis contains `1` or reduces to linear polynomials with at most `2¹²` common zeros, every zero verified against the original system and lifted to points |
| branching | summand 1 coordinates in index order, then summand 2; `b = 0` before `b = 1`; a summand whose assigned ones reach `w` has its remaining coordinates forced to 0 (the paper's residual-weight bound); stop at the first verified decomposition |
| per-target cap | `2¹³` `GroebnerSafe` calls; exceeding it is a retained timeout |

## Reference and boundary

- **Reference:** the exhaustive oracle, the Prange analogue: for each
  `P₁ ∈ F` test `x(R − P₁) ∈ F`. Cost `|F|` group subtractions and `|F|`
  membership tests per target (a popcount for `F_w`, an `O(l)` subspace
  test for `F_V`). Both oracles may quotient by Frobenius; the factor `n`
  is common and is not applied to either.
- **Floor:** none new. The free-oracle floor of the decomposition note §5.3
  and the orbit-collapse floor of §6 are derived for any base of this shape
  and are unchanged by this experiment; the only quantity the weight base
  moves is the `2^l` storage term of §6, which becomes 0.

## Unit and table

The paper's unit is `GroebnerSafe` calls per target. That is the primary
column here, next to the exhaustive candidate count. Elimination row-word
XORs per call and per target are the second unit, reported without a
conversion to curve operations, because none is calibrated. Wall time is a
practicality note. One table, rows = `(n, base, encoding)`, columns:
variables, generators, maximum generator degree, targets, correctness
(oracle answer = exhaustive label on every target), mean tame calls, mean
wild calls, mean tame depth, mean XORs, exhaustive candidates `|F|`, and the
ratio `calls / |F|`. This is a **stage diagnostic**: no `S`, no rho ratio,
no end-to-end speedup is claimed or computed.

## Success and stop conditions

Success for H1 requires all of:

1. correctness on every target of every cell (zero disagreements with the
   exhaustive oracle; every "yes" carries a verified decomposition);
2. mean tame depth at most `n − 3` on the weight cells at `n ≥ 13`;
3. the exponent of mean `GroebnerSafe` calls per target in `|F_w|`, fitted
   over the four sizes `11, 13, 17, 19`, below `0.8` for at least one of
   `C`, `FC`, `QFC`.

Otherwise H0 stands. Stop conditions: a correctness disagreement stops the
run until the bug is found (it is a bug in this code, never a result); a
cell whose targets exceed the per-target cap on more than half of the
targets is marked `timeout` and the next larger `n` is not run for that
encoding, the fit then using the sizes that completed.

Inadmissible: changing `w`, `d`, `τ`, the target seeds, or the branching
order after seeing results; dropping timed-out targets from the means;
comparing the algebraic oracle's *first-solution* cost against the
exhaustive oracle's *full-scan* cost (the exhaustive scan also stops at the
first hit, and both are reported as full-tree and first-hit).

## Cost accounting

Charged per target, cold: system construction, every `GroebnerSafe` call,
solution verification, point lifting. Reported separately as common work:
curve enumeration, normal-basis search, factor-base enumeration, exhaustive
labelling. Nothing is amortised.

## Classification

By §3 of `AGENTS.md` the expected class is **accounting**: the decomposition
note's `2^l` storage term for a Frobenius-stable set is removed by a weight
base in a normal basis, and no operation count changes. An **advance** would
need success conditions 1–3 and, beyond this protocol, an end-to-end run.
