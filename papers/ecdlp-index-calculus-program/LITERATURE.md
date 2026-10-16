# Literature notes for the survey rewrite

Working notes for §3–§5 of the survey. **Status of each entry: abstract or
landing page checked on the date given; full paper not yet read** unless
marked. The 2025–26 literature is already covered by
`research/literature_search_20261011/` (61 papers, abstract level) and is
not repeated here.

## Surveys

- Galbraith, Gaudry. *Recent progress on the elliptic curve discrete
  logarithm problem.* Designs, Codes and Cryptography 78(1):51–72, 2016;
  ePrint 2015/1022. Survey focused on summation-polynomial index calculus
  and the claims about its complexity. (checked 2026-10-15)

## Binary fields: the subexponential claims and the first-fall-degree debate

- Petit, Quisquater. *On polynomial systems arising from a Weil descent.*
  ASIACRYPT 2012, pp. 451–466; ePrint 2012/146. Conjectures degrees of
  regularity only slightly above the equations' degrees; heuristic
  O(2^(c·n^(2/3)·log n)), c < 2, for Diem's variant over F_{2^n}; estimates
  generic methods beaten only for n ≳ 2000. (checked 2026-10-15)
- Faugère, Perret, Petit, Renault. EUROCRYPT 2012 (binary-field index
  calculus exploiting the structure of the Weil-descended system).
  (landing page seen; bibliographic details to verify)
- Semaev. *New algorithm for the discrete logarithm problem on elliptic
  curves.* ePrint 2015/310; arXiv 1504.01175. Splits summation polynomials
  into cubic Boolean systems; under a first-fall-degree assumption the
  regularity degree is ≤ 4; heuristic 2^(c·√(n ln n)), c ≈ 1.69; claims
  several FIPS binary curves fall faster than Pollard, a conclusion the
  revision says does not hold under the weaker o(√(n/ln n)) assumption.
  (checked 2026-10-15)
- Karabina. *Point decomposition problem in binary elliptic curves.*
  ePrint 2015/319; arXiv 1504.02347. Auxiliary variables trade degree for
  variables. (checked 2026-10-15)
- Huang, Kosters, Yeo. *Last fall degree, HFE, and Weil descent attacks on
  ECDLP.* CRYPTO 2015; ePrint 2015/573. Introduces the last fall degree;
  computational evidence against the first-fall-degree assumption in
  characteristic 2; constructs summation-polynomial Weil descent systems
  where it is unlikely to hold. Follow-up arXiv 1505.02532.
  (checked 2026-10-15)
- Huang, Kosters, Petit, Yeo, Yun. *Quasi-subfield polynomials and the
  elliptic curve discrete logarithm problem.* J. Math. Cryptology, 2020.
  (checked 2026-10-15)
- Faugère, Huot, Joux, Renault, Vitse. *Symmetrized summation polynomials:
  using small order torsion points to speed up elliptic curve index
  calculus.* (EUROCRYPT 2014; venue to verify.)

## Decomposition without Gröbner bases (prior art for the root-finding line)

- McGuire, Mueller. *A new index calculus algorithm for the elliptic curve
  discrete logarithm problem and summation polynomial evaluation.* ePrint
  2017/1262. Index calculus with fast summation-polynomial evaluation that
  avoids Gröbner bases; a second algorithm avoids both Gröbner bases and
  summation polynomials; aimed at prime fields; slides (CAEN 2018) give
  O(p) under parameter conditions. Slides (CAEN 2018, read 2026-10-15):
  factor base F = {a_i P + b_i Q} (after Amadori–Pintore–Sala 2017);
  S_m evaluated by iterated resultants with S₃ in O(log² p); and "another
  variant": *choose m−1 points from the factor base and check if the sum of
  those points is in the factor base*, complexity O(p) for s ≥ (m−2)log₂ p.
  **This is the same mechanism as our two-summand root-finding oracle
  (subtract/add in the group, test membership), at least for prime
  fields.** The survey must present the S₃-root identity as a known
  consequence of the definition and the oracle as an engineering
  realisation with prior art (McGuire–Mueller; Amadori–Pintore–Sala),
  not as new. Full ePrint PDF blocked by the archive from this container
  (HTTP 403); read it before submission.

## Complexity analyses

- Yokoyama et al., J. Math. Cryptology 2020: linear-algebra analysis of the
  Gröbner computations in Semaev's naive index calculus under statistical
  assumptions. (abstract seen; authors to verify)
- Kudo, Yokoyama, ePrint 2024/528: maximal Gröbner degree via the
  saturation exponent; affine and homogenised runs agree up to the first
  degree fall. (from the 2026-10-11 search)

## Experimental studies and prime-field variants (checked 2026-10-15, abstracts)

- Shantz, Teske. *Solving the ECDLP using Semaev polynomials, Weil descent
  and Gröbner basis methods – an experimental study.* Number Theory and
  Cryptography (Buchmann Festschrift), LNCS 8260, pp. 94–107, 2013; ePrint
  2013/596. Experimental check of the Petit–Quisquater heuristics.
- Galbraith, Gebregiyorgis. *Summation polynomial algorithms for elliptic
  curves in characteristic two.* INDOCRYPT 2014; ePrint 2014/806. Binary
  Edwards variables lowering summation-polynomial degree; symmetry-breaking
  factor base and SAT vs Gröbner experiments (inconclusive); concludes
  Pollard rho remains much faster for prime-degree binary fields of
  reasonable size. (Our SAT and symmetry-breaking negatives agree.)
- Amadori, Pintore, Sala. *On the discrete logarithm problem for
  prime-field elliptic curves.* Finite Fields Appl. 51:168–182, 2018;
  ePrint 2017/609. Factor base built from the instance; one relation
  search among factor-base points; beats Semaev's and PKM's methods on
  prime fields.

## Generic lower bounds

- Corrigan-Gibbs, Kogan. *The discrete-logarithm problem with
  preprocessing.* EUROCRYPT 2018; ePrint 2017/1113. Generic algorithms
  with S bits of advice, online time T, success ε satisfy
  S·T² = Ω̃(ε·N), tight up to logs. Governs the retained-set-up (online)
  boundary of the frontier: the n=61 online result must be read against it
  with S measured in bits.

## Still to find and check

Diem's binary and extension-field theorems (2011 Compositio, 2013 J. Crypto
"On the discrete logarithm problem in elliptic curves II"); Joux–Vitse
(EUROCRYPT 2012/2013, cover and decomposition attacks over F_{p^6});
Gaudry–Hess–Smart and Menezes–Qu; Wiener–Zuccherato and
Gallant–Lambert–Vanstone; Trimoska et al. CP 2020 and the WDSat thesis.
