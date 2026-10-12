# Binary summation formula audit

Run `python3 audit.py --output results.json` with Python 3.10 or newer.
Standard library only. Fields are deliberately fixed to GF(8), GF(16), GF(32).

## Question

Does varying the constant coefficient in the binary three-point summation
polynomial simplify its Boolean equations, when the coordinate representation
and subspace are held fixed?

For `f_c(x,y,z) = (xy+xz+yz)^2 + xyz + c`, we have `f_c = f_0 + c`.
After substituting binary coordinates for x, y, z in a fixed linear subspace,
the Boolean algebraic normal forms therefore have identical nonconstant terms.
The experiment verifies this exact identity and audits its consequences.

## Results

798 systems across 34 coordinate subspaces were exhaustively evaluated.
All 798 constant-only comparisons passed. All 219,968 pointwise Frobenius
covariance checks and the same number of permutation checks passed.
The maximum Boolean algebraic degree was 3 for every sampled subspace.
Zero counts varied with c in all 34 subspaces.

| Field | Subspace dimension | Systems | Polynomial zeros per system, range |
|---|---:|---:|---:|
| GF(8) | 2 | 21 | 0–15 |
| GF(8) | 3 | 7 | 57–81 |
| GF(16) | 2 | 90 | 0–15 |
| GF(16) | 3 | 60 | 10–55 |
| GF(32) | 2 | 310 | 0–12 |
| GF(32) | 3 | 310 | 0–58 |

These are counts of ordered x-coordinate triples, including repetitions and
zero coordinates. They are NOT counts of verified rational-point relations.
The c=0 polynomial is used only as an algebraic reference; c=0 does not define
a nonsingular curve in this model. All 798 evaluated curve coefficients are
nonzero.

## Interpretation

Coefficient selection changes the right-hand side of a fixed nonlinear map.
It can change the solution count without changing any nonlinear coefficient.
This rules out an explanation based on coefficient-induced cancellation of
the original nonconstant terms for this fixed model and fixed subspace.
It does not rule out changes in elimination behavior or consequences of
changing the subspace itself.

Frobenius transport satisfies
`f_(c^2)(x^2,y^2,z^2) = f_c(x,y,z)^2`.
Squaring is bijective, so transporting both c and the subspace preserves the
zero count. Comparing conjugate coefficients against an untransported subspace
would compare different constraints, not disprove that equivalence.

## Limits and next scientific gate

This is an algebra audit, not an ECC2K-130 attack or performance improvement.
It does not construct isogenies, classify endomorphism rings, hold the point
count fixed, implement F4/F5, or measure degree of regularity. Coefficient
sweeps mix isogeny classes. No conclusion about volcano levels follows.

The next bounded experiment should compare small-field curves in a verified
common isogeny class, recording the exact subspace and equation encoding.
Record polynomial zero counts separately from rational lifts. For any
elimination measurements, distinguish input Boolean degree, first fall degree,
maximum processed solver degree, and a formally justified degree of regularity.
Keep unsatisfiable instances separate: proving inconsistency quickly is not
evidence of faster production of a solution. Use transported Frobenius pairs
as equivalence controls before interpreting cross-curve differences.

## Sources

- Semaev, Summation polynomials and the discrete logarithm problem on elliptic
  curves: https://eprint.iacr.org/2004/031
- Some relations between Semaev's summation polynomials:
  https://api.lib.kyushu-u.ac.jp/opac_download_md/19584/JMI2011A-9.pdf
- Kosters and Yeo, Notes on summation polynomials:
  https://arxiv.org/abs/1503.08001

`results.json` contains all field moduli, subspace bases, coefficients,
truth-table sizes, zero counts, coordinate degrees, and SHA-256 fingerprints
of nonconstant field-valued ANF coefficients. Fingerprints are comparable
within the same fixed coordinate encoding only. The script verifies the
ANF transform by its involution property and checks multiplicative inverses
for every nonzero element in each tiny field.
