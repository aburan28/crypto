# Decision: admit higher arity before a native-leaf full-rank campaign

**PASS_BOUND; reject the m3 route at the frozen m10-scale physical-choice
budget.** The exact integer [result](RESULT.json) and [independent replay](VERIFY.json)
use the [preregistered protocol](PROTOCOL.md), source/input
[lock](FROZEN.json), and the challenge subgroup order
`q=680564733841876926932320129493409985129`. The protocol was committed
as `bfa61af45cab9c05e10c85117369e9fb3005fb22`, the final calculator
and verifier as `1b590c23b588523388c1761d905cd603db8ae20f`, and the
lock was pushed before the result. No n=131 PDP instance was solved.

The balanced m10 preflight has 7,977 physical points per slot. Keeping its
total of 79,770 **slot choices** but using three ordered m3 slots gives the
most favorable allocation `[26,590,26,590,26,590]` and at most
18,799,877,179,000 labelled triples. On a uniformly marginal target in one
q-size subgroup coset, the support probability is at most the exact fraction
`18799877179000/q`, approximately `2.7624e-26`. Even one million such
probes, with arbitrary correlations, have probability at most
`18799877179000000000/q ≈ 2.7624e-20` of *any* hit. A single-row-per-probe
relation campaign at this budget therefore has no plausible path to rank.
The parent `admission.py` independently reproduced the exact hit and
rank-one upper fractions. These are bounds, not measurements of where any
triple sums actually land.

| Desired **necessary** support ceiling | Minimum combined physical slot choices | Balanced slot sizes at that minimum | Conditional minimum nonzero projected signed-orbit log classes |
| ---: | ---: | --- | ---: |
| 1% | 5,685,182,383,138 | `[1,895,060,794,380, 1,895,060,794,379, 1,895,060,794,379]` | 5,424,792,351 |
| 50% | 20,944,390,974,995 | `[6,981,463,658,332, 6,981,463,658,332, 6,981,463,658,331]` | 19,985,105,893 |

The minimum budgets are certified by the exact tuple-product values at
`S-1` and `S` in `RESULT.json`. The log-class floor assumes every chosen
nonzero projected factor point remains a required column, allows up to four
rational preimages under `[4]`, and gives the leaf a *free* transported
signed Frobenius orbit of length up to 262. It is deliberately favorable to
m3; actual leaf action and matrix costs are unmeasured. The implicit
`[44,44,43]` domains in PR #976 have a separate combined physical **upper**
bound of 87,960,930,222,077 points, not an observed factor count or a
rank-ready base. Passing a representation-size cap does not resolve this
support-versus-rank tradeoff.

This is a quantitative no-go **only** for the frozen m3, 79,770-choice,
one-million-uniform-target, one-row-per-probe route. It does not exclude a
larger or differently compressed m3 base, a better relation generator,
multiple independent witnesses per probe, or a higher-arity solver. The
evidence-ranked next attack experiment is the existing review-gated m10
representation-capacity gate [PR #937](https://github.com/aburan28/crypto/pull/937):
its balanced preflight has a tuple ceiling above q with 3,988 signed columns,
while actual natural PDP yield, solved CNFs, rank, and all charged costs remain
unknown. PR #937's independent-review and dispatch conditions must be met as
written. After representation admission, freeze natural-target solver/yield
and physical-support/rank panels; only successful panels justify paired
original/transported/native/pullback full-log costs against matched rho.
Keep ECC2K-130 speedup, n=131 solver runtime and rho crossover unset.
