# Extension-field torsion-coordinate DLP comparison: frozen protocol

Status: protocol frozen before instance screening or paired solves. This is a
follow-on to [the cubic coordinate certificate](../torsion_orbit_20261013/REPORT.md)
and [its next gate](../torsion_orbit_20261013/NEXT_GATE.md). No observation in
this file is a measured advantage.

## Question and boundaries

Does the nonfactorable coordinate `chi=(y+1)/(y-1)` and its public
order-three action `f(z)=1/(1-z)` reduce the **complete cost of a single
ECDLP** on `E:y^2=x^3+1` over a proper extension field? The null is that
torsion-orbit bookkeeping and the degree-three isogeny leave relation yield
unchanged or add cost. The generic two-class count floor, fixed before tuning,
is `S_floor=sqrt(pi/4)=0.886226925`, where `S` is all counted work divided by
`sqrt(r)`. The reference is same-target negation rho; include a stronger
automorphism rho where its action stays inside the selected subgroup. An
extension-field multiplicative DLP route must be screened by `ord_r(p)` and
priced separately if a usable pairing reduction is constructed. Do not call a
phase ratio an end-to-end ratio.

## Instance selection, frozen before the paired run

Use the existing native `Fp2` arithmetic: `F_(p^2)=F_p[s]/(s^2-nu)`, where
`nu` is the smallest quadratic nonresidue selected by `Fp2::new(p)`.
Screen primes `31 <= p <= 4093`, `p=1 mod 3`, in ascending order. Count
`n1=#E(F_p)` exactly by the Legendre symbol, set `t=p+1-n1`, and compute
`n2=(p+1)^2-t^2`. Let `r` be the largest prime divisor of `p+1+t`.
Accept only if `r>=11`, `gcd(r,3)=1`, `r` does not divide `n1`, and
`gcd(r,n2/r)=1`. Verify a nonidentity generator `G=[n2/r]P`, `[r]G=O`,
`G` has a non-base-field coordinate, and that the explicit degree-three map
`phi:E -> E':y^2=x^3-27` preserves the group law on the checked samples.
Record `ord_r(p)` and the base-field isomorphism test for `E,E'`.

Selector amendment v2, dated 2026-10-14 **after the v1 screen and before any
paired DLP run**: the v1 selector returned subgroup orders `19,31,19,73` as
`p` increased. Retain [SCREEN_V1.jsonl](SCREEN_V1.jsonl) and require the
minimum subgroup orders `17,53,211,811` in the four respective intervals.
This changes the instance screen to produce actually increasing DLP sizes;
it changes no solver or cost threshold. Choose the first accepted instance
in each of the fixed prime intervals `[31,127]`, `[128,511]`,
`[512,2047]`, `[2048,4093]` meeting those minima and giving at
least eight distinct nonidentity projected base elements from some line
`x=c*t`, `t in F_p`, where `c` is the first of `s,1+s,2+s,...` meeting this
condition. Enumerate signed lifts in increasing `t`, then packed y order.
For the first run use eight projected columns; later runs use the same rule
with `min(32,ceil(sqrt(r)))` columns. Deduplicate by the subgroup projection
`[n2/r]P`, accounting for every tested candidate, including rejected ones.
Freeze the resulting exact `(p,nu,r,G,c,base,targets,seeds,caps)` as JSON and
commit it before any paired solve. Selection is a search, not a cost result.

The quadratic extension is deliberately **ordinary** (`p=1 mod 3`) to avoid
the supersingular `p=2 mod 3` source from the prime-field screen. Do not
substitute an anti-trace x-line whose target and entire prime subgroup also
have `x in F_p`: then the descended Semaev system loses an equation and its
apparently cheap oracle is not the intended extension-field test.

## Paired solves

Use the same public source target `Q`, point order `r`, factor-base
projection, deterministic random stream, trial cap, and verified scalar in
each arm. Solver code receives only `G,Q,r,seed`, never the known scalar.

1. Direct source `S3(x1,x2,xR)` line resultant, independently selected base.
2. Source torsion-orbit canonicalisation with corrections checked in the
   full group, then projected through the cofactor. Charge all translated
   points and corrections.
3. Nonfactorable `chi` polynomial/decomposition oracle with the order-three
   Mobius action; certify returned points by exact group addition.
4. Invariant `x(phi(P))` quotient control; price map construction/evaluation,
   destination relation collection and transport/recovery. Include an
   isomorphism-only destination-coordinate control if an isomorphism exists.
5. Same-target matched negation rho and rho A/A control. Include stronger
   automorphism rho when applicable.

For each oracle, first exhaustively compare its yes/no and returned sums with
an independent finite pair enumeration on the smallest instance. Report
distinct projected base size, signed points, equal-memory and equal-call
comparisons separately. Recover factor-base logarithms from verified rows,
then the target logarithm; independently replay `[k]G=Q` in a second native
process. Do not accept a collector hit or planted answer as a solve.

## Cost and stop rules

Record setup, field multiplication/squaring/inversion, roots, polynomial
construction, resultant or SAT/F4 work, table lookups, map/translation,
relation checks, unsuccessful queries, rank/linear algebra, target descent,
recovery verification, memory peak, and per-phase CPU/wall time. Count native
work that lacks a conversion as explicitly unpriced. Publish the cold
one-target total, plus reusable setup and the exact reuse count for any
amortized column. One table reports each arm's total in the **same** unit,
`S`, ratio to the floor and matched rho, and correctness. If a primitive
cannot be converted reproducibly to that unit, leave `S` incomplete and
present a separate whole-process wall/CPU comparison with its isolation
grade. No portable wall claim from a shared-host L0 run.

Caps: at most `10^6` relation queries and `10^7` polynomial/root operations
per arm per instance; keep every exhausted or failed run. Stop a row after a
verified answer or the cap. Do not infer absence of a mathematical route from
a cap. Before scaling claims require all four sizes, isolation L2 or better,
and independent-host replay. The first small complete solve is a diagnostic
against the same-target extension-field ECDLP; it does not discharge the
four-size scaling gate.

Native Rust only. Source revisions, input hashes, exact commands, raw JSON,
replay checks, a Markdown report, SVG and PDF land in this study; timing runs
also copy the report and graph into `benchmarks/` and update its index. Every
new algebraic identity needs an `IDC1h...` certificate and a Lean proof map.
