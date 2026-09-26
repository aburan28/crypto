# Preregistered finite-fibre S3 exactness proof and exhaustive controls

Status: frozen before any new small-field control outcome. The starting
observation is the independently replayed zero-spurious affine census in
[PR #774](https://github.com/aburan28/crypto/pull/774), which also found an
exact witness with a prefix at O. This PR tests and proves a narrower
mathematical statement; it does not export an ANF, run a solver, measure
relation yield, or infer ECC2K-130 attack speed.

The curve is `E: y²+xy=x³+1` over a finite field of characteristic two.
Write `L(a)={P in E(K):P affine,x(P)=a}`. A factor domain is **sign-complete**
when it contains all of `L(a)` for every admitted x. The formula frozen here
is `S3(a,b,c)=(ab+ac+bc)²+abc+1`. [Kosters–Yeo, Proposition 2.1 and
Remark 4.8](https://arxiv.org/pdf/1503.08001) provide general summation-
polynomial context and warn that higher-degree algebraic ideals need not
inherit every rational trace condition. The fibre-and-chain exactness claim
below is derived here and is not attributed to that paper.

## Hypothesis and falsification boundary

For rationally liftable `a,b in K`, the **distinct finite** roots `c in K`
of `S3(a,b,c)=0` equal
`{x(P+Q):P in L(a),Q in L(b),P+Q != O}`. In particular every such root
rationally lifts. For a sign-complete sequence of factor x domains, every
all-affine recursive-S3 path corresponds to one globally consistent signed
point tuple: if the next S3 step needs the opposite prefix sign, negate
all preceding factors, preserving every earlier prefix x. At an affine
terminal x, global negation gives either desired full target sign. The
converse holds for signed tuples with all prescribed prefixes and final sum
affine. A prefix or final sum at O has no x and needs an explicit branch.

The proof must cover four exhaustive local cases: `a=b=0` (S3=1, no finite
root, only the 2-torsion plus itself gives O); `a=b!=0` (one linear root
from doubling, opposite signs give O); `a!=b,ab=0` (one distinct double
root, because the x=0 point is self-inverse); and `a!=b,ab!=0` (two
distinct roots from P+Q and P−Q). It must distinguish the smooth S3
polynomial itself from any future Boolean exporter that clears a denominator.
The claim is falsified by a single rational liftable pair with a mismatched
root set or a single affine path lacking a signed full-target witness.

## Frozen exhaustive controls and negative controls

For degrees `n=1,2,3,4,5,6,7`, use exactly these bit-polynomial field
moduli: `0x3,0x7,0xb,0x13,0x25,0x43,0x83`. Verify their irreducibility
without trusting the table, enumerate **every** affine curve point and
every rational x fibre, check the group order from `t0=2,t1=-1,
t_n=-t_(n-1)-2t_(n-2)`, and for every ordered rational x pair evaluate
S3 at every field element. Compare the complete root set with every signed
point sum and count each of the four cases. Preserve the lexicographically
first `(n,a,b,c)` mismatch, if any; never drop a field or substitute one.

For complete chain controls enumerate all rational x tuples and all signed
point tuples for the exact panels `(n,m)=(3,4),(4,4),(5,3)`. For each x tuple,
enumerate all finite intermediate roots and terminal roots, then compare the
complete `(x tuple, intermediate x path, terminal x)` set with the group-law
set of signed tuples whose prefixes and final sum are affine. Check for
**every** candidate path that its realised full-point set equals the entire
rational fibre of its terminal x. Count the signed tuples with O prefixes
and those ending at O separately; do not call these an affine-chain error.

Two ordered negative controls document the assumptions. Scan `n=1..7`, then
increasing x, for the first non-rational x; use `a=b=x` and record its finite
S3 root despite `L(x)=empty`. Separately scan `n=1..7`, increasing `a<b`
among nonzero two-lift rational fibres, restrict each to its smaller-y point,
and retain the first S3 root absent from that restricted single-sign sum.
These controls show why unrestricted polynomial roots or sign-incomplete
factor sets can be spurious; neither is an accepted PDP factor domain.

## Independent replay and decision

The producer uses bit-serial field multiplication with Fermat inversion and
its own affine point law. The verifier uses the already independently checked
#762 polynomial-product/reduction and Euclid-inversion field/point law,
rebuilds all enumerations, and compares the full compressed pair/chain rows
and summary. Parent source, protocol, producer, verifier and runner bytes are
SHA-pinned in `FROZEN.json` before the first outcome. Run one producer and
one verifier as sequential cold children, each with a 180-second wall and
512-MiB peak-RSS acceptance cap. Retain raw rows, all stdout/stderr,
exit/UTC/CPU/wall/RSS and native field/curve counters, source/input hashes,
and any failure/censor receipt. A cap or replay mismatch leaves the theorem
control gate failed/censored; the mathematical proof must still stand on its
own and any counterexample must be reported.

Pass only if every pair and chain equality, both assumption controls,
field/order checks and independent replay pass within caps. A pass promotes
the *affine-chain semantic exactness theorem* for this curve and
sign-complete rational fibres, **not** a direct S6/S7 resultant claim,
denominator-cleared Boolean exporter, solver UNSAT, n131 runtime, or
end-to-end ECDLP gain. Update the canonical scoreboard with this scope and
the measured control counts, leaving full-DLP S/rho unset.
