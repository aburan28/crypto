# Isogeny-class trio: three 61-bit roots for a whole-class experiment

Date: 2026-10-17. Protocol: [`PROTOCOL.md`](PROTOCOL.md), frozen before the B
and P searches; its 2026-10-17 amendment records the user's choice of the
61-bit trio as the primary picks. PDF: [`REPORT.pdf`](REPORT.pdf). Data:
[`data/`](data/). Tool: [`selector/`](selector/) (native Rust, no
dependencies) and [`identity/`](identity/) (the repository's native ICV1
producer).

## Question

The user asked for three curves, one binary, one Koblitz and one prime-field,
between 41 and about 66 bits (preferably above 45), with shapes that
emulate standardized curves. Each should have Frobenius conductor `f > 1`, so
that its isogeny class has volcanoes with more than one level, as the root of
an experiment that iterates over every curve in the isogeny class. After
seeing both candidate trios the user chose the 61-bit one.

## Picks (61-bit trio, chosen by the user)

| | Koblitz (K-shape) | Binary (B-shape) | Prime (P-shape) |
|:--|:--|:--|:--|
| ICV1 slug | `icv1-f2m61-t158598901-ab42b6c5` | `icv1-f2m61-tm2657393969-b9ff0fb0` | `icv1-fp61-tm813574577-a3b71255` |
| Emulates | K-233/283/409/571 (`a = 0`); cofactor is not standard (none exists, finding 1) | B-163…B-571 (`a = 1`, seeded `b`, cofactor 2) | P-521 (`p = 2^k − 1` Mersenne, `a = −3`, seeded `b`, prime order) |
| Field | `GF(2^61)`, `x^61 + x^5 + x^2 + x + 1` | same | `GF(p)`, `p = 2^61 − 1` |
| Model | `y^2 + xy = x^3 + 1` | `y^2 + xy = x^3 + x^2 + b`, `b = 0x9c60f836d06fc8` | `y^2 = x^3 − 3x + b`, `b = 0x54aa66b458b77ec` |
| Seed | — (fixed by the family) | `isogeny-class-trio/2026-10-17/B-61/b/12079` | `isogeny-class-trio/2026-10-17/P-61/b/11526` |
| `#E` | `4 · 3539 · 162888033982417` | `2 · 1152921505935543961` | `2305843010027268529` (prime) |
| Prime subgroup | 48 bits (cofactor 14156) | 61 bits | 62 bits |
| Trace `t` | 158598901 | −2657393969 | −813574577 |
| `d_K` | −7 | −900303761090047 | −6764617042579 |
| Conductor `f` | `1951 · 587551` | `7^2` | `3^2 · 5^3` |
| Volcanoes | `l = 1951`, `l = 587551`: height 1, inert | `l = 7`: height 2, split crater | `l = 3`: height 2, inert; `l = 5`: height 3, split |
| Strata (`c \| f`) | 4 | 3 | 12 |
| Class size `H` | **1,147,491,009** | **851,960,648** | **920,783,750** |
| Root's position | crater (`End = Z[τ] = O_K`) | floor (`End = Z[π]`, conductor 49) | floor of both volcanoes (`End = Z[π]`, conductor 1125) |
| Generator | `(0x5cc1d6ba365013e, 0x171d9a59a2d31eb9)` (registry, `EC1N61Ce0hd8c15448081e`) | `(0xc17a1d9ae3996ee, 0x9bcf1bc6ebcb05)` | `(0x112c917afe3d4fd3, 0xef721929b08bfaa)` |

Stratum sizes (curves with `End = O_c`, exact):

- K-61: `c = 1`: 1; `1951`: 1,952; `587551`: 587,552; `1951·587551`: 1,146,901,504.
- B-61: `c = 1`: 17,386,952; `7`: 104,321,712; `7^2`: 730,251,984.
- P-61: `c = 1`: 433,310; `3`: 1,733,240; `5`: 1,733,240; `9`: 5,199,720;
  `15`: 6,932,960; `25`: 8,666,200; `45`: 20,798,880; `75`: 34,664,800;
  `125`: 43,331,000; `225`: 103,994,400; `375`: 173,324,000;
  `1125`: 519,972,000.

**K-61 caveat for the whole-class pass.** The Koblitz curve is the only
curve of its class on the crater, and every small-`l` isogeny from it is a
loop (finding 2). The `O_1951` stratum (1,952 curves) is one degree-1951
step away, with kernel x-coordinates in `F_{q^975}`. The other 1,147,489,056
curves (99.9998% of the class) need a degree-587551 step whose kernel
x-coordinates lie in `F_{q^58755}`, a field of about 3.6 million bits. That
step is the open obstruction for iterating over this class; B-61 and P-61
have no such step (degrees 7, 3 and 5).

![Strata of each class](figures/strata.svg)

![Volcano tiers](figures/volcanoes.svg)

### 41/47-bit trio (the first proposal, kept for reference)

| | K-41 | B-47 | P-47 |
|:--|:--|:--|:--|
| ICV1 slug | `icv1-f2m41-tm2308219-7f48b14a` | `icv1-f2m47-t4276051-826c5755` | `icv1-fp47-tm8785349-43e36ce5` |
| Field | `x^41 + x^3 + 1` | `x^47 + x^5 + 1` | `p = 2^47 − 115` |
| Model | `a = 0`, `b = 1` (cofactor 4, the only standard-cofactor K in range) | `a = 1`, `b = 0x52e3db551eb3` (seed index 28948) | `a = −3`, `b = 0xda8a554f2ba` (seed index 117828) |
| Subgroup | 40 bits | 46 bits | 48 bits |
| `f` | `409 · 1721` | `31^2` (inert crater) | `3^5` (split crater) |
| Strata, `H` | 4, 708,153 | 3, 11,754,700 | 6, 13,461,228 |
| Root | crater | floor | floor |

Their stratum sizes: K-41 1 / 410 / 1,722 / 706,020; B-47 11,468 /
366,976 / 11,376,256; P-47 55,396 / 110,792 / 332,376 / 997,128 /
2,991,384 / 8,974,152.

## Findings

1. **No prime-degree Koblitz curve has the standard cofactor between 2^42
   and 2^66** (derived, exhaustive). The prime degrees in range are 41, 43,
   47, 53, 59 and 61; of the 12 curves only `K_0` over `GF(2^41)` has
   `#E = 4r` with `r` prime ([`data/koblitz.jsonl`](data/koblitz.jsonl),
   [figure](figures/koblitz-table.svg)). Above 45 bits the largest
   prime-degree subgroup is `K_0` over `GF(2^61)` (cofactor 14156, 48-bit
   subgroup), the K pick.
2. **A prime-degree Koblitz class has height-1 volcanoes at large `l`
   only** (derived; classical). `f = |U_m|`, the Lucas sequence of
   `τ^2 − μτ + 2`. For prime `m ≠ 7` every prime `l | U_m` has rank of
   apparition `m`, so `l ≡ ±1 (mod m)` and `l ≥ 2m − 1`; the table bears
   this out (1951, 587551 ≡ −1, −1 mod 61; 409, 1721 ≡ −1, −1 mod 41). The
   Koblitz curve is the only curve in its class with `End = O_K`
   (`h(−7) = 1`), so for every small `l` its isogeny graph is one vertex
   with loops: **the rest of the class is reachable only through the
   vertical steps of degree `l | f`.** For K-61 those are degree 1951 and
   587551, with kernel x-coordinates in `F_{q^975}` and `F_{q^58755}`
   (computed: the order of `±t/2` in `F_l^*/±1`).
3. **Binary curves over odd-degree fields have conductor primes
   `l ≡ ±1 (mod 8)` only** (derived; confirmed by the searches). `l | f`
   needs `t^2 ≡ 2^(m+2) (mod l)`, i.e. 2 a square mod `l`. So 3- and
   5-volcanoes never occur for B- or K-shape curves; the smallest possible
   vertical degree is 7. Every conductor prime in the B-61 and B-47 searches
   was in {7, 17, 23, 31, 41, 47, 71, 73, 79, …}
   ([`data/tally-B61.json`](data/tally-B61.json),
   [`data/tally-B47.json`](data/tally-B47.json)). For the same reason only
   2 of 1,053 standard-cofactor B-61 curves have a height-2 volcano,
   against 24 qualifying of 473 for P-61.
4. **Both seeded roots sit on the floor** (proved by the tests below), as
   most curves do: the floor holds 85.7% of B-61's class and 56.5% of
   P-61's. Their crater representatives are 2 (B-61) and 2 + 3 (P-61)
   ascending isogenies away.

![Search outcomes](figures/search-depth.svg)

## Method

- **Koblitz.** `#E(GF(2^m)) = 2^m + 1 − V_m` and `f = |U_m|` from the Lucas
  sequences of `τ`. The selector asserts `t^2 − 2^(m+2) = −7 U_m^2` for every
  row. The orders of K-41 and K-61 match their registry entries
  (`icv1-f2m41-tm2308219-7f48b14a`, `icv1-f2m61-t158598901-ab42b6c5`), and the native ICV1 producer reproduces their slugs.
- **B and P searches.** Coefficients come from SHA-256 of public seed
  strings, as in X9.62/FIPS 186 seeding. The selector counts points by BSGS
  over the Hasse window on seeded points, using `#E ≡ 2 (mod 4)` for B.
  A curve is kept only when the match is unique and `r = #E/h` is prime.
  `r > 4√q` then makes `#E` the only multiple of `r` in the window, so each
  kept order is proved. Not one count was ambiguous across the 327,680
  curves.
- **Selection rule (frozen).** `f > 1`, `f` 61-smooth, a volcano of height
  ≥ 2. Rank by height, then number of strata, then lowest seed index.
- **Class invariants.** `Δ = t^2 − 4q` is factored to `f^2 d_K`.
  `h(d_K)` is computed exactly by counting reduced forms (`O(√|d_K|)`;
  cross-checked against brute force for every fundamental `|D| < 4000`).
  Each stratum's `h(c^2 d_K)` comes from the conductor formula, and
  `H = Σ_{c|f} h(c^2 d_K)` (Schoof 1987, Thm 4.6).
- **Root levels.**
  - **K:** `τ ∈ End(E)` gives `End = O_K` (crater).
  - **B:** for odd `l | f`, above the floor `π` acts on `E[l]` as the
    scalar `t/2`, so all `(l^2 − 1)/2` x-coordinates of `E[l]` lie in
    `F_{q^k}`, with `k` the order of `±t/2` in `F_l^*/±1`. On the floor
    `π − t/2` is a non-zero nilpotent, so only the eigenline's `(l − 1)/2`
    x-coordinates do. The test is therefore
    `deg gcd(ψ_l, x^(q^k) − x)` using characteristic-2 division
    polynomials. B-47 gives 15 of 480 (floor). Control: over `GF(2^7)`
    the Koblitz class has `f = 7`, and the test flags exactly the Koblitz
    curve as above the floor and the other 7 class members as floor
    (`floor_test_control_koblitz_m7`).
  - **P:** a 3-isogeny BFS (roots of `ψ_3`, Vélu) finds that the root has
    exactly one rational 3-isogeny. The repository walker independently
    agrees: `src/bin/isogeny_walk.rs`, with the same trace, `Δ`, a 3-volcano
    of depth 5, and one kernel-certified 3-edge from the root
    ([`data/walker-class-P47.txt`](data/walker-class-P47.txt)). For P-61
    the walker records one 3-edge and one 5-edge from the root, so the root
    is on both floors.
- **Tests.** `cargo test --release` in `selector/` runs 15 tests:
  - BSGS against naive counts, in both characteristics;
  - characteristic-2 3-isogenies against brute-force orders, plus their
    duals;
  - division polynomials against rational torsion;
  - class numbers against brute force;
  - SHA-256 vectors;
  - the floor-test control.

## Status of each statement

| Statement | Status |
|:--|:--|
| Orders, traces, `f`, `d_K` of all six curves | computed (native); order proved by `r > 4√q` (B, P) or exact formula (K) |
| `h(d_K)` and stratum sizes, `H` | computed exactly (native form count); formula classical |
| Root levels | B, P: proved by the tests above; P cross-checked by the walker; K: theorem |
| No standard-cofactor prime-degree Koblitz in 42–66 | exhaustive enumeration |
| `l ≡ ±1 (mod m)` for Koblitz, `±1 (mod 8)` for binary | classical arguments, written above; consistent with all data |
| Seeded picks are optimal under the rule | within the frozen seed ranges only |

No theorem in this report is new, so there is no Lean obligation (AGENTS.md,
"Mathematical results should include Lean proofs"). Every structural fact
used is classical and its argument is written out above.

## What the whole-class experiment will need (not done here)

- **P-61** (921 M curves, 12 strata): the repository walker supports
  `p = 2^61 − 1` and every odd `l`, and the vertical steps have degree 3
  and 5. It holds about 87 KB per curve, so a full pass needs a compact
  class-group enumeration such as the `p61_class` segment stream
  (`research/p61_isogeny_class_20261014/`, same field, another class).
- **B-61** (852 M curves): there is no binary isogeny walker
  (`docs/isogeny-walk/README.md`: "Prime fields and odd l only").
  - It needs characteristic-2 `l`-isogenies: degree 7 for the vertical
    steps (kernel x-coordinates in `F_{q^3}`), and split primes
    (`(d_K/l) = 1`) for horizontal steps.
  - This selector already has `ψ_l` and 3-isogeny Vélu in characteristic 2.
  - 2-isogenies are horizontal here (2 splits, since `d_K ≡ 1 mod 8`).
- **K-61** (1.15 G curves): a degree-1951 step (kernel field `F_{q^975}`)
  reaches 1,952 curves. The remaining 1,147,489,056 need a degree-587551
  step with kernel field `F_{q^58755}`, which no tool here can compute yet.
  That is the open obstruction for this class.

## Reproduce

```sh
cd research/isogeny_class_trio_20261017/selector
cargo test --release && cargo build --release
B=target/release/isogeny-class-trio-selector
$B koblitz --min 41 --max 64                       > ../data/koblitz.jsonl
$B binary --m 47 --count 131072 --threads 4        > ../data/search-B47.jsonl 2> ../data/search-B47.summary.json
$B prime --bits 47 --count 131072 --threads 4      > ../data/search-P47.jsonl 2> ../data/search-P47.summary.json
$B binary --m 61 --count 32768 --threads 4         > ../data/search-B61.jsonl 2> ../data/search-B61.summary.json
$B prime --bits 61 --count 32768 --threads 4       > ../data/search-P61.jsonl 2> ../data/search-P61.summary.json
$B analyze-binary --m 61 --a2 0 --b 1 --koblitz    > ../data/pick-K61.json
$B analyze-binary --m 61 --a2 1 --b 0x9c60f836d06fc8 --seed isogeny-class-trio/2026-10-17/B-61/b/12079 > ../data/pick-B61.json
$B analyze-prime --p 2305843009213693951 --b 0x54aa66b458b77ec --seed isogeny-class-trio/2026-10-17/P-61/b/11526 > ../data/pick-P61.json
$B check-point --m 61 --a2 0 --b 1 --r 162888033982417 --gx 0x5cc1d6ba365013e --gy 0x171d9a59a2d31eb9
$B tally --in ../data/search-B61.jsonl             > ../data/tally-B61.json
# (41/47-bit trio: same commands with --m 41 / 47, --p 140737488355213 and the seeds in the table)
```

Wall time on the 4-core session host: B-47 search 1 m 59 s, P-47 search 4 m 18 s,
B-61 and P-61 searches about 27 min together (sharing the host with a build),
analyses under 5 s each. Figures: `python3 figures/plot.py` (plots recorded
values only). The curves are registered in `docs/curves/registry.json` from
[`curves.json`](curves.json).
