# Isogeny walk from P-192

Status: **tooling and enumeration; no ECDLP cost measured**

Date: 2026-10-04

Result class: not an index-calculus measurement.  No IC/rho ratio and no
scoreboard row.  P-192 is registered in this change as
`icv1-fp192-t31607402316713927207482677199-52e4af59`
(`EC1P192Cp192h5531c4a08bdb`), and the leaderboard roster and lab browser
are regenerated.

This uses the walker, trait detection and S3 storage of PRs #1336 and #1355
unchanged, plus `isogeny_walk publish` (below).

## Class invariants

From `walk.json#class` and `class_audits.json`; they hold for every curve in
the class.

| | P-192 |
|:--|:--|
| trace `t` | 31607402316713927207482677199 |
| `t² − 4p`, trial factors `< 2^20` | `−5 · 11 · 31 · C`, `C` a 184-bit **probable prime** |
| `t² − 4p` squarefree and `≡ 1 (mod 4)` if `C` is prime | yes |
| twist order | `23 · c`, `c` a 188-bit composite, unfactored |
| embedding degree | `> 1000` |
| `ℓ ≤ 61` dividing `t² − 4p` | 5, 11, 31 (one horizontal edge each) |
| Elkies `ℓ ≤ 61` | 13, 23, 37, 43 |
| Atkin `ℓ ≤ 61` | 3, 7, 17, 19, 29, 41, 47, 53, 59, 61 |
| `ecc_safety` | 9 of 10 pass; **`order_size` fails** |
| PKM overall, special-prime score, Solinas weight | 0.647, 0.25, 3 |
| max twist leak bits | 5 |

**`order_size`.** The check fails because the 192-bit group order is below
`ecc_safety`'s 200-bit default, the SafeCurves recommendation.  This is
P-192's known security level of about 96 bits.  It is a property of the
whole class, not a finding about any walked curve.

**Endomorphism ring.** P-192 is the one root of the three walked so far
where trial division leaves a single probable-prime cofactor.  If `C` is
prime, then:
- `t² − 4p` is a fundamental discriminant;
- `Z[π]` is the maximal order;
- every curve in the class has `End(E) = Z[π]`, with discriminant `D_π`;
- every `ℓ`-volcano has depth 0, so the class is a single crater per `ℓ`.

The walk's root counts are consistent with that.  On all 11,700 expanded
curves there was 1 root for each of `ℓ = 5, 11, 31` and 2 for each Elkies
prime.  `C` passes Miller–Rabin only, so `curves.yaml` records
`endomorphism_discriminant` with status `unproved_in_this_registry` and
says so.

## Runs

| run | curves | expanded | edges | failures | orders proved | replay |
|:--|--:|--:|--:|--:|--:|:--|
| [`p192-ell61-radius2`](runs/p192-ell61-radius2/) | 71 | 12 | 132 | 0 | 71 | pass |
| [`p192-ell61-20k`](runs/p192-ell61-20k/) | 20,000 | 11,700 | 128,700 | 0 | 20,000 | pass |

**The radius-2 run** commits all of its outputs.

**The 20k run** is stored in S3 as
`s3://crypto-autoresearcher/isogeny-walk/runs/p192-71c04205135fcf8f`:
- the walk is 85 MB stored, 378 MB raw;
- `STORE.json` lists every object's key and hashes;
- it was fetched back with every hash checked, and `verify` passed
  (`VERIFY.json`, routes `a3df1ed1139e…`).

Its 8 trait shards and their merge are stored under
`runs/p192-71c04205135fcf8f/traits/8/`.  The merged traits SHA-256 is
`b240c5acd9f0…`, recorded in `collect.json`.

Curve detectors over the 20k curves:

| | count |
|:--|--:|
| `non_singular`, `generator_valid` | 20,000 / 20,000 each |
| `a_minus_3_model` | 10,012 (50.1%; expected ½ for `p ≡ 3 mod 4`) |
| smallest `b`, bits of `min(b, p − b)` | 177 (expected minimum over 20,000 draws ≈ 176.7) |
| `qr_prefix_64` range | 17–47 |

The class audits, re-run on 8 walked curves, gave the root's verdicts on
every check.

Practicality note: the walk took 206 s on a 14-core M4 Pro, unisolated.

## Offline, then upload

`isogeny_walk publish --dir DIR --store s3://…` uploads a walk or trait
shard made without `--store`, exactly as `--store` would have.  It is
refused if `curves.yaml` or `isogeny_routes.json` no longer match the hashes
in `walk.json`, and it finds an existing identical run instead of uploading
again.  Tested against S3 under `isogeny-walk/_e2e-1791164808/`; see
`docs/curves/ic/README.md#running-on-a-laptop-offline`.

## Not established

- Primality of `C`, and with it `End(E)` for the class.  That needs a
  primality certificate.
- Any DLP weakness or per-curve IC difference.
