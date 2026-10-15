# Small-x factor base on the SECG curves: exact two-point relation yield and the cost of a trial

Run `2026-10-14`, seed 7, release `ic`. Report: [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)).

## Answer

For a small-x base of `B` abscissae (`2B` points) on a SECG curve, the probability that a
random element of the prime-order subgroup is a signed sum of two base points is exactly
`2B²/(n − 1)` (`n` the subgroup order): no pair sums coincide on 27 of the 33 curves;
`sect163k1`, `sect233k1`, `sect239k1`, `sect283k1`, `sect409k1`, `sect571k1` have 12 to 22
collisions among 2²⁴ sums, from Frobenius images in the base (`x`, `x²`, `x⁴` are all small
abscissae and `τ² = μτ − 2`), a relative effect of order 10⁻⁶. Cofactor `h > 1` is handled by
taking `[h]P` of each small-x point, so the base and the targets lie in the subgroup.

| curve (B = 4096) | `log₂ n` | exact yield | trials for `B + 1` relations | rho (`log₂ √n`) | S₃ per target | `B + 1` relations by S₃ vs rho |
|:--|--:|--:|--:|--:|--:|--:|
| `sect113r1`, `sect113r2` | 112.0 | 2⁻⁸⁷·⁰ | 2⁹⁹·⁰ | 2⁵⁶·⁰ | 5.4 to 7.1 ms | 52 to 53 bits slower |
| `sect163r1`, `sect163r2` | 162.0 | 2⁻¹³⁷·⁰ | 2¹⁴⁹·⁰ | 2⁸¹·⁰ | 8.7 to 9.6 ms | 77 bits slower |
| `secp160k1`, `secp160r1`, `secp160r2` | 160.0 | 2⁻¹³⁵·⁰ | 2¹⁴⁷·⁰ | 2⁸⁰·⁰ | 42 to 46 ms | 79 bits slower |

All 33 curves (`secp112r1` to `secp521r1`, `sect113r1` to `sect571r1`) are in `REPORT.pdf`; B = 8192 on
30 of them (the other three, `secp521r1` and `sect571k1/r1`, are refused by the pair-sum memory guard).
Two-point small-x sums lose to rho on every curve and run, by at least 2⁴⁰·⁹ in trials, and
deciding a trial by solving Semaev's S₃ (the cost the yield leaves out) costs 2600 (`sect113r1`) to 313 000 (`secp224r1`) times
the table lookup per target. Monte Carlo finds 0 hits in 91.6 million trials, as the yields
require. On 52 toy curves (binary and prime, `h` = 2 and 4) the Monte Carlo rate matches the
exact count (standardized differences: mean 0.02, standard deviation 1.04) and the S₃ decision
equals the pair-sum table on every field element of toy fields. The raw small-x base (no `[h]`)
keeps 1/h of the sums on nine of ten curves sampled; `sect233k1` is 3.8 standard errors above.

"61 or 63, 160 bit" was read as `sect163r1/r2` (163-bit) and, as the alternative reading,
`secp160k1/r1/r2`; both are measured.

## Use

```sh
ic list                                       # every SECG curve, plus ecc2k-95/130
ic inspect --curve sect113r1                  # (also: ic sect113r1)
ic prime --curve sect113r1 --size 4096        # binary, cofactor 2
ic prime --curve sect163r2 --size 4096 --s3-trials 20   # plus the S3 cost of a trial
ic prime --curve secp160r1 --size 4096        # prime, 160 bits
ic prime --curve sect233k1 --size 64 --raw-samples 4000 # the raw small-x base
ic prime --curve toy-binary:17,7,4 --size 128 # toy binary curve, cofactor 4
```

Code: `src/cryptanalysis/secg_standards.rs` (typed loader over
`docs/curves/standards/parameters.json`), `secg_small_x.rs` (entry point),
`binary_small_x.rs` (`GF(2^m)` with 2 to 9 limbs, binary curve, small-x base, exact pair-sum count,
toy curves), `prime_small_x.rs` (generic Montgomery prime field with 3 to 9 limbs, Tonelli–Shanks,
cofactor, toy curves), `semaev_s3_cost.rs` (S₃ decision and its timing), `src/bin/ic/prime.rs`,
`src/bin/ic/params.rs`, `examples/secg_koblitz_collisions.rs`.

## Files

| path | what it is |
|:--|:--|
| `REPORT.pdf`, `REPORT.typ` | the results report |
| `data/*.json` | one `ic prime` report per run (255) |
| `data/COMMANDS.txt` | every command that produced them |
| `graphs/plot_secg.py`, `graphs/fig{1_secg,2_s3_cost,3_raw_base}.{svg,png}` | the figures and their source |
| `graphs/make_tables.py` | prints the report's table rows from `data/` |

## Limits

- Two-point sums only; no linear algebra over the relations and no descent.
- S₃ times: 20 targets per run, one thread, single runs on a shared machine; rho charged √n affine additions.
- The base is `[h]` of the small-x points; the raw variant is sampled (4000 pairs, B = 64, ten curves) and the `sect233k1` excess is unexplained.
- Toy curves validate the method, not any SECG curve.
- Outside this study: two `curve_traits::subfield` tests failed at the merge base `45efe3b33` too, from a `debug_assert` comparing elements of different degree; fixed (one line in `subfield.rs`).
