# Small-x factor base on the SECG curves: exact two-point relation yield

Run `2026-10-14`, seed 7, release `ic`. Report: [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)).

## Answer

For a small-x base of `B` abscissae (`2B` points) on a SECG curve, the probability that a
random element of the prime-order subgroup is a signed sum of two base points is exactly
`2B²/(n − 1)` (`n` the subgroup order), because no pair sums coincide on 18 of the 21
curves measured (`sect163k1`, `sect233k1`, `sect239k1` have 2 to 20 collisions among 2²⁴
sums, a relative effect of order 10⁻⁶). Cofactor `h > 1` is handled by taking `[h]P` of each
small-x point, so the base and the targets lie in the subgroup.

| curve (B = 4096) | `log₂ n` | exact yield | trials for `B + 1` relations | rho (`log₂ √n`) |
|:--|--:|--:|--:|--:|
| `sect113r1`, `sect113r2` | 112.0 | 2⁻⁸⁷·⁰ | 2⁹⁹·⁰ | 2⁵⁶·⁰ |
| `sect163r1`, `sect163r2` | 162.0 | 2⁻¹³⁷·⁰ | 2¹⁴⁹·⁰ | 2⁸¹·⁰ |
| `secp160k1`, `secp160r1`, `secp160r2` | 160.0 | 2⁻¹³⁵·⁰ | 2¹⁴⁷·⁰ | 2⁸⁰·⁰ |

All 21 curves (`secp112r1/r2`, `secp128r1/r2`, `sect131r1/r2`, `sect163k1`, `secp192k1/r1`,
`sect193r1/r2`, `sect233k1/r1`, `sect239k1` besides the above) are in `REPORT.pdf`. Two-point
small-x sums lose to rho on every curve: the trials for a full system are about `n/(2B)`,
which reaches `√n` only for a base of about `√n/2` abscissae. Monte Carlo finds 0 hits in
66 million trials, as the yields require. On toy binary curves the exact count is inside the
95% interval of the Monte Carlo rate in 14 of 15 runs (the miss is 2.1 standard errors), and
all planted decompositions were detected and recovered in all 81 runs.

"61 or 63, 160 bit" was read as `sect163r1/r2` (163-bit) and, as the alternative reading,
`secp160k1/r1/r2`; both are measured.

## Use

```sh
ic list                                       # every SECG curve, plus ecc2k-95/130
ic inspect --curve sect113r1                  # (also: ic sect113r1)
ic prime --curve sect113r1 --size 4096        # binary, cofactor 2
ic prime --curve sect163r2 --size 4096
ic prime --curve secp160r1 --size 4096        # prime, 160 bits
ic prime --curve toy-binary:17 --size 128     # toy binary curve, hits are common
```

Code: `src/cryptanalysis/secg_standards.rs` (typed loader over
`docs/curves/standards/parameters.json`), `src/cryptanalysis/binary_small_x.rs` (four-limb
`GF(2^m)`, binary curve, small-x base, exact pair-sum count, toy curves),
`src/cryptanalysis/prime_small_x.rs` (generic ≤192-bit Montgomery prime field, cofactor),
`src/bin/ic/prime.rs`, `src/bin/ic/params.rs`.

## Files

| path | what it is |
|:--|:--|
| `REPORT.pdf`, `REPORT.typ` | the results report |
| `data/*.json` | one `ic prime` report per run (81) |
| `data/COMMANDS.txt` | every command that produced them |
| `graphs/plot_secg.py`, `graphs/fig1_secg.{svg,png}` | the figure and its source |
| `graphs/make_tables.py` | prints the report's table rows from `data/` |

## Limits

- Two-point sums only; the cost of deciding a target (a Semaev S₃ solve) is not included.
- The base is `[h]` of the small-x points; the variant that keeps raw small-x points is not measured.
- The cofactor-4 path has structural checks only, no toy statistical check.
- `ic prime` handles primes up to 192 bits and binary fields up to degree 256; larger SECG curves are named and inspectable only.
- At most 8192 abscissae (pair-sum memory); 8192 was run on three curves only.
