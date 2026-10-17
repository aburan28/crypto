# Small-x factor base on NIST P-192: exact two-point relation yield

Run `2026-10-14`, seed 7, release `ic`. Report: [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)).

## Answer

For a small-x base of `B` abscissae (`2B` points) on NIST P-192, the probability that
a random point is a signed sum of two base points is exactly `2B²/(N − 1)`, because
none of the up to 67 million pair sums coincide. It is a count, not an estimate.

| `B` | pair sums | collisions | exact yield per trial | trials for `B + 1` relations |
|--:|--:|--:|--:|--:|
| 256 | 2¹⁶ | 0 | 2⁻¹⁷⁵·⁰ | 2¹⁸³·⁰ |
| 1024 | 2²⁰ | 0 | 2⁻¹⁷¹·⁰ | 2¹⁸¹·⁰ |
| 2048 | 2²² | 0 | 2⁻¹⁶⁹·⁰ | 2¹⁸⁰·⁰ |
| 4096 | 2²⁴ | 0 | 2⁻¹⁶⁷·⁰ | 2¹⁷⁹·⁰ |
| 8192 | 2²⁶ | 0 | 2⁻¹⁶⁵·⁰ | 2¹⁷⁸·⁰ |

Pollard rho is about 2⁹⁶ group operations. Monte Carlo on P-192 finds nothing (0 hits
in 2 000 000 trials, rate below 1.84 × 10⁻⁶), as it must. The method is validated on
toy primes, where the exact count equals the Monte Carlo rate to within 0.5% (17 of 17
inside the 95% interval) and every planted decomposition is detected and recovered.
The earlier estimate for a 2¹³-point base, 2⁻¹⁶⁶, was a factor of two high: it counted
a sum and its negative as separate abscissae.

## Use

```sh
ic secp192r1                       # inspect the named curve (also: ic inspect --curve secp192r1)
ic prime --size 4096               # exact yield, Monte Carlo and planted controls on P-192
ic prime --curve toy:20 --size 512 # a toy prime, where hits are common
ic --json prime ... --out report.json
```

Code: `src/cryptanalysis/prime_small_x.rs` (field trait, small-x base, exact pair-sum
count, measurement, toy curves, 6 unit tests), `src/bin/ic/prime.rs`, and the
`secp192r1` entry in `src/bin/ic/params.rs`.

## Files

| path | what it is |
|:--|:--|
| `REPORT.pdf`, `REPORT.typ` | the results report |
| `data/*.json` | one `ic prime` report per run (22) |
| `data/COMMANDS.txt` | every command that produced them |
| `graphs/plot_prime.py`, `graphs/fig1_yield.{svg,png}` | the figure and its source |

## Limits

- Two-point sums only; three-point sums are not measured.
- Yield, not solving: the cost of deciding a target (a Semaev S₃ solve) and any linear algebra are not included.
- The small-x base is the only family; at most 8192 abscissae (memory of the pair-sum set).
- Toy curves have odd, not necessarily prime, order. They validate the counting, not P-192.
