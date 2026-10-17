# The scan's stages on v8: R04's probes after R06, R08, R11 and R12

**What this is.** A stage diagnostic: where a scanned summand's time goes on
v8, R12's accepted baseline, measured with R04's instrument. It ran on
2026-10-17 from 03:37 to 03:43 UTC on the Emerald Rapids host R12 ran on.
It is a record, not a round's measurement: it prices stages, not cold time,
and decides nothing by itself.

## The instrument

**The probe commit:** R04's probes (`--features scan-probes`), ported onto
v8 (`23138d20`) as `786c06f7` and frozen as [`probes-v8.patch`](probes-v8.patch).
The port started from R08's exploration's
[`candidate-probes.patch`](../R08-fold-kernel-20261006/candidate-probes.patch),
which conflicted in two places:
- **`Cargo.toml`:** v8's newer features are kept, and `scan-probes` is added.
- **`koblitz_index_calculus.rs`:**
  - The probe module is kept, beside v8's `MAX_N = 127`.
  - The filter's lap is placed before R12's run-prefetch pass, so **v8's
    `admitted` stage includes that pass**. At v3 there was no such pass.

**Default builds are unchanged:** without the feature, every probe is a
`#[cfg]`-gated statement.

**The binary** is `ic-r04probes-on-v8-786c06f7`, hashed in
[`binaries.sha256`](binaries.sha256).

## What ran

- **The rows:** `M1`'s two rows at the three target sizes, four rounds: 24 processes
  ([`probe.sh`](probe.sh), [`probe.log`](probe.log)).
- **Each process:** `ic price --single-target` with the suite's rho seed and
  `RAYON_NUM_THREADS=1`, isolated on CPU 2.
- **Refusals:** the busy gate refused a run 13 times before it started. Each
  refused run waited and was asked for again.
- **Contended runs:** one process, at 2^47.2 in round 1, was marked contended
  and is left out.
- **Logarithms:** every process recovered its row's logarithm, the one v8's
  own runs recover in [A15](../A15-filter-slots-20261017/pairs.tsv).

## The stages

The figures are nanoseconds per scanned summand, each the mean over the clean processes
([`stages.tsv`](stages.tsv), from [`stages.sh`](stages.sh) on
[`runs.tar.xz`](runs.tar.xz)). Each v3 figure is in brackets, from
[the v3 diagnostic](../scan-stages-v3-20261006/README.md).

| curve | `log₂ r` | subtract | key | filter | admitted | scan | admitted keys a summand |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 4.35 (12.0) | 4.70 (11.3) | 4.95 (6.1) | 2.19 (5.7) | **16.2** (35.2) | 0.028 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 4.41 (11.3) | 4.77 (11.6) | 5.11 (6.0) | 2.50 (5.8) | **16.8** (34.8) | 0.032 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 4.40 (12.1) | 5.11 (12.7) | 6.56 (7.4) | 1.57 (4.3) | **17.6** (36.4) | 0.020 |

Where the time goes, as the mean over the clean processes
([`shares.sh`](shares.sh), [`shares.tsv`](shares.tsv)):

| curve | scan / collection | collection / cold | build / cold | selection / cold |
|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 0.974 | 0.560 | 0.303 | 0.058 |
| `icv1-f2m59-tm943548413-98844ecc` | 0.975 | 0.555 | 0.292 | 0.064 |
| `icv1-f2m61-t158598901-ab42b6c5` | 0.989 | 0.815 | 0.131 | 0.023 |

## What it shows

- **The scan costs 16–18 ns a scanned summand, against about 35 at v3.**
  - The subtraction fell from 11–12 to 4.4 ns (R08's folding kernel).
  - The key fell from 11–13 to 4.7–5.1 ns (R06).
  - The admitted keys fell from 4.3–5.8 to 1.6–2.5 ns (R12), even though
    that stage now includes R12's prefetch pass.
- **The four stages are now about equal, and the filter is the largest:** 5.0
  ns at 2^44.3 and 6.6 ns at 2^47.2. Its share grows with the filter's
  size, 4 MiB at 2^44.3 and 8 MiB at 2^47.2.
- **[A15](../A15-filter-slots-20261017/README.md) shows the filter is bound by memory, not by its hash.**
  Hashing each key once instead of twice made collection slower, not faster.
- **The scan is nearly all of collection,** and collection is 56% of cold
  time at 2^44.3–2^44.5 and 82% at 2^47.2. The build is 30% at the two
  smaller sizes and 13% at the largest.

## Its limits

- **The probes' overhead is not measured.** No same-binary control ran, so
  the shares are reported, not trusted, as in R04.
- **The v3 figures come from a different binary, day and stage
  boundary.** v8's `admitted` stage includes R12's pass.
- **Seven or eight processes a size** give means, not intervals.

## Graph and report

- **Graph:** [`report/stages.svg`](report/stages.svg), with its source
  [`report/plot.py`](report/plot.py).
- **Report:** [`report/report.pdf`](report/report.pdf), with its source
  [`report/report.typ`](report/report.typ).
