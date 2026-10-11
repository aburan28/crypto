# The scan's stages on v3: R04's probes on main's head

**What this is.** A stage diagnostic: where a scanned summand's time goes
on v3, main's head `995ea207`, measured with R04's instrument. It ran on
2026-10-06 from 05:52 to 05:57 UTC, after R07's last step and R02b's
go/no-go. It is a record, not a round's measurement: it prices stages,
not cold time, and it decides nothing by itself.

## The instrument

- **R04's probe commit** (`--features scan-probes`), cherry-picked onto
  `995ea207` as `65cc4f6e` without change: its diff adds and removes R04's
  commit `9a48b389`'s 175 lines, line for line.
- **The probes** read the time-stamp counter at the scan's stage
  boundaries: subtract, key, filter, admitted. `trial` is the rest of each
  collection trial, outside the scan. `ic price` adds the totals, with the
  counter's rate measured in-process, to its report as `scan_probes`.
- **The binary,** `ic-r04probes-on-995ea207-65cc4f6e`, is hashed in
  [`binaries.sha256`](binaries.sha256).

## What ran

- `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
  `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`,
  two rounds: 12 processes ([`probe.sh`](probe.sh)).
- Each was `ic price --single-target`, `RAYON_NUM_THREADS=1`, isolated on
  CPU 2.
- **Two were marked contended,** both at `2^44.3` in round 1, and are
  left out. That leaves two processes at `2^44.3` and four at each other
  size.
- **Every process recovered** its row's logarithm, the one v3's own
  runs recover.

## The stages

Nanoseconds a scanned summand: each figure is the mean over the clean
processes ([`stages.tsv`](stages.tsv), from [`stages.sh`](stages.sh) on
[`runs.tar.xz`](runs.tar.xz)). In brackets, the stage's share of the
scan.

| curve | `log₂ r` | subtract | key | filter | admitted | scan | admitted keys a summand |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 12.0 (34%) | 11.3 (32%) | 6.1 (17%) | 5.7 (16%) | 35.2 | 0.028 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 11.3 (33%) | 11.6 (33%) | 6.0 (17%) | 5.8 (17%) | 34.8 | 0.032 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 12.1 (33%) | 12.7 (35%) | 7.4 (20%) | 4.3 (12%) | 36.4 | 0.020 |

**Against R04 on v0′** ([results](../../rounds/R04-scan-probes/README.md)),
read with care: a different binary, kernel build and day, though the
same host class.

| curve | R04 on v0′: subtract, scan | v3: subtract, scan |
|:--|:--|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 13.1, 54.9 | 12.0, 35.2 |
| `icv1-f2m59-tm943548413-98844ecc` | 32.5, 78.6 | 11.3, 34.8 |
| `icv1-f2m61-t158598901-ab42b6c5` | 33.3, 82.8 | 12.1, 36.4 |

## What it shows

- **The wide-tail premium is gone.** At `2^44.5` and `2^47.2` the
  subtraction fell from 32.5 and 33.3 ns to 11.3 and 12.1 ns, the narrow
  kernel's level. That is main's #1242, and it is why R02b's kernel no
  longer gains ([its go/no-go](../R02b-go-no-go-20261006/README.md)).
- **The admitted keys fell** from 0.17–0.21 a summand to 0.020–0.032,
  which is R05's sharper filter.
- **On v3 the scan is about 35 ns a summand at all three sizes,** in
  three roughly equal parts:
  - the subtraction, 11–12 ns;
  - the key, 11–13 ns;
  - the filter and the admitted keys together, 10–12 ns.
- **The levers left in the scan follow.** The key is R06's target: its
  GFNI and funnel kernels price the key at about 3.2–3.7 ns alone. The
  subtraction is the folding kernel's.

## Its limits

- **The probes' overhead is not measured here.** No same-binary control
  ran, so the shares are reported, not trusted, as in R04.
- **Two processes a size, at most four,** give means, not intervals.
- **Wall time on this container** carries the ±5–10% residual noise
  AGENTS.md §10 describes.
