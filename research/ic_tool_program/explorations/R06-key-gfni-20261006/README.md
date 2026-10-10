# R06 exploration: the key's basis change by GFNI, with the funnel shifts

**What this is.** Stage timings and whole-process diagnostics for the
plan's A11, the scan's canonical key
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8). They ran on
2026-10-06, before R06 was declared again, to decide what R06 covers
and what to predict. It is a record, not a round's measurement. R06's
protocol discloses it, and none of its runs is pooled with R06's.

It follows the funnel exploration
([`../R06-funnel-key-20261005/README.md`](../R06-funnel-key-20261005/README.md)),
which priced the key's other half.

## The key's two halves

`FrobeniusCanon::canon_in_place` keys sixteen words at a time:
- **the basis change**, `coords`: eight byte-table lookups and XORs a key,
  one key at a time;
- **the least rotation** of the coordinates: `n − 1` chained
  rotate-and-min steps on eight keys a register, or, with the funnel
  candidate, one `vpshldq` and one minimum a rotation.

With the funnel kernel in place the basis change costs more than the
rotation, so the candidate takes both.

## The candidate's second commit: GFNI

`34ad98f9` (on the funnel commit `a809324f`, both on main's
`995ea207`).
- **The map is `F_2`-linear,** so byte `k` of a key's coordinates is the
  sum over the element's bytes `i` of an 8 × 8 block applied to byte
  `i`. `vgf2p8affineqb` applies a lane's 8 × 8 matrix to every byte of the
  lane.
- **Eight keys are transposed as bytes** (`vpermb`), so lane `i` holds
  byte `i` of all eight. Each lane is broadcast and multiplied by the
  column of blocks into every output byte; the products are XORed; the
  same transpose puts the keys back.
- **Sixteen keys** take two transposes in, sixteen broadcasts, sixteen
  affine products and XORs a byte, and two transposes out, where the
  tables take a lookup, a shift and an XOR per byte per key.
- **The blocks come from the byte tables,** which stay the basis's
  description: the GPU fold and the external tables read them.
- **Where it runs:** AVX-512BW, AVX-512 VBMI and GFNI, detected at run
  time; `KIC_CANON_COORDS=tables` keeps the tables as the same-binary
  control.

**The coordinates are identical.**
`the_gfni_basis_change_matches_the_tables` compares them with the
tables' at every degree from 2 to 63 on the field
`find_irreducible_sparse` gives: 2,000 random words, every single bit
(a column of the map each) and the edge words. It passes with the other
30 `koblitz_fast` tests, the funnel kernel's included
([`key-tests.log`](key-tests.log), 31 passed).

## 1. The basis change alone (`gfnibench/`)

`gfnibench` builds a random invertible map at each degree, applies it to
2^20 random words by the tables and by GFNI, compares every output, and
reports the best of seven passes in ns a key. One isolated run
(`isolated_bench run --cpus 2`, uncontended):

| `n` | tables | GFNI | identical |
|--:|--:|--:|:--|
| 41 | 2.88 | 0.82 | yes |
| 53 | 3.79 | 0.94 | yes |
| 59 | 4.32 | 0.93 | yes |
| 61 | 4.31 | 0.89 | yes |
| 63 | 4.31 | 0.89 | yes |

**GFNI saves 2.1–3.4 ns a key,** about as much as the funnel kernel
saves on the rotation (2.2–2.5 times faster, 3–4 ns a key).

## 2. Whole processes on main (`995ea207`)

**What ran.**
- `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
  `icv1-f2m59-tm943548413-98844ecc` and
  `icv1-f2m61-t158598901-ab42b6c5`, three rounds.
- Three arms, the order rotating by round:
  - `base`: main's head, `995ea207` (v3's binary);
  - `key`: the candidate, `34ad98f9`;
  - `funnel`: the candidate's binary with `KIC_CANON_COORDS=tables`,
    the funnel kernel alone.
- Each process isolated (`isolated_bench run --cpus 2`), between R07's
  timed processes under the benchmark lock.
- A pair with a contended process is left out.

**Every logarithm matched** across the three arms on every row and round.
Five pairs were left out as contended (`pairs3.sh` names them), all from
the first two rounds, while ad-hoc checks ran beside them.

**Base over candidate, geometric mean over the clean pairs:**

| curve | pairs | cold, key | cold, funnel alone | collection, key | build, key | per-pair spread (sd of log), key |
|:--|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 5 | 1.118 | 1.024 | 1.175 | 1.050 | 0.084 |
| `icv1-f2m59-tm943548413-98844ecc` | 5 | 1.099 | 1.040 | 1.170 | 1.055 | 0.088 |
| `icv1-f2m61-t158598901-ab42b6c5` | 6 | 1.188 | 1.083 | 1.247 | 1.032 | 0.052 |

- **GFNI carries most of it.** Over the funnel alone, the key reads
  1.09, 1.06 and 1.10 at the three sizes.
- **The build gains little** (1.03–1.06): it keys every stored pair, but
  most of its time goes to moving the pairs into place, not to keying
  them.
- **The base arm's shares** (v3, medians): collection is 60.0%, 59.9% and
  83.2% of cold time, and the build 31.7%, 31.1% and 13.7%. A cold run
  keys 17.4 M, 17.5 M and 80.6 M words: every scanned summand and every
  stored pair.
- **The funnel alone reads lower than in the funnel exploration** on main
  (1.015, 1.086 and 1.068 there, six pairs each). Both are small samples;
  neither is R06's measurement.

## What it sets in R06's protocol

- **The candidate is both commits,** the funnel and GFNI: the key is one
  function and the round measures it whole.
- **The target sizes are all three,** each predicted at about 1.10–1.19.
- **The prediction:** about 1.12×, 1.10× and 1.19× cold time at the three
  sizes; nothing predicted to regress elsewhere, where the key is a
  smaller share.
- **The power:** a spread of 0.05–0.09 a pair gives a 95% half-width of
  about 1.6–2.8% at forty pairs a set.

## Files

- `gfnibench/`: the bench, a crate of its own (`cargo build --release`),
  its output and its isolation record.
- `key-tests.log`: the candidate's `koblitz_fast` tests.
- `runs.tar.xz`: the whole processes, each with its isolation record,
  and the script that ran them.
- `SHA256SUMS`, `binaries.sha256`: the files and the binaries.

The host is the one in R07's manifest: `6.18.44-fc-v70`, the 4-vCPU
Xeon with AVX-512, VBMI2, VPOPCNTDQ and GFNI.
