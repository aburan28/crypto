# The index calculus on two-word binary fields (B3b's design)

**Written 2026-10-01, before any B3b code.** The plan is
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9). B3b's protocol,
[`../rounds/B3b-two-word-kic/PROTOCOL.md`](../rounds/B3b-two-word-kic/PROTOCOL.md),
declares its measurement.

## 1. Where the tool stands

B3 lifted the rho reference to two-word fields:
- `GfWide`: binary fields up to degree 127 in a `u128`;
- `WideCurve` and `WideCanon`: the Koblitz curve and its signed-Frobenius
  classes on that field;
- `WideRho`: the matched walk, step for step `ParallelRho`'s.

The index calculus (`kic`) still stops at `n = 62`. Its fast paths are
one word throughout:
- **the curve**: `FastCurve` and `FastPoint`;
- **the table**: the pair table's `u64` keys;
- **collection**: the scan and the collector's `u64` scalars;
- **the descent**: its walk;
- **the linear algebra**: the sparse solver's `u64` residues, which need
  `r < 2^63`.

Its slow paths are already general: the curve's data, the factor base's
points, the orbit walk's BigUint version, the relations' rows and the
logarithm's recovery are BigUint.

## 2. What B3b adds

Two-word kernels for every phase `kic` runs:

| phase | two-word kernel |
|:--|:--|
| select | the subgroup-orbit sampler: its draw, the half-trace lift, the cofactor multiplication, and a `u128`-packed orbit walk |
| build | the folded pair table, its 128-bit canonical keys hashed into the one-word table's layout |
| collect | the `m = 3` walked scan, with a lazy batched subtraction, and the collector with `u128` scalars |
| logs | the log system for any `r < 2^127`; the sparse solver's residues in two limbs from `2^63` up |
| descent | the `m = 2` walk and the `m = 3` scan |

With them, `kic` runs **end to end at F0** on two-word fields wherever its
cost is affordable, which is where the subgroup is small. It does not
make the gate's subgroup (`E_0`, `r ≈ 2^81`) affordable at F0: that needs
F1, in B7b.

## 3. Two pipelines, one algorithm

**The one-word pipeline is left byte for byte as it is.**
- Its kernels are tuned to one word throughout.
- Its outputs are pinned on the suite.
- Every round's F0 timing depends on it.

Making it generic would touch every loop the suite times.

**B3b adds a second pipeline over the same algorithm**, as `WideRho` did for
rho. It makes the same choices everywhere:
- **the draw**: the one-word `u64` draw below degree 64; a `u128` draw
  only above it;
- **the lift**: the half-trace root `(x, x·H(c))` first, then its
  negation. `GfWide::solve_quadratic` may return the other root, so it is
  not used;
- **the order of things**: the same orbit order, row order and fold;
- **the canonical key**: `WideCanon::least_rotation`, whose normal
  element and tie rule are `FrobeniusCanon`'s up to 64 bits;
- **the table**: the same placement. The one-word table reads a key only
  through `pair_filter_hash`. The two-word hash equals it on every key
  below `2^64`, so the two tables are identical where both exist: the
  same bucket offsets, words, tags and filter;
- **the walks**: the same steps.

**The test of sameness is counter for counter.** At one-word degrees, where
both pipelines run, they must give:
- the same base and the same table;
- the same relations, in the same order;
- the same logarithms and the same counters.

The check runs on every declared instance, as `WideRho`'s test does
against `ParallelRho`.

## 4. Instances

The orders are `#E(GF(2^n)) = 2^n + 1 − t_n` by the trace recurrence;
`r` is the largest prime factor and `h = #E/r`.

| curve | `n` | `log₂ r` | `r` | `h` | role |
|:--|--:|--:|--:|--:|:--|
| `icv1-f2m67-tm19346764963-82c84cca` | 67 | 26.19 | 76589041 | 1926828573412 | F0 |
| `icv1-f2m67-t19346764963-e760df38` | 67 | 35.99 | 68352708293 | 2159006662 | F0 |
| `icv1-f2m79-tm420247971347-2a24b892` | 79 | 33.74 | 14377452373 | 42042421294532 | F0 |
| `icv1-f2m71-tm48653080717-f25c4638` | 71 | 49.06 | 588353361747061 | 4013206 | F0 |
| `icv1-f2m83-t6151469093347-cdcc5432` | 83 | 52.93 | 8569786107849059 | 1128547018 | F0, on the gate's field |
| `icv1-f2m83-tm6151469093347-debefd74` | 83 | 81.00 | 2417851639230796216685689 | 4 | the gate (`E_0`): F1 only, B7b |

**On the gate's field, not the gate's curve.** The fifth row is `K_1`
under the gate's modulus `z^83 + z^45 + z^2 + z + 1`: another Koblitz
model on the gate's field (AGENTS.md §8b). It shares the field and its
Frobenius module (`ord_83(2) = 82`), and neither the curve nor the
subgroup. Every report on it says so. Nothing measured on it discharges
the m = 83 gate (§8a).

All five F0 instances have prime degree, so no proper intermediate
subfield is involved. Their `r` are below `2^63`, so the existing `u64`
linear algebra serves the F0 runs. The two-limb residues are tested on
their own (§5), for F1 at `r = 2^69.7` and `2^81`.

## 5. Tests the code must carry

- The two-word kernels against the general arithmetic (`binary_ecc`):
  every field degree from 63 to 127, and the gate's modulus.
- **Sameness at one word**, as §3 defines it, on the declared instances.
- The two-limb residues against BigUint. A synthetic sparse system is
  solved modulo a 70-bit and an 81-bit prime, and the solution checked.
- Each two-word F0 instance at a reduced size, within the test budget.

## 6. What it does not do

- **No speed claim.** The two-word pipeline is a new path, measured as
  new rows. Its premium over one word is a stage diagnostic for B7b.
- **Not three words.** `n ≥ 128`, with ECC2K-130 at `n = 131`, waits for
  B4.
- **Not F1 at two words.** Wiring F1's samples to these kernels, and
  running them at `n = 83`, is B7b.
