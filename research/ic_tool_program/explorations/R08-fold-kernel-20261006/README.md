# The folding kernel: the scan's subtraction by carry-less folds, on R06's candidate

**What this is.** Timed explorations of the folding kernel, the lever on
the scan's subtraction (A2 in `research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§8). They ran on 2026-10-06, before any round was declared on the kernel,
on R06's candidate (`34ad98f9`: v3 plus the key's GFNI basis change and
funnel shifts). They are a record, not a round's measurement. None of their
runs is pooled with any round's, and R08's declaration discloses them
([protocol](../../rounds/R08-scan-fold/PROTOCOL.md)).

**The outcome.** The candidate cuts cold time by 18–27% at the three
largest sizes, with every logarithm unchanged. Its subtraction costs
4.7–4.9 ns a scanned summand, against 11.3–12.1 ns on v3. **R08 is
declared on it,** on v4, R06's candidate, if R06 accepts it.

## The candidate

Two commits on `34ad98f9`, frozen as R08's
[`candidate.patch`](../../rounds/R08-scan-fold/candidate.patch):
- `21f8609a`: the folding kernel, its two chains and the coordinate
  slices, with the fused scan;
- `a4c19ed7`: the build's rows past their own orbit on the same kernel.

The branch was named `r09-sub` before the round had a number, so the
binaries' file names say `r09`. They are R08's candidate.

**What changes, on a CPU with AVX-512F and VPCLMULQDQ:**
- **The reduction.** Each product reduces by two carry-less folds, eight
  lanes at once: six `vpclmulqdq`, two `vpternlogq` and an unpack for
  eight products. That is the scalar `Gf2`'s own reduction, vectorised.
  The 128-bit products carry every bit, so the kernel serves any field
  that folds, the wide tails of `2^59` and `2^61` included.
- **Montgomery's trick in two chains.** The groups' first and second
  halves invert as two chains, so the forward pass does not wait a
  multiply's latency per group.
- **Coordinate slices.** The scan hands the kernel the negated base's `x`
  and `y` as two arrays. Each rest's abscissa goes straight into the key
  buffer, and only the rests the filter admits are completed.
- **The build's rows past their own orbit** take the same kernel.

**The arms,** all from two binaries:

| arm | binary | what runs |
|:--|:--|:--|
| `base` | `ic-r06key-on-995ea207-34ad98f9` | R06's candidate |
| `sub` | `ic-r09sub-on-34ad98f9-a4c19ed7` | the candidate |
| `fold` | the same, `KIC_SCAN_SOA=0` | the folding kernel on the point list everywhere: no slices, no fused keys |
| `shift` | the same, `KIC_SCAN_FOLD=shift` | the shift kernel, which takes no slices: the base's paths in the candidate's binary |

The binaries are hashed in [`binaries.sha256`](binaries.sha256). Every
run's report carries its binary's BLAKE3, and each arm's runs carry one
hash.

## What ran

Every process was `ic price --single-target` with the suite's rho seed,
`RAYON_NUM_THREADS=1`, isolated on CPU 2 by the isolation tool. The rows
were `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
`icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`.

1. **The smoke run** ([`smoke.sh`](smoke.sh)): one process an arm, on
   `M1-T01` at each size, 12 processes, 04:45 to 04:49 UTC.
   - R07's holdout step was running then. Every smoke process waited for
     the isolation tool's exclusive lock, so none overlapped a timed R07
     process, and none was marked contended.
   - It decided only that the four-arm run was worth its time.
2. **The four-arm run** ([`explore.sh`](explore.sh)): four rounds, the
   arms' order rotating as a Latin square, 96 processes.
   - It ran from 06:26 to 08:02 UTC.
   - The container restarted three times (07:05, 07:23 and 07:53 UTC). A
     process killed by a restart left no isolation record and no report.
     The script, resumed, ran it again from the start.
3. **The scan's stages on the candidate** ([`probe.sh`](probe.sh)): R04's
   probes on the candidate, two rounds, 12 processes, 08:02 to 08:04 UTC.
   - The probe commit `843d03ce` is R04's probe commit, cherry-picked onto
     `a4c19ed7`. Its diff, [`candidate-probes.patch`](candidate-probes.patch),
     adds and removes R04's lines, line for line.
   - The binary is `ic-r09probes-on-34ad98f9-843d03ce`.
4. **The portable path** ([`portable.sh`](portable.sh)): `base` and `sub`
   with `KIC_SCAN_SIMD=0`, three rounds, 36 processes, 08:21 to 08:29
   UTC. Its section is below.

**Contended runs:**
- In the four-arm run, 10 of the 96 processes were marked contended: two
  `base`, two `sub`, four `fold` and two `shift`. Pairs with a contended
  run are left out of the clean figures, since contended runs are not
  pooled with clean ones.
- In the probe run, two of the 12 were marked contended and are left out.

**Every pair recovered the same logarithm** in both arms, all 72 of them
([`pairs-sub.tsv`](pairs-sub.tsv), [`pairs-fold.tsv`](pairs-fold.tsv),
[`pairs-shift.tsv`](pairs-shift.tsv)).

## Results

**Cold time,** the median repetition's set-up plus online interval. Each
figure is `base` over the arm, so above 1 the arm is faster. The intervals
are 95% t intervals of the geometric mean, on the clean pairs, with the
number of clean pairs in brackets:

| curve | `log₂ r` | `sub` (the candidate) | `fold` | `shift` |
|:--|--:|:--|:--|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | **1.293** [1.259, 1.328] (5) | 1.057 [0.978, 1.142] (6) | 0.995 [0.927, 1.068] (6) |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | **1.223** [1.158, 1.292] (7) | 1.107 [1.012, 1.211] (4) | 1.021 [0.969, 1.075] (6) |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | **1.362** [1.268, 1.464] (8) | 1.137 [1.063, 1.215] (8) | 0.949 [0.908, 0.992] (8) |

On all eight pairs, contended included, `sub` reads 1.302, 1.246 and
1.362.

**The phases,** the median over each arm's clean processes, in ms
([`phases.tsv`](phases.tsv)):

| curve | arm | processes | cold | build | collect |
|:--|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | `base` | 7 | 716.4 | 227.4 | 407.5 |
| | `sub` | 6 | 543.7 | 200.3 | 268.3 |
| | `fold` | 7 | 665.0 | 244.6 | 346.0 |
| | `shift` | 7 | 727.7 | 220.4 | 434.7 |
| `icv1-f2m59-tm943548413-98844ecc` | `base` | 7 | 732.4 | 244.0 | 401.2 |
| | `sub` | 8 | 601.2 | 222.5 | 289.4 |
| | `fold` | 5 | 717.2 | 270.2 | 348.9 |
| | `shift` | 7 | 747.1 | 232.4 | 430.5 |
| `icv1-f2m61-t158598901-ab42b6c5` | `base` | 8 | 2,758.0 | 442.4 | 2,237.4 |
| | `sub` | 8 | 2,033.7 | 382.6 | 1,511.9 |
| | `fold` | 8 | 2,513.6 | 524.4 | 1,900.7 |
| | `shift` | 8 | 2,913.0 | 420.2 | 2,385.7 |

**The scan's stages on the candidate,** in ns a scanned summand: the mean
over the clean processes ([`stages.tsv`](stages.tsv), from
[`stages.sh`](stages.sh)). Beside them, v3's, measured by the same probes
on main's head
([record](../scan-stages-v3-20261006/README.md)).

| curve | processes | subtract | key | filter | admitted | scan | v3: subtract, key, scan |
|:--|--:|--:|--:|--:|--:|--:|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 3 | 4.89 | 4.81 | 6.24 | 5.95 | 21.89 | 12.03, 11.32, 35.18 |
| `icv1-f2m59-tm943548413-98844ecc` | 4 | 4.84 | 4.75 | 6.30 | 6.64 | 22.53 | 11.32, 11.63, 34.80 |
| `icv1-f2m61-t158598901-ab42b6c5` | 3 | 4.68 | 5.11 | 6.78 | 4.40 | 20.98 | 12.06, 12.68, 36.35 |

A cold run scans 13.7 M and 13.6 M summands at the two smaller sizes
and 74.4 M at `2^47.2`. The probes count over a process's three repetitions,
41.2 M, 40.7 M and 223.1 M.

## Reading

- **The candidate is the gain, at every size.** Collection falls by
  28–34%, from 401–408 ms to 268–289 ms at the two smaller sizes, and
  from 2,237 to 1,512 ms at `2^47.2`. The build falls by 9–14%.
- **The subtraction is now as cheap as the key.** Both cost about 5 ns a
  summand. The key's fall from v3 is R06's, which is in the base; the
  subtraction's is this candidate's. No probe ran on the base, so the
  base's own subtraction is v3's figure, by construction, not by
  measurement.
- **The slices carry more than half the gain.** The folding kernel on the
  point list (`fold`) gains 6–14%. The kernel, the slices and the fused
  keys together (`sub`) gain 22–36%. Collection says the same: `fold`
  collects in 346–1,901 ms, and `sub` in 268–1,512 ms.
- **The point list costs the folding kernel its build gain.** In `fold`
  the build is slower than the base's at every size (245 against 227 ms
  at `2^44.3`, 524 against 442 ms at `2^47.2`). The build's rows take
  their slices only in `sub`.
- **The filter and the admitted keys are now the larger half of the
  scan:** 11.2–12.9 ns of 21.0–22.5 ns. Each admitted key costs
  210–225 ns, the price of one trip to memory. That is A13's target.
- **The fallback costs a little at `2^47.2`.** `shift` runs the base's
  paths in the candidate's binary. At the two smaller sizes it matches
  the base: 0.995 and 1.021, both intervals across 1. At `2^47.2` it reads
  0.949 [0.908, 0.992]. Its collection is above the base's on every one
  of the eight pairs, by 0.1–16%.
  - There the shift kernel does not apply, so both binaries run the
    scalar `Gf2` path, whose code the candidate does not change.
  - Why the candidate's binary runs that path more slowly is not
    established here.
  - The portable-path check below measures what a CPU without the vector
    kernels would see.

## The portable path

**What a CPU without the vector kernels runs,** measured on this host by
turning them off: R06's candidate (`base`) against the candidate (`sub`),
both with `KIC_SCAN_SIMD=0`, which turns off the scan's vector subtraction
and the vector key alike ([`portable.sh`](portable.sh)).
- It ran from 08:21 to 08:29 UTC: `M1`'s two rows at the three sizes,
  three rounds, the order alternating, 36 processes, none contended.
- Every pair recovered the same logarithm
  ([`pairs-portable.tsv`](pairs-portable.tsv)).

Cold time, `base` over `sub`, 95% t intervals of the geometric mean, six
pairs a size, with the medians over each arm's six processes:

| curve | `log₂ r` | cold time | 95% interval | cold, ms (`base`, `sub`) | collect, ms (`base`, `sub`) |
|:--|--:|--:|:--|:--|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 0.926 | [0.827, 1.036] | 1,260, 1,373 | 879, 960 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.059 | [0.992, 1.131] | 1,441, 1,321 | 986, 926 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 0.984 | [0.948, 1.021] | 5,927, 6,007 | 5,229, 5,215 |

- **No size shows a portable regression.** Every interval contains 1,
  and the point estimates fall on both sides of it.
- **Nor does any rule out one of about 5%.** Six pairs a size cannot
  resolve that, and the `shift` arm's 7% at `2^47.2` above is not
  explained. A round that claims no cost on other CPUs would need its own
  portable rows.
- **The portable path is 2.2–2.9 times slower** than the vector path on
  this host: 1,260–1,441 ms against 540–600 ms at the two smaller sizes,
  and 5.9 s against 2.0 s at `2^47.2`.

## Its limits

- **Two rows a size,** in four rounds, with five to eight clean pairs a
  size: enough to declare, not to accept. R08's rows and holdouts are the
  measurement.
- **The probes' overhead is not measured.** As in R04, no same-binary
  control ran, so the stage shares are reported, not trusted.
- **One host class.** Everything here ran on this host's Xeon, which has
  AVX-512F and VPCLMULQDQ. No other class was measured, apart from the
  portable-path check on this host.
- **Wall time on this container** carries the ±5–10% residual noise that
  AGENTS.md §10 describes.

## Files

- [`runs.tar.xz`](runs.tar.xz) holds each process's report, isolation
  record and logs, under the suite's directory names:
  - `runs/` holds the four-arm run;
  - `smoke/` holds the smoke run;
  - `probes/runs/` holds the probes;
  - `portable/` holds the portable check;
  - `logs/` holds the scripts' logs.
- **The tables come from it:**
  - [`pairs-tsv.sh`](pairs-tsv.sh) gives the `pairs-*.tsv` files, as
    `pairs-tsv.sh runs base sub` and `pairs-tsv.sh portable base sub`;
  - [`pairs2.sh`](pairs2.sh) gives the intervals;
  - [`phase-medians.sh`](phase-medians.sh) gives [`phases.tsv`](phases.tsv),
    as `phase-medians.sh runs base sub fold shift`;
  - [`stages.sh`](stages.sh) gives [`stages.tsv`](stages.tsv), as
    `stages.sh probes/runs`.
- [`SHA256SUMS`](SHA256SUMS) checks every file here.
