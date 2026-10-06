# R06: the scan's canonical key, by GFNI and funnel shifts

**Declared 2026-10-06, before any R06 timed run.** This is Track A, the
collection scan's canonical key (A11) in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8. Two explorations
ran first and are disclosed below. Nothing below changes after the first
timed run, except by a dated amendment appended at the end.

## Hypothesis

**The key is two kernels, and each costs about as much as the other.**
Each scanned summand, and each stored pair, is keyed by the least
rotation of its normal-basis coordinates (`FrobeniusCanon`):
- **The basis change** is eight byte-table lookups and XORs a key: about
  32 instructions, 2.9–4.3 ns a key alone (the exploration's
  `gfnibench`).
- **The least rotation** is `n − 1` chained steps of two shifts, a mask
  and a minimum on eight keys: 5.4–7.1 ns a key alone (the funnel
  exploration's `kbench`).

**The candidate replaces both, with the same keys:**
- **GFNI for the basis change.** The map is `F_2`-linear, so byte `k` of
  the coordinates is a sum over the element's bytes of 8 × 8 blocks,
  which `vgf2p8affineqb` applies to every byte of a lane. Eight keys are
  transposed as bytes (`vpermb`), each byte lane is broadcast and
  multiplied by its column of blocks, and the same transpose puts the
  keys back: 0.82–0.94 ns a key alone.
- **Funnel shifts for the rotation.** Each rotation's window is one
  `vpshldq` (VBMI2) of the coordinates written twice, compared
  top-aligned: 2.4–2.8 ns a key alone.

**The keys are identical,** so every count and every answer is
unchanged:
- `the_gfni_basis_change_matches_the_tables` compares the GFNI
  coordinates with the tables' at every degree from 2 to 63, on random
  words, every single bit and the edge words;
- `the_funnel_kernel_matches_the_chained_one` compares the two rotation
  kernels and the scalar search at every `n` from 2 to 63;
- `canon_many_is_canon_elementwise` compares the bulk keys, which now
  take both kernels, with the scalar key.

## The phase and its share

- **Collection is 60.0%, 59.9% and 83.2% of cold time** at
  `icv1-f2m53-tm56619371-dac20a85`, `icv1-f2m59-tm943548413-98844ecc` and
  `icv1-f2m61-t158598901-ab42b6c5`, and the build 31.7%, 31.1% and
  13.7%: v3's binary in the key exploration's base arm (medians, `M1`'s
  rows). R07's own figures, from many more runs, are the record.
- **Every scanned summand is keyed once,** and every stored pair once,
  so a cold run keys 17.4 M, 17.5 M and 80.6 M words at the three sizes.

## Class

Engineering, if it gains: the keys, and with them every count and every
answer, are unchanged by construction.

## The candidate

**Frozen here as [`candidate.patch`](candidate.patch)** (SHA-256 in
`candidate.patch.sha256`): two commits against v3 `995ea207`.
- `a809324f`, the funnel kernel: the exploration's `funnel.patch`,
  unchanged.
- `34ad98f9`, the GFNI basis change.

**Where each runs, detected at run time:**
- the funnel kernel where the CPU has AVX-512 VBMI2;
- the GFNI basis change where it has AVX-512BW, AVX-512 VBMI and GFNI;
- everywhere else, the code that runs today.

**Same-binary controls:** `KIC_CANON_KERNEL=chain` keeps the chained
rotation, and `KIC_CANON_COORDS=tables` the byte tables.

**What does not change:**
- the byte tables stay the basis's description: the GPU fold and the
  external tables read them, and the GFNI blocks are made from them;
- rho's canonicalisation (`canon_with_shift`, one point at a time), so
  the reference arm runs exactly what it ran.

## Hardware class

The result holds for x86-64 with AVX-512 VBMI2 and GFNI (this host's
Xeon). On any other class the candidate runs the code it runs today, and
no claim is made there.

## Arms

- **The base:** v3, `995ea207`, if R07 accepts it.
- **The candidate:** the base plus `candidate.patch`, built with
  `IC_BUILD_COMMIT`.
- **The manifest** records both binaries' SHA-256 and the commits they
  were built from.

## Before any timed step

1. **The tests,** on the candidate's tree:
   `cargo test --release --lib -- koblitz_fast`, which includes the three
   tests above.
2. **The pin, untimed:** the candidate on all 90 suite rows. Every output
   must equal v0's.
3. **The A/A:** the base against a byte-identical copy on `M1`'s 22 rows,
   five rounds.

## Rows

- **The suite rows:** `M1`'s 22 rows, and the other six suite rows at
  each target size: 40 rows, five rounds, ABAB, isolated, 400 processes.
- **The fresh holdouts:** eight rows at each target size, by suite v1's
  own construction (`icprog holdouts`):
  - recipe seeds 218 to 221;
  - two targets a seed, `T127` to `T134`;
  - `public_hash_seed` 23127 to 23134, and rho seeds `0x230000 + 127` to
    `+ 134`.

  They are committed in [`holdouts/`](holdouts/) with their
  `SHA256SUMS`. No round and no exploration has run them. 24 rows, five
  rounds, ABAB, isolated: 240 processes.
- **The extension.** A set whose interval half-width exceeds 3% after
  five rounds gets rounds 6–10, pooled.

## Power

The key exploration's per-pair spread of the log ratio was 0.05–0.09.
At forty pairs a set that gives a 95% half-width of about 1.6–2.8%.

## Prediction

| curve | `log₂ r` | key exploration, base over key | base over funnel alone | predicted |
|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.118 (5 pairs) | 1.024 | **1.08–1.15×** |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.099 (5 pairs) | 1.040 | **1.07–1.13×** |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.188 (6 pairs) | 1.083 | **1.15–1.22×** |

Elsewhere the key is a smaller share of cold time; no size is predicted
to regress.

## Success and stop

**Accepted** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At each target size, the paired cold-time ratio's 95% interval lies
   above 1.03, on the suite rows and on the fresh holdouts separately.
3. No size regresses beyond its A/A band.

**Rejected** otherwise. **Stopped** if a test fails, an output differs,
or a verification fails.

**If accepted**, the candidate becomes v4, with its binary hash, host and
`S` from this round's runs.

## Inadmissible

- Reusing either exploration's runs, pooling them with R06's, or
  choosing between them.
- Changing either kernel, the suite, the recipes or the targets, or
  drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting the key, the scan or collection as the speedup: the speedup is
  cold time.
- Crediting the gain to one kernel: the round measures them together;
  the exploration's same-binary control splits them, as a diagnostic.

## Cost

- **The tests, the pin and the A/A:** about 90 untimed and 220 timed
  processes, about an hour and a half.
- **The suite rows and holdouts:** 640 timed processes, about three and
  a half hours, and up to 640 more if extended.

## Order

- **R06 runs on v3, after R07's decision** and the rule's comparison at
  v3. If R07 rejects v3, R06 waits for the baseline that replaces it.
- **The commands** are R07's, with `r06` for the round:
  `icprog run r06 <step>` for `manifest`, `pin`, `aa`, `compare`,
  `holdout` and `extend`, then `icprog analyse r06`.

## Explorations before this declaration (disclosed)

- [`../../explorations/R06-funnel-key-20261005/README.md`](../../explorations/R06-funnel-key-20261005/README.md):
  the funnel kernel alone, the invariant filter set aside, and whole
  processes on v2 and on main.
- [`../../explorations/R06-key-gfni-20261006/README.md`](../../explorations/R06-key-gfni-20261006/README.md):
  the GFNI basis change alone (`gfnibench`), and whole processes on main
  with three arms: v3, the candidate, and the candidate with
  `KIC_CANON_COORDS=tables` (the funnel alone).
