# R08: the scan's subtraction by carry-less folds

**Declared 2026-10-07, before any R08 timed run.** This is Track A, the
collection scan (A2) in `research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§8. Explorations ran first and are disclosed below. Nothing below changes
after the first timed run, except by a dated amendment appended at the
end.

## Hypothesis

**The scan's subtraction is a third of the scan on v3, and most of it is
the reduction.** R04's probes on main's head put the subtraction at
11.3–12.1 ns of a scanned summand's 35–36 ns
([record](../../explorations/scan-stages-v3-20261006/README.md)). The
8-lane kernel reduces each product by shifts and XORs over the tail's
terms. A shift whose count comes from a register is two micro-ops on this
core, and a pentanomial's tail spends sixteen of them a product.

**The candidate reduces by Gf2's own trick, eight lanes wide,** with the
same field elements lane for lane:
- **Two folds a product.** With one operand shifted up by `64 − n`, each
  fold is one `vpclmulqdq`, so eight products take six carry-less
  multiplies, two three-way XORs (`vpternlogq`) and an unpack.
- **Every tail, every degree.** The 128-bit products carry every bit, so
  the kernel serves any field that folds (`2·deg t ≤ n + 1`). That takes
  in `GF(2^59)` and `GF(2^61)`, whose tails reach `z^66` and which ran the
  scalar path before #1242.
- **Montgomery's trick in two chains,** the groups' first and second
  halves, so the forward pass is not one multiply's latency per group.
- **Coordinate slices of the negated base:** the scan hands the kernel
  `x` and `y` as two arrays, writes each rest's abscissa straight into the
  key buffer, and completes only the rests the filter admits.
- **The build's rows too.** A row's suffix past its own orbit, where no
  abscissa is the representative's, takes the same kernel. A row's own
  orbit keeps the point list, since its sums are degenerate.

**The outputs are identical,** so every count and every answer is
unchanged:
- `the_folding_kernel_matches_the_scalar_addition_on_slices` compares the
  kernel with the scalar addition at degrees on both sides of the old
  `n + deg t ≤ 65` bound;
- `the_folding_kernel_s_scans_hand_the_sink_what_the_point_list_does`
  compares the scans' witnesses, in order, with the point list's;
- the folded build's two-pass test now covers `2^59` and `2^61`, and the
  table is word for word the same.

## The phase and its share

- **Collection is 56.9%, 54.8% and 81.1% of cold time** on the base at
  `icv1-f2m53-tm56619371-dac20a85`, `icv1-f2m59-tm943548413-98844ecc` and
  `icv1-f2m61-t158598901-ab42b6c5`, and the build 31.7%, 33.3% and 16.0%.
  These are the medians of the exploration's `base` arm, `M1`'s rows; the
  round's own figures are the record.
- **A cold run scans 13.7 M, 13.6 M and 74.4 M summands** at the three
  sizes, and its build adds 3.7 M, 3.9 M and 6.3 M pairs, each of which
  takes the subtraction once.

## Class

Engineering, if it gains: the field elements, and with them every count and
every answer, are unchanged by construction.

## The candidate

**Frozen here as [`candidate.patch`](candidate.patch)** (SHA-256 in
`candidate.patch.sha256`), a plain diff from the base's tree, R06's
candidate `34ad98f9`, to `a4c19ed7`, two commits:
- `21f8609a`, the folding kernel, its two chains and the coordinate
  slices, with the fused scan;
- `a4c19ed7`, the build's rows past their own orbit on the same kernel.

Applying the patch to `34ad98f9`'s tree gives `a4c19ed7`'s tree exactly.

**Where it runs, detected at run time:** AVX-512F and VPCLMULQDQ. Every
other CPU, and any field that does not fold, runs the code that runs today.

**Same-binary controls:**
- `KIC_SCAN_FOLD=shift` keeps the shift kernel, which takes no slices, so
  the base's paths;
- `KIC_SCAN_SOA=0` runs the folding kernel on the point list everywhere.

## Hardware class

**R08's class is R06's** (R06's amendment 1): the programme's reference
class, x86-64 with AVX-512F, PCLMULQDQ, VPCLMULQDQ and GFNI, with AVX-512
BW, VBMI and VBMI2 for the base's key.
- **The result holds for that class.** The folding kernel needs
  AVX-512F and VPCLMULQDQ.
- **The runner enforces the class.** `icprog run r08` refuses every
  timed step on a host outside it, and `icprog host-class r08` reports
  the class against the host.
- **Every other CPU runs the code it runs today.** On this host, with the
  vector kernels turned off, the exploration's portable check found no
  size regressing, and could not rule out a cost of about 5%. No claim is
  made for any other class.

## Arms

- **The base:** v4, R06's candidate `34ad98f9`, if R06 accepts it, the
  binary R06 measured. If R06 rejects it, R08 does not run on this
  declaration: the candidate is ported to the baseline that replaces v4,
  and an amendment discloses the port before any run.
- **The candidate:** the base plus `candidate.patch`, built with
  `IC_BUILD_COMMIT`.
- **The manifest** records both binaries' SHA-256 and the commits they
  were built from.

## Before any timed step

1. **The tests,** on the candidate's tree:
   `cargo test --release --lib -- koblitz_fast koblitz_index_calculus`,
   which includes the three tests above.
2. **The pin, untimed:** the candidate on all 90 suite rows. Every output
   must equal v0's.
3. **The A/A:** the base against a byte-identical copy on `M1`'s 22 rows,
   five rounds.

## Rows

- **The suite rows:** `M1`'s 22 rows, and the other six suite rows at
  each target size: 40 rows, five rounds, ABAB, isolated, 400 processes.
- **The fresh holdouts:** eight rows at each target size, by suite v1's
  own construction (`icprog holdouts`):
  - recipe seeds 222 to 225;
  - two targets a seed, `T135` to `T142`;
  - `public_hash_seed` 23135 to 23142, and rho seeds `0x230000 + 135` to
    `+ 142`.

  They are committed in [`holdouts/`](holdouts/) with their
  `SHA256SUMS`. No round and no exploration has run them. 24 rows, five
  rounds, ABAB, isolated: 240 processes.
- **The extension.** A set whose interval half-width exceeds 3% after five
  rounds gets rounds 6–10, pooled.

## Power

The exploration's per-pair spread of the log ratio was 0.022, 0.059 and
0.086 at the three sizes. At forty pairs a set, that gives a 95%
half-width of about 0.7%, 1.9% and 2.8%. The interval clears 1.05 with
80% power when the true ratio is at least about 1.06, 1.08 and 1.09.

## Prediction

| curve | `log₂ r` | exploration, base over candidate | 95% interval | predicted |
|:--|--:|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.293 (5 pairs) | [1.259, 1.328] | **1.25–1.33×** |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.223 (7 pairs) | [1.158, 1.292] | **1.16–1.29×** |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.362 (8 pairs) | [1.268, 1.464] | **1.27–1.46×** |

Elsewhere the subtraction is a smaller share of cold time, and no size is
predicted to regress.

## Success and stop

**Accepted** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At each target size, the paired cold-time ratio's 95% interval lies
   above 1.05, on the suite rows and on the fresh holdouts separately.
3. No size regresses beyond its A/A band.

**Rejected** otherwise. **Stopped** if a test fails, an output differs,
or a verification fails.

**If accepted**, the candidate becomes v5, with its binary hash, host and
`S` from this round's runs.

## Inadmissible

- Reusing an exploration's runs, pooling them with R08's, or choosing
  between them.
- Changing the kernel, the suite, the recipes or the targets, or drawing
  new holdout seeds after a run.
- Pooling contended runs.
- Quoting the subtraction, the scan or collection as the speedup: the
  speedup is cold time.
- Crediting the gain to one part of the candidate: the round measures the
  kernel, the slices and the build's rows together; the exploration's
  same-binary controls split them, as a diagnostic.

## Cost

- **The tests, the pin and the A/A:** about 90 untimed and 220 timed
  processes, about an hour and a half.
- **The suite rows and holdouts:** 640 timed processes, about three
  hours, and up to 640 more if extended.

## Order

- **R08 runs on v4, after R06's decision,** on a host of its class. If
  R06 rejects, see Arms.
- **The commands** are R06's, with `r08` for the round:
  `icprog run r08 <step>` for `manifest`, `pin`, `aa`, `compare`,
  `holdout` and `extend`, then `icprog analyse r08`.

## Explorations before this declaration (disclosed)

All are in
[`../../explorations/R08-fold-kernel-20261006/README.md`](../../explorations/R08-fold-kernel-20261006/README.md),
on R06's candidate, on 2026-10-06:
- **The smoke run:** one process an arm, `M1-T01` at each size.
- **The four-arm run:** `base`, the candidate, and the candidate's two
  same-binary controls (`fold`: the kernel on the point list; `shift`:
  the base's paths). It gave the prediction above.
- **The scan's stages** on the candidate, by R04's probes: the
  subtraction at 4.7–4.9 ns a scanned summand, against v3's 11.3–12.1.
- **The portable path:** both arms with `KIC_SCAN_SIMD=0`.

Two levers explored on the candidate are not part of it, and are declared
separately if at all:
- the build's partition streams and per-run filters
  ([record](../../explorations/A1-build-on-fold-20261006/README.md));
- the admitted keys' run prefetch
  ([record](../../explorations/A13-admitted-prefetch-20261006/README.md)).
