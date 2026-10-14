# Scope: making index calculus faster end to end, and where it could beat rho

2026-10-13.  Planning note: nothing here is measured by this note.  Every
number is quoted from a committed record (cited) or is arithmetic on them.
It follows the ecbench ladder of
[`../pdp_root_finder_ecbench_20261013/`](../pdp_root_finder_ecbench_20261013/README.md)
and the extension-field ledger of
[`../notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md`](../notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md)
§11.

## 1. Where the end-to-end cost of index calculus sits today

From the ladder's sessions (`E_0` over `F_{2^n}`, standard subspace base of
dimension `d`, `rootfind-mitm:m=3`, verified runs only, means):

| `n` | `d` | base points | columns | trials | relations | hit rate | total GAE | cofactor | relations phase, share of `S` | oracle set-up share | `S` ÷ strong rho |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 13 | 5 | 33 | 17 | 30 | 15 | 0.51 | `1.85·10³` | 4 | 82 % | 16 % | 2.9 |
| 19 | 7 | 139 | 70 | 112 | 65 | 0.58 | `1.71·10⁴` | 4 | 71 % | 29 % | 15 |
| 23 | 8 | 275 | 138 | 373 | 129 | 0.35 | `1.01·10⁵` | 4 | 81 % | 19 % | 61 |
| 31 | 10 | 1,053 | 527 | 4,064 | 345 | 0.085 | `4.33·10⁶` | 1,492 | 94 % | 6 % | about 3,000 |

Three facts decide what can be sped up:

- **The relation phase is the cost.**  Factor base and linear algebra are
  under 2 % of `S` at every size.  Each trial costs `2A` point additions
  (`A` = abscissae), and the number of trials is `relations / hit rate`.
- **The hit rate is set by the whole group, not by `r`.**  A trial point
  decomposes with probability about `(base points)³ / (6·#E)`; that is 0.085
  at `n = 31`, as measured.  Rho costs `√r`.  On a curve with a large
  cofactor (`n = 31`: `#E = 1,492·r`) index calculus pays for the whole
  group and rho does not, which is most of the jump from 61× to about
  3,000× between `n = 23` and `n = 31`.  Read the ladder per size, never
  pooled.
- **Index calculus with any linear-scan oracle has total cost about
  `2^{2n/3}` at the best `d ≈ n/3`, against rho's `2^{n/2}`.**  No change
  inside the oracle alters that; it changes the constant.

## 2. Faster end to end: the levers, ranked

Each lever names what it changes (`S`, the counted unit, or wall time only),
its expected size, and its first falsifier.  None changes an exponent.

| # | lever | changes | expected | how to test | status of the code |
|--:|:--|:--|:--|:--|:--|
| L1 | **Choose `d` per size by the cost model** `T ≈ (cols/hit)·2A + A²/2`, not by the F6 ladder's convention | `S` | at `n = 31`, `d = 11` doubles the base and raises the hit rate about 8× (to about 0.68): about 1,000 trials of 2,106 additions plus a table of about `1.1·10⁶`, against 4,064 trials of 1,054 today, so about 1.3× lower `S`; small, but free | ecbench grid `d ∈ {n/3 − 1, n/3, n/3 + 1}` at `n = 23, 31`, same seeds | ready: a spec change |
| L2 | **Prefer curves with a small cofactor** | `S` | at `n = 31` the cofactor costs about `2^{10.5}` in hit rate | compare `E_0`/`E_1` at sizes with cofactor 2 or 4 | ready: curve choice; register the curves (AGENTS.md §11) |
| L3 | **Lazy ordinates in the scan**: compute `x(R ∓ P_k)` only, finish `y` on a key hit (`FastCurve::add_many_lazy`, which has an AVX-512 path) | wall | one multiplication of about seven per point saved, plus SIMD lanes: 1.3–2× on the relation phase | new oracle id `rootfind-mitm-lazy` (a changed algorithm is a new id), ladder rerun, wall admitted at L2 | primitive exists |
| L4 | **Frobenius-folded pair table** (`FrobeniusPairTable`, `mitm-frobenius`) under the batched scan | `S` (set-up), wall | table `2n×` smaller: set-up share 6–29 % falls to under 1 %; probes stay in cache | `rootfind-mitm-frobenius` id, same ladder | table exists; scan to port |
| L5 | **Exact block size against hit rate**: the block of 32 overshoots on hits (1.05–1.26× `mitm`'s charge); a block of `⌈1/hit⌉·k` or an adaptive one removes most of it | `S` | brings `rootfind-mitm` to within 1–2 % of `mitm` in `S` while keeping the wall gain | unit test of charges, ladder rerun | trivial |
| L6 | **Interleave the trials of the one target**: the trial points `R_t = [a_t]G + [b_t]Q` of a single solve are independent, so their scans can share SIMD lanes and memory traffic | wall | 1.5–3× on the relation phase | same as L3 | to build; strictly one target per run (the records already verify `[d]G = Q` for one `Q`) |
| L7 | **Threads** for the relation phase | wall only | near-linear in cores | ecbench runs one thread by contract; report separately | to build |

L1, L2 and L5 change the counted cost and are cheap to try; L3, L4 and L6
are wall-time engineering on the measured hot loop; L7 is a deployment
choice the harness deliberately excludes.  A realistic combined effect is
1.5–3× in `S` and 5–20× in wall time on this ladder.  The ratio to rho still
grows with `n`.

## 3. Where index calculus can beat rho: the extension-field programme

For `E(F_{p^k})` with the factor base on `F_p`-abscissae, Gaudry and Diem's
analysis gives `p^{2−2/k}` against rho's `p^{k/2}`: strictly better for
`k ≥ 3`.  The repository has measured this, not just cited it
(`RESEARCH_RESIDUAL_WALKS.md` §11):

| `k` | what is measured | where index calculus crosses rho | what blocks it |
|--:|:--|:--|:--|
| 3 | the whole pipeline, `C₃` (cost of one `S₄` solve) lowered in three rounds; double large prime | past `2^{230}` (§11.7) | `C₃` and the double-large-prime merge constant |
| 4 | `C₄` measured at 17× its floor, about `1.2·10¹²` `F_p` multiplications a solve (§11.17); linear algebra alone at 0.52× rho (§11.19) | about `2^{150}` (§11.17) | `C₄`, the cost of one decomposition |

At `k = 4` the linear algebra already beats rho (0.52×, eight curves, all
verified).  The relation phase is all that stands in the way, and its cost
is `trials × C₄`.  The crossover moves linearly in `C₄` in `p`, so as the
fourth power in group size: **every 10× off `C₄` moves the crossover down
by about `2^{13}`**.  That makes `C₄` the single number this programme
should attack, and it is exactly the requirements document's problem: an
exact solver for the zeros of one fixed-shape system, specialised enough
to be fast.

What the root-finding work brings to it, and what it does not:

- **Not the meet-in-the-middle oracle.**  Over an `F_p` base its per-trial
  cost is linear in `p` (the MITM relation phase measured 1,054–3,371× rho
  in §11.19).  Gaudry's per-trial cost is constant in `p`, which is why the
  exponent wins.  A fast oracle here must stay constant in `p`.
- **A compiled straight-line solver for the one shape (C1 + R8 of the
  requirements).**  Every trial's Weil-restricted `S₅` system has the same
  monomials and only new coefficients.  Record the elimination once (which
  rows, which pivots, which reductions to zero), then replay it per trial as
  a branch-free straight-line program over `F_p` (Joux–Vitse's F4Remake
  idea, taken to a compiled kernel with SIMD across trials of the one
  target).  Literature (Joux–Vitse) shows 2–3× over their own F4 and up to
  two orders over Magma on this shape.  Target: `C₄` 10–100× lower, that is
  the crossover from `2^{150}` to `2^{137}`–`2^{124}`.
- **Decomposition into `k − 1` points** (Joux–Vitse): a cheaper system per
  trial at a lower hit rate; worth one measurement against the compiled
  solver.
- **Symmetry**: the `D₄` invariants at `m = 4` (requirements C8; Bézout 4096
  → 512 on paper) cut the system before it is solved.

Honest scope: these curves are not deployed curves, and a crossover at
`2^{120}`–`2^{150}` group size is still a research result, not an attack
on anything in use.  The value is the first measured index-calculus win
over rho in this harness, on curves where the theory says it must exist.

## 4. Proposed programme

| step | what | output | size |
|--:|:--|:--|:--|
| 1 | L1 + L2 + L5 on the binary ladder (spec changes and a block-size fix) | ecbench sessions at `n = 23, 31, 37`, `S` per size against strong rho | 1 day |
| 2 | L3 + L4 as new oracle ids, ladder rerun | admitted wall comparisons against `ic-rootfind` | 2–3 days |
| 3 | `E(F_{p^4})` in ecbench: a curve construction for the §11.17 curves, the MITM and Gaudry oracles registered, a strong-rho reference | the §11.19 result reproduced inside ecbench with replays | 3–4 days |
| 4 | compiled straight-line `S₅` solver (record once, replay per trial; R8 + C1) | `C₄` measured against §11.17's `1.2·10¹²`; crossover re-derived | 1–2 weeks |
| 5 | if `C₄` falls 10× or more: the `k = 4` ladder in ecbench at the largest `p` that fits, against strong rho | the first index-calculus-below-rho ratio measured in the harness, or a measured miss | follows 4 |

Decision rule for step 4, fixed now: success is `C₄` at least 10× below
§11.17's measured value on the same `p` and seeds, exactness checked
against the MITM oracle on every solve; anything less is recorded as an
engineering constant and the programme stops at step 3.
