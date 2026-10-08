# Which isogeny classes over F_{p⁶} hold a Joux–Vitse weak curve, and where the weak curves sit

**Exhaustive census, toy sizes, 2026-10-07.**  Protocol frozen before the
run in [`PROTOCOL.md`](PROTOCOL.md); Addenda A and B registered on held-out
sizes before those sizes were read.  **Class: boundary for the registered
W3, advance for the walk stage only** (see §5).  `S`, end-to-end cost and
speedup are **unset**; nothing here concerns a prime-field or deployed curve.

## 1. Answer

1. **Depth-1 exclusion, exact at every size run.**  A class whose
   2-volcano has depth 1 (`v₂(f) = 1`, `f` the conductor of `Z[π]`) holds
   **no** weak curve.  That is half of all classes and 37% of all
   full-2-torsion curves, with zero exceptions over 2,065
   such classes and 1,929,837 curves (§3, A1).
2. **Depth ≥ 2 admission.**  Of the classes with `v₂(f) ≥ 2`, 86 to 95%
   hold a weak curve, rising with `p` (§3, A2).  The exceptions are not small classes: the
   random model gives them a chance of a few percent at most.
3. **Weak curves concentrate at height.**  Pooled over all classes, the
   weak fraction times `q` is about 2 at height 1, 4 at height 2 and 13.5
   at height 3, at every size (§3, A3).  For `p ≡ ±3 (mod 8)` nothing sits
   above height 3; for `p ≡ ±1 (mod 8)` (`p = 7, 17`) heights 4 to 6 are the
   densest of all, 15 to 27 per `q`, as registered in advance for `p = 17`
   (§3, A4, B1).  B2 (`p = 19`) is untested: that run hit the protocol's
   wall limit.
4. **The registered law W3 is falsified**: the triple
   `(v₂(f), D_K mod 8, |t| mod 16)` predicts but is not a function.  W4 as
   written fails at `p ≥ 11` because the deepest levels are empty for
   `p ≡ ±3 (mod 8)`; the non-uniformity it was reaching for is real and
   is item 3.

## 2. Method (`census.rs`, standalone Rust)

`F_{p⁶} = F_{p²}[s]/(s³ − c)` with `F_{p²} = F_p[i]/(i² − n)`.  Every
`λ ∈ F_{p⁶} ∖ {0, 1}` gives the Legendre curve `y² = x(x − 1)(x − λ)`;
one `λ` is kept per `j`, which is complete for full-2-torsion curves up to
`F_{p⁶}`-isomorphism and quadratic twist.  Per `j`: the weak flag
(`N_{F_{q³}/F_q}(c) = 1` for a cross-ratio `c ∈ {λ, 1 − λ, λ/(λ − 1)}`,
`q = p²`); the three rational 2-isogenies by Vélu on each 2-torsion point,
with the codomain's `j` and whether it keeps full 2-torsion; the trace by
baby-step giant-step on the order of random points, once per connected
component of the 2-isogeny graph and checked again on a second member of
every component (zero mismatches at every size); the height as the
shortest path to a floor node.  Classes are keyed by `|t|`, since a twist
has the same `j`, flag and graph.  `D = t² − 4q³ = f² D_K` is factored
exactly.  A 20-curve self-test agrees with brute-force point counting at
`p = 5` and `p = 7`.  Seed `20261007`.  Everything below is printed by
`./census tabulate results 5 7 11 13 17 19`.

## 3. Results against the registered predictions

| p | p mod 8 | full-2-torsion j | weak fraction × q (W1: 3) | in a weak class (W2) | classes | depth-1 classes, weak (A1) | depth ≥ 2 classes weak (A2) | weak × q at heights 1 / 2 / 3 (A3) | heights ≥ 4 (A4, B) | W3 conflicting triples | W5 |
|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|:--|
| 5 | 5 | 2,605 | 2.95 | 0.593 | 51 | 25 of 25 empty | 23 of 25 (0.920) | 2.12 / 4.39 / 10.00 | 0 of 20 | 1 of 17 | 0 mismatches, 50 of 50 |
| 7 | 7 | 19,609 | 3.00 | 0.562 | 148 | 74 of 74 empty | 63 of 73 (0.863) | 1.96 / 3.94 / 11.19 | 60/154, 12/23, 3/6 at heights 4, 5, 6 (19 to 26 × q) | 6 of 24 | 0, 147 of 147 |
| 11 | 3 | 295,261 | 2.99 | 0.588 | 606 | 303 of 303 empty | 271 of 302 (0.897) | 1.98 / 4.02 / 13.77 | 0 of 2,353 | 13 of 32 | 0, 605 of 605 |
| 13 | 5 | 804,469 | 3.00 | 0.591 | 1,015 | 507 of 507 empty | 464 of 507 (0.915) | 2.03 / 3.96 / 13.45 | 0 of 6,372 | 14 of 35 | 0, 1,014 of 1,014 |
| 17 | 1 | 4,022,929 | 3.00 | 0.608 | 2,313 | 1,156 of 1,156 empty | 1,099 of 1,156 (0.951) | 2.00 / 4.03 / 10.29 | 2,196/27,673, 321/3,442, 18/355 at heights 4, 5, 6 (23, 27, 15 × q); 0 of 24 at 7 | 8 of 39 | 0, 2,312 of 2,312 |
| 19 | 3 | 7,840,981 | 3.00 | partial | partial | partial | partial | partial | partial (B2 untested) | partial | partial |

Reading against the registration:

- **W1 passes** at every size: the weak fraction is `3/q` to two decimals.
- **W2 passes** at every size: `0.56`–`0.61`, inside `[0.45, 0.70]`, matching the sampled census.
- **W3 is falsified** at `p = 7, 11, 13, 17`: 6 to 14 of the 24 to 39 triples contain both weak and non-weak classes.  The conflicts are lopsided (for example `(2, 4, 6)` at `p = 11`: 36 weak, 1 not) but a function they are not.
- **W4 fails as written** at `p = 11, 13`: the height-1 fraction is `0.8×` the pooled fraction, not below `0.5×`, and the deepest level is empty.  The non-uniformity it reached for is A3.
- **W5 passes**: zero second-BSGS mismatches at every size, and the maximum height equals `v₂(f)` in every class with a floor.
- **A1 passes** on the held-out sizes `p = 13, 17` as registered, and on every size: 2,065 depth-1 classes and 1,929,837 curves with no weak curve.
- **A2 passes**: `0.915` and `0.951` on the held-out sizes.
- **A3 passes** on both held-out sizes, with `p = 17`'s height-3 value `10.29` at the bottom of its interval.
- **B1 passes**: `p = 17` has weak curves at heights 4, 5, 6, the densest levels of the census.
- **B2 is untested**: the first `p = 19` run was stopped at the protocol's 30-minute wall limit on a host at load average 411 from other sessions' jobs, after it had completed the enumeration (7,840,981 `j`, weak fraction `0.00831 = 3/q`, 1,963,180 components) but not the traces.  A second attempt, labelled `results_attempt2/`, was launched and is reported there if it finishes; it is not cited here.

Per-size tables (`v₂(f)`, `D_K mod 8`, classes, weak classes, nodes,
weak nodes, weak fraction × `q`; then height, nodes, weak nodes, weak
fraction × `q`) are in [`results/TABLES.md`](results/TABLES.md); the
per-class records are `results/census_p{p}.jsonl`.

Two further regularities, not registered, reported as observations:

- **2 unramified in `K` makes a class four times denser.**  At
  `v₂(f) = 3`, classes with `D_K ≡ 1, 5 (mod 8)` have weak fraction
  `7`–`10/q`; classes with `D_K ≡ 0, 4 (mod 8)` have `2`–`2.5/q`
  (`p = 11, 13`).
- **Every `v₂(f) = 1` class has `D_K ≡ 0 (mod 8)`**, and the maximum
  height in a class equals `v₂(f)` in every class with a floor (W5), which
  is the volcano structure read off the data.

## 4. Proof status

A1 is an empirical law with zero exceptions; it is not proved here.  The
derivation in the protocol gives the mechanism: on a weak curve the two
conjugate 2-torsion points `(α, 0)` and `(σα, 0)` are halvable together or
not at all, which ties the curve's 4-torsion to its height, and a depth-1
volcano has no level at which that pattern fits.  Turning that into a
proof, and explaining the `p mod 8` split, is the open item.

## 5. What it changes for the walk, and what it does not

The cover route's walk (ledger §17) starts from a curve of order divisible
by 4, moves by 2- and odd-degree isogenies, and stops at a weak curve or
when its component is exhausted; it succeeded in 20 to 54% of runs,
falling with `p`, in `[q/9, q]` distinct curves when it did.

With this census the walk becomes:

1. **Decide before walking.**  `t` is public.  Compute `v₂(t² − 4q³)` and
   `D_K`; if `v₂(f) = 1`, the class holds no weak curve and the walk is
   refused at zero cost.  Half of all classes, and every one of the walk's
   exhaustion failures in those classes, go away.
2. **Ascend, then walk level.**  Odd-degree isogenies keep the 2-volcano
   height; vertical 2-isogenies change it by one, and the direction is
   decided by halvability tests on the 2-torsion points.  Ascend to
   height 3 (or to `v₂(f)` if that is 2), then walk with odd degrees only.
   Every curve visited is then at the densest level, about `13.5/q`
   instead of the `2/q` of height 1, where a random walk spends 75% of its
   steps.  Expected curves visited fall from about `q/2.7` to about
   `q/13.5`, a factor near 5, and the class is never exhausted without an
   answer.
3. For `p ≡ ±1 (mod 8)` the deeper levels are denser still (`p = 7`: 19 to
   26/q at heights 4 to 6) and the ascent should continue to the crater.

What this does not change: the route's cost.  At `p ≥ 251` the walk is
below 1% of the end-to-end cost (ledger §18.4), so the route's `S / rho`
moves by less than that.  What moves is the walk's **reach**: which
curves of order divisible by 4 the route can reach at all, and whether it
can say so without running.  Nothing here touches a prime-field curve,
and the weak class exists only over an extension field of degree
divisible by 3.

## 6. Files

| file | role |
|:--|:--|
| `PROTOCOL.md` | frozen before the run; Addendum A registered after `p ≤ 11`, Addendum B after `p = 13`, each before the next sizes were read |
| `census.rs` | the instrument; `rustc -O census.rs -o census`; `./census P results`; `./census P --selftest`; `./census tabulate results P…` |
| `results/census_p{p}.jsonl` | one record per class: `t`, counts, `D`, `D_K`, `f`, `v₂(f)`, residues, by-height histogram |
| `results/summary_p{p}.json` | the run's own summary line, with the W3 conflicts listed |
| `results/TABLES.md` | the `tabulate` output for every size |
| `results/run_p{p}.log` | stderr of the larger runs |

Related: ledger §§13, 17, 18 on branch `research/jv-cover-end-to-end-20261007`
(`research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md`),
whose §18.3 registered this question as exploratory.  A parallel session on
this host was running `jv_isogeny_walk --exact-census` and an
`iso1_class_census` at `p = 37` while this census ran; their results were not
read and are not used here.
