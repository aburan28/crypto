# n83 F6 flat x-index gate

Registered before changing the lookup implementation or running the candidate.

## Hypothesis and scope

The exact four-summand `F6SignedPairIndex` at n83 hashes 4,108,723 signed
pair-sum representatives through a standard `HashMap<u128, usize>`. A flat
open-address table of 32-bit indices, with the full 83-bit x coordinate
checked against the existing sum array at every probe, may reduce random
memory traffic and improve the target query. Keep the first representative
per x, the separate identity slot, insertion and query order, sign tests,
portable fallback, and independent group replay unchanged. The table must
reject indices that cannot fit a `u32`. Linear probing must terminate before
the table capacity is exhausted. No probabilistic key truncation is allowed.

This is an **exploratory four-summand component experiment**, not a complete
F6 or IC speedup. At the frozen K0 base, the uniform-target four-summand
coverage ceiling is `4.662e-12`. Query misses therefore do not estimate
natural relation yield. Reusable index construction is reported separately
from the target-dependent query. Complete F4/F5/F6, one-target IC online and
paired rho ratios remain unknown.

## Frozen inputs and reference

- Parent branch: `codex/f6-n83-compact-xmap-20261005` at
  `558cfc1dc23d91c4a38ed5d5a3051d78e045e1da` (#1421).
- Baseline geometry source SHA-256:
  `feb569917be8b526eee743e988bc01850c3bba32a312f272bd50e42313c3df25`.
- Frozen baseline full, small and planted binaries SHA-256:
  `75009df1bfe412f00390b1d7e707ba5d4cb64d196a143f7888e6bfe235a15d95`,
  `8d2bed511af624b35e5bd48d16b9b42c1920ddaea6ef6163e0c01867a2b6a0e4`,
  `4dda20b3331d799036c48bc8362f78ff773c3b78aadfb48c090e48a4ab7987b6`.
- Curve: registered `icv1-f2m83-tm6151469093347-debefd74`, exact K0
  subgroup. Factor base: standard dimension 12, cofactor projected,
  4,054 actual usable points, 2,027 signed columns, 8,219,485 unordered
  pairs, 4,108,723 signed sum representatives. The smaller dimension
  8/10 probes use 258/1,048 usable points.
- Frozen public T001 target:
  `(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`.
- Build candidate with the same local release compiler and flags as baseline.
  Hash the source, binaries, scripts and raw output. No build during paired
  timing. Record physical CPU, OS, Rust version, and `RAYON_NUM_THREADS=1`.

## Exactness and decision

Run the nine focused `f6_wide_geometry` tests, including exhaustive small-base
membership, identity/same-x exceptions and n83 packed arithmetic. Run the
full-base planted `[0,2,4,6]` control through portable and PMULL paths and
replay the returned witness in the group. For the small and full ordinary
probes, require identical representative counts, exact outcome and group
verification between arms. Preserve failed processes, timeouts and OOMs.

Run a full-base A/B/B/A/A/B/B/A/A/B panel (A = frozen #1421 baseline,
B = candidate), with five processes per arm, a 120-second limit each, and
no build between runs. This gives five adjacent pairs, alternating which
arm runs first, plus repeated A/A and B/B controls. Report each process's
index build, exact query, peak RSS and status. A small-base A/B/B/A screen
checks dimensions 8 and 10 before the full panel. Retain the candidate only
if exactness passes; full-base median target-query time falls at least 10%;
peak RSS falls at least 15%; median index-build time rises no more than 10%;
and neither smaller case regresses more than 20% in its two-process median.
Otherwise restore the parent runtime code and publish the negative receipt.

This Mac has no auditable exclusive CPU partition. Any wall ratios are
exploratory, and a noisy or inconclusive result does not promote a speed
claim. The cost boundary remains the complete one-target verified IC online
interval and same-point rho reference under equivalent resources; both are
unmeasured for this component and recorded as unknown.
