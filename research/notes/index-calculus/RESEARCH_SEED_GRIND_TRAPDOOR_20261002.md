# Seed-grind rigidity defeat for the GHS magic-number trapdoor

**Class:** preregistered protocol. No run, no variant, no scoreboard row yet
(all measurements below are marked `PENDING`; this PR is the plan, §"What
this owes the scoreboard").
**Status:** experiments not started. This document is the hypothesis, the
boundary, the frozen inputs, the success/stop conditions, and the cost
accounting, written *before* measuring, as the repository requires.
**Adopted weakness oracle (do not re-derive):**
[`src/cryptanalysis/ec_trapdoor.rs`](../../../src/cryptanalysis/ec_trapdoor.rs)
(Teske 2006: `magic_number_full`, `ghs_genus`, `audit_curve`,
`construct_trapdoor_curve`) and
[`src/cryptanalysis/ghs_descent.rs`](../../../src/cryptanalysis/ghs_descent.rs)
(`solve_via_descent_m1`, the recovery path), demoed by
[`examples/ghs_attack_demo.rs`](../../../examples/ghs_attack_demo.rs).
**Implementation language:** Rust, extending the two modules above. No Python
in any algorithm, sampler, auditor, or result generator (repository rule).

> **The question, in one line.** The Teske trapdoor is a binary curve
> `E/F_{2^N}` that looks generic but admits a cheap GHS/Weil-descent known to
> whoever chose the right factorisation `N = n·l`. The public auditor
> `audit_curve` already *catches* it by trying every factorisation. So the
> only interesting trapdoor left is a **hidden** one: a curve whose weakness
> survives a *published verifiably-random generation*, because the generator
> gave its author enough seed freedom to grind toward the weakness while the
> output still passes a weaker, realistic rigidity check. This protocol asks
> whether that is possible at composite `N`, and prices the entropy that makes
> it either possible (an attacker's advance) or detectable (the defender's
> win). **The deliverable is the detector and its false-negative bound, not an
> attack.**

## 0. Why both directions collapse into one experiment

The user's two directions are the two halves of one mechanism:

- **Direction 2 — the hidden weak-field structure.** Already built. The
  weakness is a small Hess magic number `m = dim_{F_2} span{σ^i(√b)}` for some
  factorisation `N = n·l`, giving GHS genus `g = 2^{m−1}` and transporting
  ECDLP on `E(F_{2^N})` (size `~2^N`) to HCDLP on a genus-`g` Jacobian over
  `F_{2^l}`. We adopt it; we do not re-derive it.
- **Direction 1 — hiding it in a rigid-looking seed, and detecting it.** Novel.
  `construct_trapdoor_curve` chooses `(a,b)` directly, so its output is
  trivially audit-visible. The open question is whether an author who may only
  publish `curve = H(seed)` for a public hash `H` can still land on a
  trapdoored curve, and how the auditor's completeness and the author's seed
  freedom trade off.

These are the same curve viewed by constructor and by auditor. One experiment
measures both.

## 1. State the boundary before measuring (§"boundary" rule)

Two boundaries, each derived, neither tuned.

### 1a. Detectability floor (a counting floor)

For a binary curve drawn by a verifiably-random process, let
`p(g*) = Pr[ min over factorisations N=n·l of the GHS genus ≤ g* ]`
under the honest (uniform-seed) distribution. A floor follows directly:

```
  grind_bits(g*) := log2(1 / p(g*))
```

is the expected seed-grinding entropy an author must spend to land a curve of
target genus `≤ g*`. If the published generator exposes `freedom_bits` of
author-chosen seed (e.g. an unexplained constant, a nonce, a "random" counter
start), then the trapdoor is **hideable** exactly when

```
  grind_bits(g*) ≤ freedom_bits          (attacker can hide it)
  grind_bits(g*) >  freedom_bits          (defender flags it by rarity)
```

`grind_bits` moves only with `N`'s factorisation lattice and the target genus
`g*`; it cannot be tuned away by any choice of seed. This is the floor the
thread tracks a **ratio** against: `grind_bits(g*) / freedom_bits`.

### 1b. Attack reference (the usual `S`)

A hidden weak curve is only a weakness if the descent it enables is actually
cheaper than the generic attack on the *same* curve. In the repository unit

```
  S = total operations / sqrt(n)              (n = prime subgroup order)
```

the reference is matched Pollard rho on `E(F_{2^N})` itself, flat at
`S ≈ 1.3`. The descended attack's `S` is the **whole** pipeline, cold:
factorisation search is free to the trapdoor owner but priced for the auditor,
then `descend_m1`/`descend_m2`, HCDLP on the Jacobian, lift back, and
`[k]P = Q` verification. A curve is a real trapdoor only when its descended
`S` sits below rho's, robustly, with a verified recovered logarithm.

## 2. One table, one unit (§"one table" rule)

Every variant is a row; `S` is over the whole method; the floor is a constant
column. Rows to be filled by the follow-on measurement PRs:

| variant | min GHS genus | descended `S` | `S`/rho (≈1.3) | grind_bits/freedom_bits | recovers `[k]P=Q`? | class |
|:--|:--|:--|:--|:--|:--|:--|
| honest baseline curve (no grind) | PENDING | PENDING | PENDING | — | PENDING | — |
| Teske construct (audit-visible) | PENDING | PENDING | PENDING | 0 (direct choice) | PENDING | — |
| seed-ground candidate (audit-hidden attempt) | PENDING | PENDING | PENDING | PENDING | PENDING | — |
| matched Pollard rho reference | n/a | ≈1.3 | 1.0 | — | PENDING | reference |

No row is a result until its correctness column is a verified yes.

## 3. Progress is the ratio, not the constant (§"classes")

- **advance** — a curve passes the published rigidity check yet has descended
  `S` below both the honest prior and rho, with `grind_bits ≤ freedom_bits`.
  This is the only outcome that would mean a hideable trapdoor exists.
- **engineering** — the auditor or sampler got faster; detection unchanged.
- **relabelling** — genus pushed down but descended `S` rose (cost paid in
  HCDLP or lift). The §"relabelling" failure mode, drawn in this setting.
- **accounting** — a correction to the counting or the `S` budget; claim no gain.

## 4. Falsification target, declared in advance (§"target")

Fix the composite field first (see §6: primary `N = 51 = 3·17`).

**Trapdoor-exists (attacker) succeeds iff** there is a curve with
`min-genus ≤ g*` such that (i) descent + HCDLP recovers the planted logarithm,
verified `[k]P = Q`, on every seed; (ii) descended `S` < matched rho `S`; and
(iii) `grind_bits(g*) ≤ freedom_bits` for a *plausible published generator*
(freedom model fixed in the protocol, not chosen after seeing the curve).

**Detector (defender) succeeds / we abandon the attacker hypothesis iff** for
every reachable `g*`, either descent recovery fails verification, or
`grind_bits(g*) > freedom_bits` so the curve is flagged by rarity under the
honest prior at a stated false-positive rate.

Inadmissible, as in the reference note: changing `N` or the freedom model after
the fact, changing the operation accounting, skipping `[k]P=Q` verification,
hand-picking favourable seeds, or quoting the HCDLP phase alone as the method `S`.

## 5. Separate the phases, and price all of them (§"phases")

Priced independently, with the dominant one named at each `N`:

1. **Honest-prior sampling** — empirical `p(g*)`: draw uniform-seed curves,
   run `audit_curve`, histogram the min genus. Fixes the floor of §1a.
2. **Grind cost** — seeds tried until `min-genus ≤ g*`; compare to
   `1/p(g*)`. This is `grind_bits`.
3. **Descended attack `S`** — end to end on one trapdoored curve, cold, all
   phases in §1b inside it. This is the only number that earns "weak".
4. **Detector operating point** — false-negative and false-positive rates of
   the rarity test against a declared freedom model.

An exponent/`S` claim covers the whole method or is not made. Fit over ≥4
field sizes where feasible before any extrapolation, and mark extrapolations.

## 6. Scope and subfield disclosure (AGENTS.md §8b — mandatory)

**The GHS magic-number trapdoor is a composite-degree phenomenon.** It needs a
proper intermediate subfield `F_{2^l} ⊊ F_{2^N}`, i.e. `N = n·l` with `l>1`.

- **It cannot exist at the ECC2K-130 challenge field `m = 131`** (prime: no
  proper intermediate subfield over `GF(2)`), nor at `m = 31, 53, 83` (all
  prime). There is no factorisation to descend through. For those fields the
  "detector" for this trapdoor class is the structural fact itself: `m` prime
  ⇒ the Teske/GHS construction is excluded. State it as a **structural
  (partly negative) defensive result for the challenge family**, not an attack.
- **Primary experiment runs at the in-family composite `N = 51 = 3·17`**
  (proper subfields `GF(2^3)`, `GF(2^17)`), the composite-degree comparison
  size named in §8b. Smaller composite degrees are smoke tests only.
- Any finding at `N = 51` is **labelled subfield-dependent** and does **not**
  transfer to `m = 31, 53, 83, 131` without a separate argument. A gain that
  uses `51`'s subfields supports nothing at the prime degrees (§8b rule).
- The prime-degree trapdoor question is **different** and is named here as a
  separate follow-on, not conflated: (a) isogeny-linked weakness, already
  partly audited elsewhere in the tree; (b) seed-grinding toward an
  index-calculus-favourable factor-base/decomposition structure in the
  repository's own IC attacks, priced in the same `S` unit. Each gets its own
  protocol and PR before any measurement.
- Manifests keep field degree `N`, subgroup order `n`, genus `g`, and
  polynomial-variable counts separate; matching bit counts is not matching
  instances.

## 7. Frozen inputs and the commands that will produce the evidence

Adopted, unchanged (the weakness oracle and recovery):
`ec_trapdoor::{FieldTower, magic_number_full, ghs_genus, audit_curve,
construct_trapdoor_curve}`; `ghs_descent::{descend_m1, descend_m2_affine,
solve_via_descent_m1, brute_force_ecdlp}` as the cross-check on recovery.

To be written in Rust in the follow-on measurement PRs (marked PENDING here):

- `src/cryptanalysis/seed_grind_trapdoor.rs` — (i) a verifiably-random curve
  generator `seed → (a,b)` via a public hash, with an explicit `freedom_bits`
  model; (ii) the honest-prior sampler calling `audit_curve`; (iii) the
  grinder that searches seeds for `min-genus ≤ g*`; (iv) the rarity detector.
  All deterministic from a recorded seed; all outputs replayable.
- A panel example + frozen JSON under
  `research/notes/index-calculus/seed_grind_trapdoor_<date>_run/`, scored by a
  Rust/CLI summariser (no Python), each run in its own directory, never
  overwritten; failures/timeouts retained.

Every concrete curve a measurement PR constructs will be registered as an ICV1
slug (`scripts/build_curve_registry.py`, and `docs/curves/sources/specs.txt`
for seed-built curves) in the PR that first names it. This protocol names only
curve *families* parametrically, which §11 permits.

## 8. What this owes the scoreboard (§7 of AGENTS.md)

Nothing yet, and that is stated, not hidden. This is a protocol: no measured
number exists, so `docs/index-calculus-scoreboard.html` gets **no new row in
this PR**. When the §5 phases are measured, their PR adds the rows of §2 in the
same unit and against the same boundaries, carries the class chip per §3, and
cites the frozen run files — the page cites, never computes. Until then the
challenge-field verdict stands unchanged; this protocol does not touch it.

## 9. Cost accounting

- Phases 1–2 (`N=51`): curve sampling + `audit_curve` over the factor lattice
  of 51, cheap; bounded by the seed count, which is itself the measurement.
- Phase 3: one genus-`g` HCDLP at `N=51`; small, intended to *fit on one host*
  as a correctness-bearing end-to-end `S`, not a record.
- Phase 4: detector sweep over the honest-prior sample; cheap.
- No GPU/EC2 needed for `N=51`. If a larger composite is added for
  distributional stability, launch under the `meow34` key pair (§9) and record
  the host manifest and isolation per §10; wall time stays secondary to the
  counted units.

## 10. Reading guide

This note is to the trapdoor thread what
[`RESEARCH_RESIDUAL_WALKS.md`](RESEARCH_RESIDUAL_WALKS.md) is to the
residual-walk thread: boundary and target before measurement, one unit, one
table, and an explicit account of the outcome that would count as a finding
versus the (more likely, and still valuable) defensive negative — that at the
prime challenge degrees the class is excluded by structure, and at composite
degrees it is hideable only when the published generator hands the author more
seed freedom than the weakness costs to grind.
