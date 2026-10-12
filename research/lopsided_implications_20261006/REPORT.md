# Wider implications of the lopsided refutation: audit report (2026-10-06)

Questions and method in the [protocol](README.md). Rounds 1–2 built the
thin-product interface, assessed four IC insertion points, and screened
twelve further targets (interface merged in PR #1485; survey in open PR
#1493). This round asks what else the discovery implies. Status labels:
**[proposed]** / **[derived]** / **[verified]** / **[measured]** (none in
this round).

## Headline

**Our hardness foundations do not move.** The repository's floor is the
unconditional generic-group bound (Shoup 1997) and its reference is
measured rho — neither cites the refuted hypotheses. The audit finds zero
load-bearing dependencies on 3SUM/APSP/Exact Triangle hardness, one
research-note sentence naming the conjecture (repaired in this PR), and a
trail of false positives. The method implication is a proposal: future
algebraic-identity claims should ship machine-replayable certificates, as
the paper's Lean verification exemplifies.

## Audit table (receipts: commands below)

| Hit | Context (quoted) | Classification [derived] |
|---|---|---|
| `3SUM over a group`, `3SUM hardness ... generic sets` in `research/notes/ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md` §§8–9 | "Finding P1+P2+P3 = R ... is 3SUM over a group, and pair enumeration realises the generic Theta(\|F\|^2) bound ... 3SUM hardness is a conjecture about *generic* sets" | **Repaired (this PR)**: the sentence named the conjecture itself. Integer 3SUM is now refuted (Alman–VVW), but that algorithm does not transfer to curve-point decomposition, so the curve-point analogue stays conjectural and the note's conclusion (need a structure-exploiting oracle) is unchanged. Edit preserves the conclusion and adds the refutation as a sharpening |
| `generic 3SUM search` in `research/nagao_relations/solver_15/contract.json` | "...an oracle that exploits the algebraic structure of F rather than re-indexing a generic 3SUM search over it" | **Already correct**: demands algebraic structure over generic search; the refutation (integer domain, no transfer) does not weaken the requirement |
| `dense_3sum_wall_bits: 80` in `docs/ic/boundary_targets.json` ("dense 3-sum wall ~80 bits") | Historical wall-clock baseline of our own dense 3-sum enumeration over the factor base (superseded/null-verdict rows) | **Historical measurement, unaffected**: data about our enumeration, not a hardness assumption; measured numbers stand |
| `matrix-vector` in `ic_boundary.rs`, `ic_framework/stages.rs`, `pq_wiedemann.rs`, `FRAMEWORK.md` | Wiedemann/matvec operation counts | **False positive**: sparse matvecs, not the OMv conjecture |
| `b3sum` in a tournament `registration.json` | `source/vendor/blake3/.../build_b3sum.py` (build file manifest) | **False positive**: BLAKE3 build tooling |
| `omV` in `src/utils/random.rs` | `crypto.getRandomValues` (web-API name) | **False positive** |
| `Sethi-Ullman`, `fine-grained token/timing` in cryptanalysis tree | Instruction scheduling; API token scopes; side-channel timing | **False positives**: unrelated meanings of the words |
| `conditional ... bound` on the scoreboard (e.g. "conditional architecture bound") | Component/counting bounds conditioned on stated architecture premises | **Already correct**: "conditional" means premise-relative accounting, not fine-grained-conditional hardness; no language change needed |
| Scoreboard verdict, `docs/ic/FRAMEWORK.md`, ledger floors | "scored against their generic-group floor and matched Pollard rho reference"; floor `S = sqrt(pi/2A)`; Shoup 1997 citation in the residual-walks note | **Unaffected (unconditional)**: generic-group lower bounds do not depend on any refuted hypothesis |

Receipts **[verified]** (2026-10-06, branch `cursor/lopsided-implications-35ee`):
tree-wide case-sensitive grep for `3SUM|APSP` outside the two lopsided
studies returned only manifest/JSON noise plus the two genuine prose hits
above; `Exact Triangle|Zero-Weight` returned one prose hit
(`zero-weight` in a rotated-row certificate note — relation-weight
language, verified below); OMv/dynamic-lower-bound patterns returned only
matvec substrings. Follow-up greps with quoted context separated the
`b3sum`/`getRandomValues`/scheduling/token/timing false friends.

The `zero-weight` hit: `research/notes/ecc2k130/rotated_row_certificate_20260925/RESULT.md`
uses "zero-weight" for relation-row weight classes, not the Zero-Weight
Triangle/k-Clique hypotheses — **false positive** (domain mismatch).

## Dependency graph

<figure>
<object data="figure.svg" type="image/svg+xml" style="width:100%">dependency graph (see figure.svg)</object>
<figcaption>Figure 1 -- what the refutation touches and what it cannot reach.
Green: intact. Red: refuted hypotheses and the one sentence repaired.
Grey: false positives. Editable source: <code>figure.mmd</code>.</figcaption>
</figure>

## Method implications [proposed]

1. **Machine-found identities need machine-checkable certificates.** The
   paper's algorithm was discovered by Claude and its main results were
   verified in Lean 4/Mathlib. The repository already lives this ethos for
   runs (`ecbench verify --replay`, sealed records); the gap is
   *algebraic-identity* claims (bilinear identities, summation-polynomial
   relations, orbit-fold counting arguments). Proposal: any future claim of
   the form "identity I holds" ships a certificate its checker replays
   without trusting the author — a Lean proof, or at minimum a
   randomised-evaluation harness with frozen seeds independent of the
   discovery code — committed beside the claim. Next action: a follow-up PR
   writing the certificate format into the relevant study template; no
   format is imposed here.
2. **Reductions are algorithms — screen them on arrival.** The paper turns
   hardness reductions into speedups by composing them with a new
   primitive. Our `screen_workload` gate (round 2) is the standing mechanism
   for this: any newly published transfer-shaped result gets screened
   before it gets cited. Next action: cite the gate (not this report) in
   future protocols that consider external algorithmic transfers.
3. **Tiny-exponent polynomial gains do not move our verdicts — by design.**
   The paper's speedups (e.g. n^2 to n^1.9992) are real and do not change a
   single `S` row here, because end-to-end measured cost against rho is
   exactly the filter such results must pass. No methodology change is
   needed; the implication is confirmatory.
4. **Watch the general regime, not the concrete one.** The paper's
   epsilon below 0.1204 regime (Section 4) is where a future transfer is
   likeliest to come from; the concrete 1/18 line encoded in the gate is
   the conservative screen. This is already recorded as round-2 open check
   3; no new action beyond keeping it open.

## Boundary, table, ratio (AGENTS.md Secs. 1–8)

- Boundaries: generic floor + matched rho, unchanged and — per this audit —
  unchangeable by the refutation (unconditional vs conditional).
- No table row produced: audit classifications are not DLP answers.
- Falsification target: a tree-wide hit showing a floor, reference, or
  verdict depending on a refuted hypothesis would overturn the headline.
  The receipts above are the search that found none.
- All phases unpriced (`null`); wall time appears nowhere.

## Graphs checked — no change [verified by inspection]

- `docs/index-calculus-scoreboard.html` + `docs/ic/progress-timeline.json`:
  verdict rests on the generic floor; no edit owed.
- `docs/ic-leaderboard.html` / `LEADERBOARD.md` / `leaderboard.json`: no
  measurement landed.
- `docs/browser/data.json`: no new identities of any kind.
- `docs/curves/registry.json`: no curve named; family notation only.
- The one prose repair (RR panel) and the `figure.svg` encoding repair are
  the only non-study file changes, both disclosed here and in the PR body.
- Figures are new visuals for this search round, not canonical-graph
  updates.

## Open checks [proposed]

1. Round-2 PR #1493 is still open; its merge does not affect this audit
   (different files), but its screening gate is cited above — re-verify the
   citation path after it merges.
2. The certificate-format follow-up (method implication 1) is unowned;
   owning it means writing the template section, not just proposing it.
3. If a future paper transfers thin-product speedups to curve-point
   decomposition, the RR-panel sentence repaired here must be revisited —
   the repair names the transfer gap explicitly so the revisit has a hook.

## Sources

- Paper: `https://arxiv.org/pdf/2610.06783` (v1); Fig. 1 (p. 9), Secs.
  5.1–5.4; Sec. 6 (future directions); acknowledgments/methodology on
  machine discovery and Lean verification.
- Rounds 1–2: `research/lopsided_thin_product_20261006/` (merged PR #1485),
  `research/lopsided_other_speedups_20261006/` (open PR #1493).
- Audit receipts: the greps quoted above, run 2026-10-06 on this branch.
