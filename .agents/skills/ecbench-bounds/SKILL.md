---
name: ecbench-bounds
description: Record what an ECDLP method costs as a sealed, re-derivable bound (constant and exponent with intervals, per phase, in ecbench's counted unit, with the field multiplications, squarings and inversions behind the group operations where the curve counts them), read the frontier of methods nobody has beaten, and challenge an incumbent with a frozen paired spec whose verdict says whether a candidate advances, trades, matches or regresses and at which level (exponent, constant or primitive). Use whenever someone asks what the established cost of a method is, whether a change beat it, or how to try.
---

# Bounds, frontiers and challenges

Read [`docs/bounds/README.md`](../../../docs/bounds/README.md) first: the
four records, the four levels, the dominance rule and the acceptance rule are
there. This skill is the procedure. Native `ecbench` only; no Python harness.

## Answer "what does method X cost here?"

1. Read `docs/bounds/FRONTIER.md` for the domain (problem, family, target
   kind, unit, tier). Quote the row's `ops` (× floor) with its interval, its
   memory, its `α` with its interval and the sizes it was measured on. Never
   quote a point without the interval, and never carry a `toy` figure to
   another tier. Where the domain's table has `field muls`, `field sqrs` and
   `field invs` columns (the primitive level, README §2 and §5: modular
   multiplications, squarings and inversions behind the group operations,
   per `√r`), quote them the same way; a row showing `unknown` there did not
   count them, and unknown is not zero.
2. If the method has no record, fit one from a committed session:

   ```bash
   cargo build --release --bin ecbench
   ./target/release/ecbench bound fit --dir research/TOPIC/sessions/S --arm ARM \
       --audit research/TOPIC/audit.json --root . --out docs/bounds/records/NAME.json
   ```

   A session with fewer than four sizes gives a `constant` bound; the exponent
   is then descriptive and you say so. Mixed tiers are refused: pass `--tier`.
   Check it and rebuild the page:

   ```bash
   ./target/release/ecbench bound check --root . --record docs/bounds/records/NAME.json
   ./target/release/ecbench frontier build --bounds docs/bounds/records \
       --out docs/bounds/frontier.json --markdown docs/bounds/FRONTIER.md
   ```

## Answer "did my change beat it?"

Only a verdict answers this. A `compare` of two arms in a session of your own
design is a measurement; it is not a frontier move.

1. Register the candidate as a method (`ecbench-extend`). A changed algorithm
   is a new registry id; never change what a registered id counts.
2. Pick the standing challenge for the domain and axis in
   `docs/bounds/challenges/`, or seal a new one (`challenge seal`) when the
   domain has none. A new challenge names the frontier holder's `bound_id` as
   its incumbent.
3. Pick an epoch nobody has used for that challenge and write the spec:

   ```bash
   ./target/release/ecbench challenge spec --challenge docs/bounds/challenges/C.json \
       --candidate '{"id":"METHOD","params":{...}}' --epoch N --out spec.json
   ```

4. Run it as any session (`ecbench-measure` §3; counts need no isolation
   level), then judge it with every run replayed:

   ```bash
   ./target/release/ecbench challenge verdict --challenge docs/bounds/challenges/C.json \
       --dir OUT_DIR --epoch N --replay-all --bounds docs/bounds/records --root . \
       --out OUT_DIR.verdict.json --bound-out docs/bounds/records/NAME.json \
       --audit-out OUT_DIR.audit.json --exit-code
   ```

5. Read the statement. `advances` names the axes and, when `ops` is among
   them, the level: `exponent` only with disjoint `α` intervals over four or
   more sizes, otherwise `constant`. When `ops` is not among them but a field
   axis is (`field_muls`, `field_sqrs`, `field_invs`), the level is
   `primitive`: the group-operation count held and the field work behind it
   fell. A field axis decides only when the challenge's `acceptance.axes`
   names it (a cheaper formula is a claim on `field_sqrs`, so name it);
   otherwise the three are reported, not deciding, whenever a run of either
   arm counted them (prime-field curves today), and left out when nothing
   counted them and nothing asked. Named but uncounted, an axis is `unknown`
   and decides nothing. `trade` is a new Pareto point, not a replacement.
   `matches` is a null result and is still committed.
   `inadmissible` names why; fix the cause, never the rule. Read `stages` to
   say which sub-algorithm moved, and `accounting` before claiming an advance
   that comes with more unpriced work.
6. Rebuild the frontier and land the session, the verdict, the receipt, the
   bound and the page in one PR (`ecbench-measure` §7). CI re-derives every
   record and fails a stale page.

## What not to do

- Do not edit a record, a verdict or the page by hand; regenerate.
- Do not call a primitive-level change (a cheaper formula) a constant or
  exponent change in `ecbench.gae`; the `gae` figure does not see it (README
  §2, §7). On prime-field curves it shows on the `field_muls`, `field_sqrs`
  and `field_invs` axes at an unchanged `ops`, and that is the `primitive`
  level; on binary and Koblitz curves and for index calculus those axes are
  unknown, never zero, and the change is not measured here at all.
- Do not average verdicts across epochs or pool a `toy` row with a `medium`
  one. One verdict is one session; one domain is one tier.
- Do not report wall time as a bound. It is a practicality note, graded by
  isolation level, elsewhere.
