---
name: ic-propose-candidates
description: Develop concrete elliptic-curve index-calculus candidate changes from measured bottlenecks, with falsifiable hypotheses and matched experiment configurations. Use for the crypto repository's next research round or a requested algorithm improvement.
---

# Propose IC candidates

Locate the crypto repository and read its `AGENTS.md`. The operational code lives
under `research/ic_candidate_tournament_20260915/`; read `OPERATIONS.md` and the
previous round's `decision.json`, stage summaries and exclusive phase costs.

## Candidate work

- Identify which complete-solve phase dominates, and distinguish a mathematical
  mechanism from a configuration or implementation improvement.
- State the applicable boundary before measuring. For a fixed signed base of B
  points and m summands, at most `binomial(B+m-1,m)` target images are possible.
  A faster solver does not remove a coverage barrier. Count Frobenius images once.
- Start from the actual selected source/configuration, including local changes.
  Preserve its snapshot. Keep factor-base support and m fixed for implementation
  comparisons; give support-changing experiments a separately declared panel.
- Describe each candidate with its parent, source/configuration change, hypothesis,
  affected phases, expected memory impact and numeric falsification condition.
  Record regressions and deduplicate configurations already tried on that source.
- First measure individual changes. A combined implementation needs a new full
  run; adding the cheapest phases from different runs is not a measured algorithm.

The current worker's implemented configuration surface is deliberately small:
pair-table/enumeration decomposition, collection batch/window, and dense/sparse
scalar linear algebra. Generate its starting candidates with:

```bash
python3 research/ic_candidate_tournament_20260915/tournament.py propose --out /tmp/ic-candidates.json
```

Inspect the generated registry and make bounded, evidence-backed edits before
`prepare --candidates FILE` freezes it. Each entry has an `id`, `hypothesis` and
`config`; retain the first `incumbent` entry. Unsupported knobs must be implemented
and tested before they become runnable candidates.

For a subsequent round, generate new configurations around the selected incumbent:

```bash
python3 research/ic_candidate_tournament_20260915/tournament.py propose \
  --from-round /absolute/path/to/finished-round --out /tmp/ic-next-candidates.json
```

Use the printed baseline source path, target count, metric gates and a new seed
for preparation. Single-target and multi-target jobs remain separate panels;
reusing logs across targets must still charge the complete job's setup. This avoids
repeating configurations from the preceding round; also inspect older rounds for
duplicates of the same source/mechanism. Stop when the agreed budget or three
rounds without promotion is reached.

For source-level research, use an isolated derivative of the frozen source. Add
`"source_root": "/absolute/path/to/candidate-checkout"` to that candidate's registry
entry. `prepare --source-root BASELINE` freezes and builds the baseline separately,
then freezes and builds each source candidate. Identical configurations are allowed
when their source hashes differ. The exact same checker and scoring rules apply.
Keep instrumentation and the evaluation boundary intact; patch the algorithm,
not the profiler or certificates. Changes outside the worker's supported interface
need an adapter and equivalence checks before performance admission.

Use only development evidence for iteration. Selection data picks a provisional
challenger; final-confirmation data is for one locked challenger. Once exposed,
those final cases must not be reused as fresh confirmation in later rounds.

Deliver concrete runnable candidates or explicitly identified implementation work,
their falsification targets, and the next bounded complete-DLP experiment. Never
describe a proposed change as a measured gain.
