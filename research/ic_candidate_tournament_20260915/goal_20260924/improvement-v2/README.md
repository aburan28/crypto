# Versioned references for the remaining bounded rounds

The tournament now supports separate qualified IC references for cold and online
metrics through `--campaign-version 2`. This consumes the reference and observer
evidence accepted in [PR 876](https://github.com/aburan28/crypto/pull/876).
The [protocol](PROTOCOL.md) fixes the exact sources, settings, budget and gates.

The cold incumbent and online IC leader differ on the qualification panel.
Pairing every metric against the cold reference would lose the stronger online
comparison. Version 2 pairs cold instructions and cold native time with
`incumbent`, and online time with `ic_online`. Both rho references remain
mandatory. IC reference arms are reported but cannot enter the challenger
portfolio or become a selected challenger.

The accepted observer evidence supports fully charged instrumented comparisons.
Its legacy mode remains diagnostic. No cost is subtracted, and no uninstrumented
or low-overhead performance claim follows from this binding.

Round one retains its own rules and frozen evaluator. Rounds two and three retain
the original eighteen-test familywise budget, 72 fresh confirmation targets and
three repetitions per point. At most ten challengers plus the incumbent fit the
unchanged 3,500-pair cap alongside the three extra reference arms; the maximum
schedule is 3,480 pairs. Six diverse challengers, including an exploration slot,
remain the selection budget.

## Prepare declarations without exposing targets

Restore the accepted archive with the existing hash-checked restorer. For
example, from the repository root:

```sh
python3.12 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-generic-reference-qualification-20260926 \
  --out /tmp/ic-reference-v2-inputs

python3.12 research/ic_candidate_tournament_20260915/reference_registry_v2.py \
  --bundle /tmp/ic-reference-v2-inputs/ic-generic-reference-qualification \
  --out /tmp/ic-reference-v2-inputs/reference-registry.json
```

The declaration command checks the prepared source bytes, generic build/worker
identity, exact qualification report and accepted observer evidence. It prints
the incumbent source, qualification and observer paths and writes the three
reference declarations in canonical order. It executes no worker and generates
no target.

A registered round uses those paths with the existing `tournament.py prepare`
interface and these versioned options:

- `--campaign-version 2 --attempt-number 2` (or the registered third attempt).
- `--reference-registry` naming the generated three-reference JSON.
- `--qualified-report` and `--qualified-observer` naming the retained raw reports.
- The original history, all preceding `--prior-round` archives, and every required
  supplemental `--exposed-fixtures` input.
- The registered candidate panel, seed, cells, 72-target allocation and resources
  in the protocol. Legacy `--rho-source-root` and `--rho-config` overrides are
  rejected for version 2.

Preparation, execution and replay all bind these references. Generic native and
profiled executions retain separate run numbers. A failed or missing reference
blocks a complete comparison; it cannot become a successful-subset estimate.

## Validation and remaining work

Controls exercise the accepted report/build/source binding, changed reference
settings, the complete schedule budget, failed and missing online references,
portfolio exclusion, and the final decision when a reference is incomplete.
Known-cost synthetic cohorts verify the original 50,000-draw familywise rule
against separate metric baselines. A replay of the committed development run
records checks the new online comparison against the observed online leader.
These are software controls, not additional research measurements.

Round two completed under [round2.json](round2.json): eleven arms on source
digest `8582e4ab4b63e98696a0ff00ee296e2902923a2c39c325ab0f3e3950ffbb2b28`, seed
`2026092552`, wrapper [`run_improvement_v2.py`](../../run_improvement_v2.py).
The frozen checker verified 3,480/3,480 pairs and retained the incumbent;
selected `stop5_word` failed promotion. See
[the evidence report](../improvement/round2/EVIDENCE.md) and archive
`ic-improvement-round2-20260928`. Do not redispatch `run_round_two=true`. One
registered attempt remains; a third round needs a new panel, seed and complete
target exclusions, and must not retune on confirmation or replay.
