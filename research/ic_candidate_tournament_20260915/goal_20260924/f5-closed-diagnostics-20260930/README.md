# Closed F5 diagnostic outcome

Status: retrospective diagnosis complete; complete F4/F5 solver admission and
fresh paired qualification remain pending. No producer or native solver ran.
The accepted source archive and every registration remain unchanged.

The measured row being inspected is
`(IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h57ee76cbeaed, 4937b2d124b2,
IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h57ee76cbeaedW4937b2d124b2R0)`.
Its curve is `EC1N17Ckb1hbbe2b5b6b1e6`: field degree 17, subgroup order
65587, cofactor two. The independent replay confirms 63 geometric points,
62 distinct nonidentity cofactor images, 29 folded columns and final rank
28. The target stage never began; verified online time and speedup remain null.
This small kb1 development instance is not an ECC2K-130 fidelity gate.

## Retained counters

Every one of the 256 ordinary attempts already retains all nine SolveStats
fields. In particular, reductions, splits and propagations do not need to be
reconstructed from the human effort label. The numeric summary below describes
statuses within the same historical candidate, not a new candidate comparison.
All counts are unconverted native diagnostics; none is a calibrated operation
unit or a measured online speed.

| Original outcome | Attempts | At the 4096-reduction cap | Reductions per attempt | Splits per attempt | Actual built degree | Oversize events |
| --- | ---: | ---: | --- | --- | --- | ---: |
| Witness | 61 | 0 | 14–4010 | 10–1543 | 3 in every attempt | 0 |
| Bounded incomplete | 195 | 195 | 4096 in every attempt | 1570–1579 | 3 in every attempt | 0 |

The incomplete attempts report exhaustion, no unsupported encoding and no
matrix oversize event. This identifies the observed stopping condition. It
does not identify why that traversal missed a particular witness. Missing
per-query matrix word counts, cache hit/miss counts and observed resolved
variable/value policy remain instrumentation gaps. Reuse the existing core
counters; add the missing observations without changing an old manifest.

## Mathematical chain controls

The auditor independently re-adds every retained full-point witness: 61 F5
witness records and 57 SAT v2 records on the original 189-query prefix. There
are 118 records, including queries witnessed by both arms; these are not 118
independent trials or a pooled natural-yield sample. All six summand orderings
for every record have a finite intermediate point and satisfy both mathematical
S3 links, giving 708 finite ordering controls. No identity intermediate occurs
in this retained witness set. Separate disclosed unit controls cover repeated
points, inverse signs, two-torsion and identity intermediates.

The already exposed missing-direction witnesses reproduce:

| Ordinary query | Probe scalar | SAT geometric indices | F5 reductions / splits | Independent coefficient in column 27 | Mathematical finite chain |
| --- | ---: | --- | --- | ---: | --- |
| 164 | 8015 | 62, 9, 43 | 4096 / 1573 | 51777 | All six orderings |
| 173 | 26268 | 2, 44, 51 | 4096 / 1575 | 13810 | All six orderings |

Column 27 is independently reconstructed as `[18941,32935]`; every retained
F5 matrix row has coefficient zero there. These two witnesses therefore have
no obstruction at the mathematical finite-chain level. The implemented Boolean
ANF, parameter specialization, variable permutation, F5 criterion and matrix
consequences are still unvalidated by this analysis. A valid mathematical
chain alone does not prove the encoder faithfully represents it. No budget or
branch policy is selected, and the F5 family is not promoted.

The next gate is to export the exact implemented Boolean systems for these
disclosed inputs, evaluate their full-point-derived chain assignments, and
cross-check F5 consequences against unfiltered F4 on identical systems. Then
freeze a bounded development-policy study before a new complete-DLP control.
Fresh qualification and calibrated rho comparison follow those gates; the
three sealed confirmation rounds stay closed.

## Reproducible retained analysis

The protocol was committed at `71e17a71e` before analysis. The initial auditor
was committed at `7d05047ff`; the corrected presentation at `724cf49fc`.
`RESULT-v2.json` is the current derived receipt and includes every per-attempt
counter row, every witness ordering, exact input/source hashes and invocation.
Its SHA-256 is
`1dda1cf39581e24e2b9fadc9e65e92ed97d4d011ab019bf40e710113fa0b6cd3`.
The 28,948,695-byte input archive remains
`039e535fa73e5c794aeff666bc02989076f04abbf751b5f609884a2785760678`,
retained in the accepted paired-evidence directory alongside its full inventory.
All archive members were checked before using the records.

The first `RESULT.json`, stdout and stderr are retained. Its generic integer
summary added built degrees, which have no additive work interpretation.
The corrected auditor reports built-degree histograms and ranges and excludes
degree from additive totals. `PRESENTATION-CORRECTION.json` hashes both receipts
and verifies that all counters, witnesses, ranks, identities and other results
are unchanged; only that presentation and its source/invocation metadata differ.
This is an accounting correction, with no changed solver run or new hypothesis.

Five focused controls pass, including malformed counter/dispatch rejection,
separate reduction and split counts, repeated/sign/identity cases, a two-torsion
case with no finite intermediate, false full-point witness rejection and an
external archive mismatch rejected before member reads. No timing or kernel
benchmark is inferred from these tests.

To replay the derivation into a new output without invoking a solver:

```sh
env PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/f5_closed_diagnostics.py \
  --bundle research/ic_candidate_tournament_20260915/goal_20260924/paired-fresh-n17a1/results-20260929 \
  --out /absolute/new-f5-diagnostic-replay.json
```

Input hashes and mathematical/counter results must match. Invocation paths and
Python executable metadata are replay provenance and can differ on another
host. Post-execution source/import inventories are analysis provenance, not a
retroactive repair of any historical measured candidate's binding.
