# Registered generic/reference comparison and observer study

The [protocol](PROTOCOL.md) freezes the next qualification before measurement.
No results exist at registration. This is a development qualification, separate
from the two remaining improvement rounds. It cannot promote a challenger or
silently replace the accepted IC/rho references.

`generic_reference_qualification.py` runs the existing tournament with the exact
qualified `pairinv` and `both` sources, generic dense/sparse IC, and six requested
rho widths from each source. It verifies the source bytes before preparation,
uses the same frozen one-point workloads across all arms, preserves every failure,
and audits the complete 1350-slot schedule with the transported evaluator.

Before any measured comparison, `generic_observer.py prepare` freezes 360 adjacent
enabled/legacy native pairs on the fifteen development points. `run` preserves
both outputs, process limits, external audit timing and failure rows. `verify`
reconstructs all certificates, semantic comparisons and descriptive statistics
without executing a worker. Legacy mode remains an unadmitted observation with
an `OBS1` identifier and null candidate ID, scientific phases and instruction cost.
Enabled observations retain the parent's candidate/reference and workload IDs,
with separate fresh execution numbers. Missing legacy matrix and dispatch fields
are not invented. Its outer timing is a diagnostic, distinct from the enabled
worker's fully charged scientific interval.

The `Optimized IC producer admission` workflow has an explicit
`qualify_generic=true` dispatch input, default false. All producer controls must
pass first. Ordinary PR CI runs only the 45-pair observer wiring control on the
existing n13a0 fixture panel, including nine expected preparation failures. The
research job retains the complete comparison/observer data and controlled build
and has a 120-minute cap. A measured failure is never retried or overwritten.

Run locally only in the registered Linux environment, with admitted artifacts
and a controlled generic build:

```sh
python3.12 research/ic_candidate_tournament_20260915/generic_reference_qualification.py \
  --artifacts /tmp/qualified-producer-inputs --generic-build /tmp/ic-generic-build \
  --out /tmp/ic-generic-reference-qualification
```

The output must be new. Reports retain separate online/cold leaders, uncertainty,
actual base sizes, folded columns, full stage evidence and unsuccessful attempts.
The observer study measures whole-mode effects; no overhead is subtracted from
any IC or rho result. Its descriptive intervals are not held-out significance
tests. Registration alone supplies neither measured costs nor a speedup.

After execution, retain raw evidence durably, independently replay it and update
the research note and canonical scoreboard. Exclude all 25 exposed points from
later confirmation through the existing supplemental-fixture mechanism. A new
reference binding and versioned round-two panel must precede competitive use.
