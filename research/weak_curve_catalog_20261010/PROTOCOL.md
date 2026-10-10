# Weak-curve catalogue: scope and evidence contract

Recorded before catalogue generation on 10 October 2026.

The requested deliverable is a growing classification of weak elliptic-curve
families, their properties, and ways to recognize and construct their members.
This round supplies a canonical family JSON catalogue, a native validated
index of existing exact synthetic models, a searchable local browser, a
source-linked report, a classification diagram, and a reviewed PDF. It extends
the existing ISO-1 publication PR without changing its frozen measurements.

Each family records the field and subgroup domain, the mechanism, exact or
conditional recognition, model versus isogeny-class scope, construction or
finite enumeration description, necessary verification, source attribution,
overlap, and current implementation/evidence status. A structural signal is
labelled separately from a reduction or geometric covering witness. Unknown
subgroups, generators, complete endomorphism orders and costs stay null.

The initial instance index imports the independently replayed models in
`two_branch_20261009/curve_records.json` and the source models in the frozen
512-source `large_population_20261009/population_run1/curve_records.json`.
The complete p7 two-torus census and per-trace CSV supply exact hyperelliptic
model and class labels only in their recorded field representation and domain.
The later 512-source two-branch search supplies bounded status and costs, with
unresolved class labels preserved. Larger odd-degree constructed controls
retain geometric witness status; their solver advantage remains unmeasured.

The native builder binds every input by SHA-256, checks exact model identity,
field representation, trace and cardinality, and rejects duplicate identities
with inconsistent records. Direct ordinary p7 membership uses the complete
parameter invariant sets; class support uses the complete per-trace oracle.
An exact class zero is accepted only after that oracle join. Other-family
membership remains unclassified until its recognizer is implemented and
independently validated. Existing cubic exclusions remain branch-scoped.

Counts in the catalogue are inventory counts. The p7 replay inventory includes
successful routes and conversions, and is postselected. The 512 large sources
use the original independent model law. Neither inventory is silently changed
into a uniform isogeny-class density estimate. Family overlaps are retained;
their counts are not added as disjoint populations.

Validation covers catalogue references and scopes, exact source bindings,
model identity, positive/negative/unknown label joins, mutation rejection for
unsafe status promotion, deterministic generation, and browser/PDF layout.
All arithmetic and research data processing use native Rust and existing
compiled PARI evidence. No new timing ratio is promoted in this round.

The canonical IC ratio graphs and the prior two-branch census figures retain
their measured values: this classification round adds organization and source
links, rather than a new paired cost measurement. The catalogue gets its own
family relation diagram and inventory summary. Requested broader classification
and efficient large-field class discovery remain open where no proved test or
verified construction cost is available.
