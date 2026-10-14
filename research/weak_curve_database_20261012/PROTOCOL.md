# Database expansion protocol, frozen 12 October 2026 UTC

## Requested objective and delivered slice

The owner asked to build toward the most comprehensive database of weak
curves, devise theorems and heuristics, and give a clear mechanism for
constructing a weak curve. This release integrates the full local
family/model/property corpus and every complete exact two-branch CSV at
p=7,11,13,17,19,29,37. It also imports the earlier p41,p43,p47 cubic
diagnostics with their weaker label grade. Published mechanisms and
research hypotheses become queryable registries. The cubic norm-one
parameter identity is formalized in Lean.

| Requirement | Current evidence | Status |
| --- | --- | --- |
| Queryable, identity-safe database | SQLite artifact, fixed-source Rust importer, schema, source lock and full reimport check | Verified for the named local corpus |
| Every available complete two-branch class census | Seven CSVs; 87,486 complete exact ordinary rows | Verified import |
| Larger-prime census evidence | p41,p43,p47 CSVs; 252,254 Hasse candidates, 246,512 ordinary rows | Verified import as probabilistic cubic diagnostics |
| Model and property inventory | 814 ICV1 models, 814 matched property rows, arithmetic/hash checks | Verified import |
| Theorems and heuristics | Ten theorem records, four heuristic records, eight prior hypotheses and 21 family proof obligations | Registry complete for current scoped entries; global coverage open |
| A construction mechanism | Cubic norm-one algorithm, proof, source and model gates, exact p7 example | Written and linked to recorded geometry |
| Faithful Lean theorem | Six cubic cocycle and fourth-power declarations, pinned build, axiom audit, fresh replay | Algebraic core checked; finite-field instance and full covering theorem open |
| Visual report | Source-linked report, quantitative SVG/TSV and compiled PDF | Validation recorded separately |
| Largest-database comparison | Public catalogues need normalized model/class/family denominators | Open |
| Exhaustive taxonomy of all published weak-cover cases | Iijima–Momose–Chao tables and later odd-prime cases need normalized import | Open |

## Import rules

Every row has an exact source path and SHA-256. Class rows also retain
their original line number. The importer checks the complete domain of
each seven exact CSVs, and all inclusive Hasse candidates in each
diagnostic CSV. It rejects duplicate traces, malformed 0/1 labels,
incorrect branch unions, impossible trace/order arithmetic, changed
ICV1/model hashes and missing property counterparts. Parent p7 oracle
flags must match the exact p7 CSV. No diagnostic zero becomes a
complete negative class label.

Family entries, individual models, isogeny-class support, covering
geometry, selected subgroup transport and measured cost stay distinct.
An unknown subgroup is stored as null. A published construction can be
a recipe while its local recognizer or covering witness remains pending.

## Expansion sequence

1. Normalize each published equation/table case against the exact
   family definitions and record the source theorem, characteristic,
   degree, genus, covering group, twist convention, exclusions and
   equivalence relation. Add a case only after overlap review.
2. Add an independently checked model witness or a complete
   finite-domain oracle. Retain no-model, timeout, probabilistic
   counter and infrastructure statuses separately.
3. Prove each new mathematical implication in Lean when a faithful
   statement is feasible. Map the declaration to the exact theorem;
   keep geometric and CM obligations open until defined and proved.
4. Extend prospective exact two-branch oracles to p71,p73,p79 only
   after an artifact-registry check and frozen validation split.
   Existing historical p53,p59,p61,p67 data are retrospective.
5. Reuse the 192–252-bit 512-model source panel only under its
   original full-rational-2-torsion law. Select at most 16 source IDs
   per available degree 6,10,14 by hash order before opening outcomes.
   Odd-prime mechanisms require a new independent input law.
6. A selected-subgroup cost claim requires r,G, map transport,
   independently checked target result, exclusive phase costs and
   a same-point one-target rho comparison under isolated CPU evidence.

The database build is metadata integration; it launches no new curve
census, Sage job, pairing computation or target solve.
