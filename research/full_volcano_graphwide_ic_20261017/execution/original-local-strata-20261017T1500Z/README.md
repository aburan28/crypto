# Local conductor certification of the original five prefixes

`local-strata certify PREFIX_DB FRESH_CERTIFICATE` now exposes the existing
conservative labeling calculation for any compatible replayed ordinary
prefix. It does not change the formula, databases or global stratum status.
The new command passes all fifteen complete-graph visibility-cut checks;
generic C output equals both the calibration output and previous frozen
C certificate byte for byte. Full source/field/coefficient identity, DB SHA
and existing map-replay metadata stay in each native certificate.

| Historical input | Source End(E) conductor | Exact prefix vertices | Observed strata | Unknown vertices |
| --- | ---: | ---: | ---: | ---: |
| A | 268435456 | 512 | 10 | 1026 |
| B | 129140163 | 729 | 7 | 887 |
| C | 967458816 | 517 | 24 | 1352 |
| D | 544195584 | 509 | 24 | 1372 |
| E | 644972544 | 509 | 24 | 1372 |

Full canonical source identities/inputs and native populations are in
`summary.json`; the names A–E remain historical fixture filenames, not new
curve aliases. [Download prefixes, complete certificates and exact binary](https://drive.google.com/file/d/1zeGSq3GWne9yTH0R37EZ3YI89xG_yR2L/view).
[Inspected report](../../reports/run-20261017T1500Z-original-local-strata/report.pdf).
Extraction/hashes and extracted-binary C reproduction passed. The original
databases are immutable. This does not establish complete A/C/D/E graphs,
all original strata, a new DLP timing result or independent arithmetic proof.
No admitted best-method graph changes: this is input certification only.
