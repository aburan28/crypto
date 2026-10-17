# Original C finite-field diagnostic

Raw results and process observations are in `runs.tar.xz`; extract with
`tar -xJf runs.tar.xz` and check `sha256sum -c runs.tar.xz.sha256`.
The AGENTS.md §13a repack at mtime 2026-10-17 00:00:00Z is byte-identical.
Both complete executions are retained; replay-match JSON binds their hashes.
The native comparator removes only elapsed fields. It is not independent replay.

3,738/3,738 prefix targets verified. Actual largest prime order is 7.
There are 15,042 equations from 9,171 queries; multiple equations/query are
allowed. The separate two-target source-model run is not pooled with these.
The four native tests pass, including all nonzero fixture scalar inputs,
all nontrivial destination generators/targets and degenerate-image rejection.

Local certificates in `../local-strata-evidence/` identify 517 vertices across
24 strata, with 1,352 vertices unknown. Native `stratum-statistics.json` keeps
unknown conductor as a null bucket and every status separate. The source
conductor is 967,458,816; the prefix is not the full 1,934,901,248-vertex class.
No full 165-stratum population or qualified runtime claim is made.

[Self-contained Drive archive](https://drive.google.com/file/d/1Xi7mjLCskYEauXlwzAfQQtOoHw6waop-/view)
includes exact binaries, inputs, source/protocol, raw runs, three complete
calibration graphs, native statistics, CSV and inspected four-page PDF.
[Report](../../reports/run-20261017T0626Z-finite-field-C-prefix/report.pdf).
Large-extension finite-field IC, full original graphs, full matrix and exact
representation registration remain open. The method here is projected-prime
subgroup two-summand relations, not NFS or integer-smoothness IC.
