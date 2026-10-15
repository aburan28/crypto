# Native CPU campaign 38766895864

The four ZIPs are unchanged GitHub Actions artifacts from
[run 38766895864](https://github.com/aburan28/crypto/actions/runs/38766895864),
attempt 1, candidate `a21306ca20418a4c4bfcc8448de9b8851b8449ed`, baseline
`79006eb85a400474e83f6a76fac0dca6b60a8b5d`. Artifact IDs: macOS ARM64
11989194184; macOS Intel 11991802290; Linux x86-64 11990867468; Linux ARM64
11990337774. Each ZIP retains executables, source, compiler receipt, complete
raw matrices, campaign attempts, stderr/stdout and isolation receipts. Linux
archives also retain exploratory legacy-family baselines.

`INDEPENDENT_AUDITS.tar.gz` contains native rechecks and the complete AMD
analysis, separately from the originals. Extract each ZIP into its platform
named directory, then extract the audit archive into their common parent.
`FULL_ARTIFACT_FILES.json` inventories that expanded tree, the original ZIPs
and ARM job log with bytes and SHA256. `FILES.json` inventories this committed
packet. Foreign executables are hashed, not executed in the independent audit.

All four 96-case coefficient comparisons pass. AMD EPYC 7763 has all twenty
valid isolated runs and all 640 paired medians independently recomputed. Apple
M1 (Virtual) and Intel i7-8700B each have twenty valid count/output runs and
zero admitted wall samples. ARM completed nineteen admitted runs after the
first A/A command was refused before execution. Its full campaign remains
CENSORED. The native audit retains its passing full matrix but rejects the
incomplete series. No partial ARM timing ratio is promoted. The overall
workflow conclusion is failure. A wrong-source negative control in this
packet also confirms that the native auditor rejects mismatched provenance.

The first audit-archive packaging attempt reported temporary-file errors
collecting macOS extended attributes on the full source volume. That attempt
remains on the build SSD. The committed archive was made afresh using
`COPYFILE_DISABLE=1`, `--no-xattrs` and `TMPDIR` on the build SSD. A fresh
extraction compares byte for byte with all four original audit directories.
Every original ZIP passes `unzip -t` and matches the downloaded copy exactly.
No raw receipt was deleted or replaced.

Read the [source-linked report](../cpu-report-20261015/REPORT.md),
[reviewed PDF](../cpu-report-20261015/REPORT.pdf) and
[all 64 cells](../cpu-report-20261015/RESULTS.csv) for the exact boundary,
A/A noise ranges, slower cases and VM limitations. Native source:
`../audit_cpu.rs`, `../summarize_cpu.rs`, `../make_cpu_report.rs`.
