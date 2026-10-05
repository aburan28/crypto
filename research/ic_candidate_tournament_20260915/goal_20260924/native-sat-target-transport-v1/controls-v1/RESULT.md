# Disclosed n17 SAT transport controls

These controls test source-instance and input-transport compatibility on three
already disclosed public points. They are **not** natural-yield observations,
fresh targets, complete IC runs or speed measurements. The source for the two
probes, prepared exporter and bounded prepared-child transport was committed
at `84b5f74069cd4c288326433173a4beab53f60d75` before the first control.
The post-hoc watchdog amendment was committed at
`e56264022` before its separately labeled run. All controls used the same
accepted macOS ARM64 CryptoMiniSat and exporter assets as the scientific
ordinary registration. Build: `cargo build --locked --offline --release`
with Rust 1.93.1 (Homebrew), Darwin 25.6.0 arm64; the Cargo.lock SHA-256 was
`f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79`.

| Control | Public point | Accepted file result | Prepared stdin result | Control outcome |
| --- | --- | --- | --- | --- |
| CMS query-00 v1, 60 s watchdog | [40991, 73355] | TIMEOUT | SOURCE_UNSAT | **FAIL**: status parity absent |
| CMS query-01 v1, 60 s watchdog | [73003, 104622] | SOURCE_UNSAT | SOURCE_UNSAT | PASS |
| CMS query-02 v1, 60 s watchdog | [59775, 2910] | SAT_MODEL | SAT_MODEL | PASS; both ANF/CNF models lift to the same full-point witness |
| CMS query-00 v2, 120 s watchdog | [40991, 73355] | SOURCE_UNSAT | SOURCE_UNSAT | PASS, post-hoc transport diagnostic only |

The first frozen failure remains a failure. Its file run reached about 375,000
conflicts before the 60-second wall deadline; stdin proved UNSAT at 415,862
conflicts. On this unisolated host, those times do not establish a transport
speed difference. The longer v2 run changed only the wall watchdog, retained
the original one-million-conflict cap and exact input bytes, and does not
retroactively pass v1. The four CMS result records linked in the table below
include the failure,
stdout/stderr hashes, binary and source hashes, geometric class, model checks
and precise status parity. The positive model in both modes independently
satisfies the retained ANF and CNF and lifts to factor-base indices [53, 1,
48]. No native child-group drain audit is claimed for the accepted-CMS
compatibility probe.

| Prepared-exporter control | ANF/CNF-XOR/Magma fixed SHA-256 parity | Deterministic manifest/source identity | Malformed/oversize rejection | Never-fed child drain |
| --- | --- | --- | --- | --- |
| query-00 | PASS | PASS | PASS | PASS |
| query-01 | PASS | PASS | PASS | PASS |
| query-02 | PASS | PASS | PASS | PASS |

The exporter controls prestarted the new binary with static field/base
arguments, observed its exact READY marker before delivering each public
point through one bounded JSON stdin request, then compared all three source
files against the accepted exporter and the fixed hashes in
[CONTROL_PLAN.md](../CONTROL_PLAN.md). They separately sent malformed JSON
and an oversized 512-byte request to fresh children and cancelled a never-fed
child, retaining every stdout/stderr file, PID ledger and receipt. The
three exporter result records linked below all report
`transport_control_pass=true`.

The pinned release binary SHA-256 values are
`676c85f7e08b9339462ce8cc1cb079ae476539cb7b0bf9103fda4cf541003dd8`
(new prepared exporter),
`a9005ea75e03d35e251880166f0f69bb59b73b6726ae042205eff7b634b23f5e`
(exporter probe), and
`68e832a8d6c40cab4e70ee4abf8223df960da249a70dc26d87e41771e632c874`
(accepted-CMS stdin probe). The accepted exporter and CMS pins remain
`4da7f5781da1aba2f76d6b9ce8f6e4ea43845f2465cf8edf0612615931d0629a`
and
`6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af`.
The complete seven raw control directories are in
[`raw-controls.tar.gz`](raw-controls.tar.gz), SHA-256
`cc8d409ea2fc8dc796e5584f94f9d3603e732ea762e1ef6d6b9a45827994a9fe`,
with 133 tar entries. Individual result SHA-256 values are:

| Record | SHA-256 |
| --- | --- |
| [exporter-00 v1](ic-exporter-prestart-control-00-v1.json) | `53fa33e8bbf43a5815a1218868ea4cbef1d9729b6af7e1e9aa747f7749191404` |
| [exporter-01 v1](ic-exporter-prestart-control-01-v1.json) | `9cef3131d88d85a3cbbc0abd4f4cb43c7e38c6b19f31cdc110cbad851be6282d` |
| [exporter-02 v1](ic-exporter-prestart-control-02-v1.json) | `d8cd085f960afa408b49dbbe258549e08f0737556c34525f0002476eb4f3fb11` |
| [CMS-00 v1 failure](ic-cms-stdin-control-00-v1.json) | `9617df6140309986613072e05ac2c722ffde54fa1e1acea02a9073c6d454e875` |
| [CMS-01 v1](ic-cms-stdin-control-01-v1.json) | `3920b14ba6c1adbae7a6f9ccbed3580cef4f608857a245725333a3fbef411966` |
| [CMS-02 v1](ic-cms-stdin-control-02-v1.json) | `324c65fab8483b438df78d65bc0013d9181f812f5badecfa2ab817d4bbf3a36a` |
| [CMS-00 v2 diagnostic](ic-cms-stdin-control-00-v2.json) | `f7329988d7fd8e40fa7ec4a14137aaca179946270daa672a932cfa08ffa5bc6b` |

The accepted CMS reader line occurs **before** stdin stream and parser
initialization. The separately proposed post-parser marker patch is still
unbuilt. Thus the accepted-stdin compatibility controls do not establish the
target-independent startup boundary required for primary online timing.
Moreover, the independently audited scientific CMS ordinary panel has only
rank 22/29 under its registered 100,000-conflict cap. It supplies no complete
factor-log table. The next SAT preparation needs a new budget or encoding,
source freeze, fixed natural panel and original independent audit before a
SAT target registration can be considered. No SAT online speedup is claimed.
