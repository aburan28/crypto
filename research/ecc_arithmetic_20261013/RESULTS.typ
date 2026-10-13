#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)
= ECC arithmetic instruction-set pass: result

*Status:* field-level candidate measured; no qualified wall-time result.
The run report is at research/ecc_arithmetic_20261013/reports/run-20261013T1837Z/report.md.
The raw paired data are at research/ecc_arithmetic_20261013/paired_mac_final.csv.

== Implementation

The ARM PMULL and Intel PCLMUL paths now use three carryless products for
two-limb elements and six for three-limb elements. The prior schoolbook
path used four and nine. The four-limb and portable algorithms, reduction,
public types and curve formulas remain the same. These are source-level
operation counts, not a point or scalar performance result.
The existing portable binary-field fallback has operand-dependent control
flow and table access; this change does not validate a constant-time
secret-scalar path.

The modified public multiply passed 17 focused release tests with reference
equivalence and boundary fields. Both frozen sources compiled for x86-64 Linux
with Rust 1.98.0; no Intel runtime measurement was made. The full release
library suite ran in the sandbox: 4,044 passed, nine local-socket tests failed
with permission denied, and 114 were ignored. The collaboration group then
passed 37/37 with loopback access. A complete loopback-enabled rerun was
stopped under host load above 150, so a single full-suite pass is unverified.

== Measurement decision

The final source-bound Apple M4 Pro panel completed 24 AB/BA pairs for each
of GF(2^113), GF(2^131), GF(2^163), and GF(2^233), 500,000 multiplies per
arm per round. All checksums matched. GF(2^131) baseline rounds ranged from
27,826,500 to 932,235,083 ns; candidate rounds ranged from 16,801,958 to
694,301,375 ns. The unchanged four-limb control also varied widely, and
host load exceeded 100. No wall-time advantage or regression passes the
controlled comparison gate.

#image("paired_mac_final.svg", width: 100%)

No Intel or GPU job completed. The sandboxed Modal attempt could not connect;
automatic approval review rejected its network-enabled source-upload retry
without specific export authorization. NIST primes, Curve25519, Ed25519,
Koblitz scalar recoding, Metal and CUDA/HIP remain open. The ECC2K boundary
and ecbench claim charts were checked but receive no new performance row.
