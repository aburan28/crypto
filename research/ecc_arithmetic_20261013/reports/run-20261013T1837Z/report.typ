#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)
= Binary field arithmetic panel: Apple M4 Pro

*Run ID:* RUN-ECC-ARITH-20261013-MAC-FINAL  \
*Status:* correctness checked; wall-time comparison rejected for host contention  \
*Scope:* F2mElement::mul in GF(2^113), GF(2^131), GF(2^163), GF(2^233).

== Work and provenance

There were 256 deterministic operand pairs per field from seed
0x0123456789abcdef, 500,000 multiplications per arm per round, and 24 paired
rounds per field. AB and BA order alternated. Every timed input was checked
against an independent bitwise schoolbook multiplication before timing. The
thin driver imports the production source file.

Host: MacBook Pro Mac16,7, Apple M4 Pro (10 performance and four efficiency
cores), macOS 26.6, Rust 1.93.1. CPU placement and competing jobs were
uncontrolled. The system load at 11:41 PDT was 101.02, 99.46, 105.70.
Base commit: 11459b07f62f2e3293669d9cd0f3c52fd07bcc19.
The working tree was dirty with the candidate and new study files; the
baseline was extracted unchanged from that commit.
Baseline source SHA-256: 2bc2a1f55887ef7db1fdb1ce5a6cd05c86b7c75e9eb8ade36497ff5686b37fd3.
Final candidate source SHA-256: 8b32c5467ee8c87752b9cb86594fb742aa9053fd5918172686c4270c9e936364.

Command: `cargo run --offline --release --manifest-path research/ecc_arithmetic_20261013/Cargo.toml --bin compare -- 500000 24`
with CARGO_TARGET_DIR set to /private/tmp/ecc_arith_thin_target.
Raw final CSV: research/ecc_arithmetic_20261013/paired_mac_final.csv,
SHA-256 cdf372b0e1840bde27c053d50a9949f1025c84ce1f09ccc4cad552869fa3dd25.
Earlier separate-arm and paired rows used a pre-format candidate file that
was not frozen exactly and are exploratory only.

== Observations and validation

All 192 timed arm rows completed. The arm checksums matched for each field.
The source uses three instead of four carryless products for two limbs and
six instead of nine for three limbs. That count is not a wall-time claim.

GF(2^131) baseline ranged from 27,826,500 to 932,235,083 ns per 500,000
products; candidate ranged from 16,801,958 to 694,301,375 ns. The unchanged
four-limb control also varied widely. No wall-time comparison passes the
measurement gate on this host. The graph displays every final paired round.

#image("../../paired_mac_final.svg", width: 90%)

The candidate's focused release tests passed 17/17. Both frozen sources
compiled for x86-64 Linux, without an Intel runtime check. In the sandbox,
the release library suite had 4,044 pass, nine local-socket permission failures,
and 114 ignored. The collaboration group then passed 37/37 with loopback
access. A complete loopback-enabled rerun was stopped under host load above
150; one full-suite pass is unverified. Automatic approval review rejected
the Modal source-upload retry, so no remote worker ran.

== Accounting and open scope

The timer includes public field multiply allocation, dispatch, carryless
product and reduction. It excludes point arithmetic, scalar multiplication,
Frobenius recoding, data transfer and the ECC2K walk. Intel, Metal, CUDA,
ROCm, prime-curve and worldwide comparisons were not measured. The existing
ECC2K boundary scoreboard and ecbench charts were checked; neither gains a
performance row because the wall-time gate failed.
