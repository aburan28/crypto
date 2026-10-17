# Binary field arithmetic panel: Apple M4 Pro

**Run ID:** `RUN-ECC-ARITH-20261013-MAC-FINAL`

**Status:** correctness checked; wall-time comparison rejected for host contention

**Scope:** `F2mElement::mul` in GF(2^113), GF(2^131), GF(2^163), and GF(2^233)

## Work and provenance

The study used 256 deterministic operand pairs per field from seed
`0x0123456789abcdef`, 500,000 multiplications per arm per round, 24 paired
rounds per field, with AB/BA order alternating. Each timed input was checked
against the independent bitwise `schoolbook_mul` first. The full source file
is imported into the thin Rust driver; the algorithm is not reimplemented in
the measurement harness. The candidate is a six-product three-limb
carryless multiply and a three-product two-limb multiply. The four-limb path
is an unchanged-workload control.

The host was MacBook Pro Mac16,7, Apple M4 Pro (10 performance plus four
efficiency cores), macOS 26.6, Rust 1.93.1. CPU placement and competing
workloads were uncontrolled. The observed system load at 11:41 PDT was
101.02, 99.46, 105.70 for the one-, five-, and fifteen-minute averages.
The base repository commit was `11459b07f62f2e3293669d9cd0f3c52fd07bcc19`.
The working tree was dirty with this candidate source and new study files;
the baseline source was extracted unchanged from that base commit.
The baseline source SHA-256 was
`2bc2a1f55887ef7db1fdb1ce5a6cd05c86b7c75e9eb8ade36497ff5686b37fd3`;
the final candidate source SHA-256 was
`8b32c5467ee8c87752b9cb86594fb742aa9053fd5918172686c4270c9e936364`.

The exact command was:

```sh
CARGO_TARGET_DIR=/private/tmp/ecc_arith_thin_target cargo run --offline --release --manifest-path research/ecc_arithmetic_20261013/Cargo.toml --bin compare -- 500000 24
```

The first separate-arm runs are retained in
[`baseline.csv`](../../baseline.csv) and [`candidate.csv`](../../candidate.csv).
The first paired run is retained in [`paired_mac.csv`](../../paired_mac.csv).
It and the separate candidate run used a pre-format candidate source with
SHA-256 `b9a2088ad7159376e5e7a2020a8618d2b66000e61e001e98d8c97e643b6428d7`;
that exact file was not frozen, so those rows are exploratory only.
The final source-bound paired run is
[`paired_mac_final.csv`](../../paired_mac_final.csv) with
[`compiler and launch log`](../../paired_mac_final.stderr). Its SHA-256 is
`cdf372b0e1840bde27c053d50a9949f1025c84ce1f09ccc4cad552869fa3dd25`.
The source to redraw the graph is [`plot.py`](../../plot.py).

## Observations

All 192 final timed arm rows completed and all checksums matched within each
field. The candidate's source-level carryless product count changes from
four to three for two limbs and nine to six for three limbs. This is an
operation count, not a measured wall-time advantage.

The final GF(2^131) baseline arm ranged from 27,826,500 to 932,235,083 ns
per 500,000 products; its candidate arm ranged from 16,801,958 to
694,301,375 ns. The four-limb control, which takes the same schoolbook path
in both versions, also swung widely. This host does not meet a wall-time
comparison gate. The raw data are plotted without filtering below.

![All paired rounds, log elapsed time](../../paired_mac_final.svg)

The candidate's focused release tests passed: 17 passed, zero failed
([test output](../../candidate_tests_final.stdout)). Both frozen sources
compiled for x86-64 Linux with Rust 1.98.0
([cross-check log](../../x86_cross_check_attempt2.stderr)); this is not an
Intel runtime test. The full release library suite compiled and ran in the
sandbox: 4,044 passed, nine collaboration tests failed because local socket
creation returned `Operation not permitted`, and 114 were ignored
([output](../../full_lib_tests_attempt4.stdout)). The collaboration group,
including those nine tests, passed 37/37 with loopback access
([output](../../pollard_collab_tests.stdout)). A complete loopback-enabled
rerun was stopped under host load above 150
([interruption record](../../full_lib_tests_attempt5_interruption.txt));
there is no single completed full-suite pass. Earlier builds hit system-volume
ENOSPC and a machine-load stop; a sparse checkout attempt lacked compile-time
fixtures, all retained in this study. The first Modal attempt could not connect from the sandbox
([stderr](../../modal_cli_attempt1.stderr)). Automatic approval review
rejected the network-enabled retry because it would upload the source bundle
to Modal without specific source-export authorization. No remote worker was
started by that retry.

## Accounting and limits

The timing includes allocation, runtime feature dispatch, carryless product,
and reduction in the public multiply. It excludes point addition, scalar
multiplication, Frobenius recoding, host-to-GPU transfer, and any ECC2K walk.
No Intel, Metal, CUDA, ROCm, or other prime-curve timing was obtained in this
run. The requested worldwide-best comparison remains open. The unchanged
four-limb control and the large within-arm variation explain why no
quantitative speedup is claimed here.

The checked canonical charts were the ECC2K boundary scoreboard and ecbench
comparisons. Neither receives a new performance row because the wall-time
gate failed; this study graph is the current visual record.
