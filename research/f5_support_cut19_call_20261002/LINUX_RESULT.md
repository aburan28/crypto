# Fixed cut-19 physical Linux result: original gate failed

GitHub Actions run [37147179024](https://github.com/aburan28/crypto/actions/runs/37147179024) rebuilt PR #1292 on physical Linux x86-64 at its `6ea405635` head (the Actions merge SHA was `00de94c667f5898843e93dca69e60b8f97c0463d`). The runner passed AVX2/BMI2 preflight, release GF(2) and F5 tests, and the native benchmark build. The paired runner preserved all attempts and selected the first isolation-qualified complete attempt for each seed.

The selective-echelon / fixed-cut-19 n24 degree-4 **complete-call** one-thread medians were 2.967, 2.965, 2.962, and 2.984 times for `frozen`, `holdout_a`, `holdout_b`, and `holdout_c`. The respective exact five-pair bootstrap 95% lower bounds were 2.938, 2.948, 2.907, and 2.950. The two-thread medians were 2.520, 2.517, 2.489, and 2.564. Exact rank, canonical row space, criterion and structural counts matched; the candidate returned original rows via the cut-19 certificate.

**The frozen overall gate failed.** One-thread smaller-case wall medians were below that seed's A/A minimum on `holdout_a` and `holdout_b`; all four two-thread seed blocks also had at least one such control failure. These cases did not take the support certificate and had equal non-timing output and route fields. Their timing failures remain failures under the original protocol, even though the primary case exceeded 2×. Do not promote this run as a passing full gate.

The complete, compressed Actions artifact, including every attempted run, per-call output, isolation receipt, CPU choice and analyzer report, is `linux_run_37147179024.tar.gz` (SHA-256 `dd6c9066cfddd10f5eea6ba4d1cd4bc70c6fd7c62d549629bd971198fbad1a0e`). Extract with `tar -xzf linux_run_37147179024.tar.gz -C DIR`.
