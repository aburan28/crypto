# Restore the routine CodeQL PR budget

At source head 2f3fad538b4c5955cd8ee0867d0022c70c94f649, the Rust CodeQL job completed successfully, but the PR budget check failed. The attempted 60-minute Rust exception violated the repository's proven 25-minute limit for routine pull_request and merge_group jobs. The original failure and full logs are preserved under raw/. The checker and its exceptions are unchanged.

The successful Rust job separated into about 10minutes of extraction/finalization and about 18minutes of query evaluation. The workflow now builds the complete Rust database in a 25-minute job and passes a relocatable, hashed database bundle to the existing analyze matrix job. Actual Rust queries keep their original Analyze(rust) check identity and /language:rust category. Other languages, push/scheduled Rust scans and the manual long lane retain their existing analysis behavior.

The receiving Rust job uses GitHub's official setup-codeql action, checks the bundle checksum, exact source commit and CodeQL version, then runs database analyze without a query override. The CLI uses the configuration retained at database creation; the bundle includes generated suites and diagnostics. Result upload remains mandatory. A missing/failed database job, artifact, checksum, version, source binding, query execution or upload fails the applicable Rust job. Both checkouts retain and prove every tracked Rust source and Cargo manifest. No source subset or query exclusion was added.

Validation is recorded in the adjacent raw outputs and validation-summary.json:

- The original checker reports 4violations before the change and 0after it in the authored head's 403-workflow private snapshot.
- The exact CI synthetic merge 921e2a69fbc519963e69031501a6dfe9a66205ec reproduces the same 4violations across 406workflows; applying only this workflow change gives 0violations. Its three additional workflows are listed in workflow-input-difference.json. Neither the source checkout nor base history was merged or rewritten.
- All 14 existing budget-policy controls pass. The original workflow parser, final actionlint run and authored-workflow whitespace check pass.
- A fresh cargo test --release --lib --locked --offline gate passes 4148 tests, with 0 failures and 117 ignored. Runtime was 422.43seconds, compilation 0.33seconds, descriptive operational values on a shared host. The actual test binary hash is bound in validation-summary.json.
- Root src/, Cargo.toml, Cargo.lock and build.rs have no runtime changes. The frozen native isogeny executable used for the separate M44497 experiment is unchanged.

The initial local snapshot contained 403 tracked workflow files; the CI log counted 406 because GitHub tested a synthetic merge into 1e03f1d0b9c77d72ee281561fd7f58405c4ae0b0. Both input sets were retained and checked. Initial metadata editing through the tool interface encountered a JavaScript parse error before any command or file write; the corrected edit and final validation are recorded here. Earlier valid intermediate actionlint/budget outputs remain alongside explicitly named final outputs. No test result was rewritten.

New-head CI and actual end-to-end execution of the split CodeQL jobs remain pending publication. Local metadata checks and the successful prior single-job analysis establish no new-head green status or measured speedup. The PR remains subject to every required applicable check before merge and official release.

Primary implementation references: [CodeQL database bundle](https://docs.github.com/en/code-security/reference/code-scanning/codeql/codeql-cli-manual/database-bundle), [database analyze and retained query configuration](https://docs.github.com/en/code-security/reference/code-scanning/codeql/codeql-cli-manual/database-analyze), and [official standalone setup action](https://github.com/github/codeql-action/blob/6f2ac81f5cd81fe351a197ba0b00d91d03cd2aee/setup-codeql/action.yml). The inspected upstream analyze implementation and its license/provenance are preserved under raw/.

