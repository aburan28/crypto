# Capacity attempts for the frozen 4+2 pair

The branch-scoped GitHub Actions run [39054543374](https://github.com/aburan28/crypto/actions/runs/39054543374)
was queued from 17:37 UTC until an ARM runner began at 17:56 UTC on
2026-10-17. Its source commit is `17f03c518377468eeb49331ecfd8c3dedc2d72d0`.
The hosted workflow, rather than any alternative attempt below, is the
primary source of the eventual matched-arm evidence.

While queued, a bounded AWS `m7g.xlarge` fallback was prepared in
`us-west-2` with a 7,600-second process cap, instance-initiated terminate
behavior, a 50 GiB delete-on-termination root volume, and an exact-commit
bootstrap. The EC2 `RunInstances --dry-run` returned `DryRunOperation`.
The real `RunInstances` call returned `Blocked`: “This account is currently
blocked and not recognized as a valid account.” A subsequent tagged-instance
query returned `[]`, and `HeadObject` for its reserved S3 artifact key
returned 404. No EC2 execution or solver result came from this attempt. The
short-lived presigned upload URL was removed from local scratch after the
rejection. Its planned non-secret launch digest and key are retained in the
private local preparation receipt, not in the repository.

A short Modal read-only capacity probe observed `x86_64`, cgroup v1,
`memory.limit_in_bytes = 4294967296`, and CPU quota/period
`1000000/1000000`. The probe could not read a cgroup swap limit or
`/proc/swaps`; it therefore did not establish the strict worker's zero-swap
condition and was not used for the frozen pair. The Modal probe links are
[first](https://modal.com/apps/a-buran28/main/ap-exRa4nMMnBzSdFf5ODfqXg),
[cgroup](https://modal.com/apps/a-buran28/main/ap-4QBj0HFyt9FDkmu5WWqTb0),
and [swap](https://modal.com/apps/a-buran28/main/ap-5CWS2jmk83BAhTYfoWduqx).

Local macOS release-library recompilation exceeded 60-, 180-, and
600-second windows without a test result. A retained temporary build cache
was moved intact to the external SSD, increasing free system space from
about 102 MiB to about 2.2 GiB. Docker daemon access still timed out after
20 seconds outside the sandbox. These are infrastructure and verification
outcomes; they do not describe the N83 solver's behavior.

During final audit validation, another `cargo test --release --lib` attempt
remained at `Compiling crypto` without reaching test execution. It was stopped
when system-volume free space fell to 1.5 GiB, before an `ENOSPC` failure.
The focused Python evidence tests, 63-test N83 suite, 46-test boundary suite,
frozen-plan tests, workflow lint and independent hosted artifact replay passed.
This release-library attempt supplies no current-head pass or test failure.
