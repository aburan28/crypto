# Candidate release build: infrastructure failure

The first `cargo build --offline --locked --release --bin ecbench` exited 101 on 2026-10-14. Its terminal error was:

```text
error: failed to write /Volumes/SSD990/llm/tmp/rustciMv5Jq/lib.rmeta: No space left on device (os error 28)

error: could not compile `crypto` (lib) due to 1 previous error
```

The immediately following `df -h` showed 99 MiB available on `/Volumes/SSD990`, 841 GiB on `/Volumes/SSD990-2`, and 903 MiB on `/private/tmp`. The focused release unit test had passed before this failure. This is a compiler temporary-storage failure and has no mathematical or performance interpretation. The retry sets `TMPDIR` to `target/tmp` on `/Volumes/SSD990-2`; the failed build is retained here rather than replaced by its retry.
