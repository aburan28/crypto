# New native panel environment

- Protocol commit: `2d3a7f3723925311e74a29ff1cdc7c74a0a97f50`
- Exporter binary: `target/release/examples/koblitz_n83_factor_base_export`
- Exporter binary SHA-256: `881ab00aa57d34f576ab5e935560beee88ee5b9d5f5a824f1f1b009bdef331f9`
- Compiler: `rustc 1.93.1 (01f6ddf75 2026-02-11) (Homebrew)`
- Cargo: `cargo 1.93.1 (Homebrew)`
- Operating system / architecture: `Darwin arm64`
- Isolation: local L0; other Rust compilers and research processes were observed on the shared host during the build and run interval.
- Build command: `cargo build --release --example koblitz_n83_factor_base_export`; completed in 21m 44s. Build time is outside the factor-base construction and verified-import intervals.
- Run command for each `N` in 1, 2, 3: `target/release/examples/koblitz_n83_factor_base_export pilot /private/tmp/n83-prebuilt-factor-base-20261016-run-N 300`, wrapped by `/usr/bin/time -l`.
- `/usr/bin/time -l` recorded an outer real interval but exited 1 after each completed child because sandboxed `sysctl kern.clockrate` was unavailable. The three `pilot.exit` values and stderr are retained; full native manifests record 54/54 completed objects in each child. The RSS field is unavailable.
- The temporary object directories were cleared across the chat continuation before replay. The three first replay attempts exited 1 on missing input files; their stderr and exit receipts are retained. A fourth, untimed-for-selection panel was constructed at `raw-reconstruction/`, inside the workspace. Its 54 compressed, plain and point-set hashes exactly match each original manifest and the historical S3 panel.
- For each original manifest, a separate recovery replay used a symlink to the single reconstructed byte-identical object set and `/usr/bin/time -p`. Every native replay returned `PASS` on 54 objects. This checks each original manifest against the reconstructed bytes; it is not an independent host or codebase.
- The process receipts, panel manifests, and replay summaries live under `runs/run-N/`. The reconstructed full point objects remain in the ignored `raw-reconstruction/` directory. Durability is claimed through the previously published S3 objects only after their exact compressed hashes matched. The local paths alone are not durable references.

The exact checksum binding to the historical S3 panel and the complete run
statuses are reported in `RESULT.md`.
