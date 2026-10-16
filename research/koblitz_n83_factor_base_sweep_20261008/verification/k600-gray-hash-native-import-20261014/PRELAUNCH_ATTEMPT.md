# First native prelaunch attempt

- Source commit: `d5ee15191c4f9adf24677a07857f528e7753dc91`.
- Command: `uv run --with blake3 python -B research/koblitz_n83_factor_base_sweep_20261008/verification/k600-gray-hash-native-import-20261014/native_import.py --archive-dir /private/tmp/n83-k600-paired-archives-20261014 --output-dir /private/tmp/n83-native-k600-import-20261014`.
- Result: process exit 1, `ValueError: native toolchain version changed` in the version gate. The output directory was not created; no Rust build or import ran.
- Observed versions: `cargo 1.98.1 (797e8a9bc 2026-08-05)` and `rustc 1.98.1 (48a229cea 2026-09-01)`. The source gate expected `1.98.0` because the toolchain directory is named `1.98-aarch64-apple-darwin`.
- Classification: prelaunch runner error, not an empirical factor-base result.
