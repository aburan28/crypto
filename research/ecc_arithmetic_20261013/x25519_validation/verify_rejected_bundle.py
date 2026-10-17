"""Verify the original rejected manifest against preserved historical bytes."""

import hashlib
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parents[2]
remap = {
    "src/ecc/f25519.rs": here / "f25519_before_barrier.rs",
    "research/ecc_arithmetic_20261013/x25519_validation/ecc_f25519_check.rs":
        here / "ecc_f25519_check_rejected.rs",
}
for line in (here / "MODAL_BUNDLE_rejected.sha256").read_text().splitlines():
    expected, name = line.split("  ", 1)
    path = remap.get(name, root / name)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual == expected, (name, actual, expected)
    print(f"{name}: OK")
