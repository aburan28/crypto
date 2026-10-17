"""Reconstruct the exact pre-fix field source for the ARM64 branch audit."""

import hashlib
from pathlib import Path

here = Path(__file__).resolve().parent
source = (here / "../../../src/ecc/f25519.rs").resolve().read_text()
barrier = """        // Keep the all-zero/all-one mask opaque to LLVM. Without this
        // barrier it turned the canonical selection into data-dependent
        // branches on ARM64 even though the Rust expression used bitmasks.
        let mask = core::hint::black_box(0u64.wrapping_sub(1 - borrow));"""
original = "        let mask = 0u64.wrapping_sub(1 - borrow);"
assert source.count(barrier) == 1
raw = source.replace(barrier, original).encode()
assert hashlib.sha256(raw).hexdigest() == "960c2a11e03b81d3b9da3563e611d8b1a55118cd555f73bbc312d401e5dd52c0"
(here / "f25519_before_barrier.rs").write_bytes(raw)
(here / "ecc_before.rs").write_text(
    '#[path = "f25519_before_barrier.rs"]\npub(super) mod f25519;\n'
    '#[path = "../../../src/ecc/x25519.rs"]\npub mod x25519;\n'
)
driver = (here / "ecc_x25519_check.rs").read_text()
assert driver.startswith("mod ecc;")
(here / "ecc_x25519_before_check.rs").write_text(
    driver.replace("mod ecc;", '#[path = "ecc_before.rs"]\nmod ecc;', 1)
)
