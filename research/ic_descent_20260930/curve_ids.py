#!/usr/bin/env python3
"""EC1 curve identities for §22's nine curves (docs/curve-identities.md).

A §20 parameter file names its curve only as `(a, n)`; the library's
`KoblitzCurve::new(a, n)` fixes the defining polynomial, the subgroup and
its generator.  `examples/koblitz_curve_records.rs` read those off the
library into `curve_records.json`; this hashes them with the repository's
reference helper, `tools/curve_identity.py`, and writes `curve_ids.json`.
It certifies nothing: equal hashes mean equal metadata.

    cargo run --release --example koblitz_curve_records -- \
        1:19 1:23 1:45 0:37 1:43 1:47 0:41 0:53 0:61 > curve_records.json
    python3 curve_ids.py > curve_ids.json

The readable tag is `e0` for `a = 0`, the family of the ECC2K-130
challenge curve (`y² + xy = x³ + 1`), and `e1` for `a = 1`
(`y² + xy = x³ + x² + 1`), a different Koblitz model (AGENTS.md §8b).
No candidate identity is assigned: this round records no factor-base
set digests, and unknown stays null.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from curve_identity import curve_identity  # noqa: E402


def proper_subfields(n: int) -> list[str]:
    """Proper intermediate subfields of GF(2^n) over GF(2): GF(2^d), 1 < d < n, d | n."""
    return [f"GF(2^{d})" for d in range(2, n) if n % d == 0]


def main() -> None:
    raw = (HERE / "curve_records.json").read_bytes()
    records = json.loads(raw)
    curves = []
    for rec in records["curves"]:
        a, n = rec["a"], rec["n"]
        ident = curve_identity(rec["field"], rec["curve"], f"e{a}")
        curves.append({
            "label": f"K_{a}/GF(2^{n})", "a": a, "n": n,
            "model": "E_0: y^2 + xy = x^3 + 1 (the ECC2K-130 family)" if a == 0
            else "E_1: y^2 + xy = x^3 + x^2 + 1 (not the ECC2K-130 family)",
            "proper_intermediate_subfields_over_gf2": proper_subfields(n),
            **ident,
            "field": rec["field"], "curve": rec["curve"], "group_order": rec["group_order"],
        })
    print(json.dumps({
        "convention": "docs/curve-identities.md, tools/curve_identity.py",
        "source": {"records": "curve_records.json",
                   "records_file_sha256": hashlib.sha256(raw).hexdigest(),
                   "producer": "examples/koblitz_curve_records.rs"},
        "candidate_uid": None,
        "curves": curves,
    }, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
