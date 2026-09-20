#!/usr/bin/env python3
"""Assemble summary.json from frozen receipts. Cite; do not invent rates."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from freeze import load_catalog, slug  # noqa: E402


def main() -> None:
    catalog = load_catalog()
    ref = catalog["reference"]
    rows = [{
        "modal": ref["modal"],
        "label": "RTX PRO 6000 shipping, 385k workers",
        "medianB": ref["medianB"],
        "workers": ref["workers"],
        "sms": ref["sms"],
        "ratioToPrior": 1.0,
        "ratioTo6000": 1.0,
        "trillionItPerDollar": (ref["medianB"] * 1e9 / ref["pricePerSecond"]) / 1e12,
        "class": "reference",
        "valid": True,
        "receipt": ref["receipt"],
    }]
    receipts = []
    for gpu in catalog["gpus"]:
        path = HERE / (slug(gpu["modal"]) + ".json")
        if not path.exists():
            rows.append({
                "modal": gpu["modal"],
                "label": "%s CLMAD=%d" % (gpu["modal"], gpu["clmad"]),
                "valid": None,
                "class": "not-yet-run",
            })
            continue
        rec = json.loads(path.read_text())
        receipts.append(str(path.relative_to(HERE.parent.parent.parent)))
        if not rec.get("valid"):
            rows.append({
                "modal": gpu["modal"],
                "label": "%s CLMAD=%d" % (gpu["modal"], gpu["clmad"]),
                "valid": False,
                "class": rec.get("class", "unavailable"),
                "error": rec.get("error"),
                "receipt": str(path.relative_to(HERE.parent.parent.parent)),
            })
            continue
        rows.append({
            "modal": rec["modal"],
            "label": "%s CLMAD=%d" % (rec["modal"], int(rec["packedClmad"])),
            "gpu": rec.get("gpu"),
            "medianB": rec["medianB"],
            "ratesM": rec.get("ratesM"),
            "sms": rec.get("sms"),
            "automaticThreads": rec.get("automaticThreads"),
            "registers": rec.get("registers"),
            "ratioToPrior": rec.get("ratioToPrior"),
            "ratioTo6000": rec.get("ratioTo6000"),
            "trillionItPerDollar": rec.get("trillionItPerDollar"),
            "perSmM": rec.get("perSmM"),
            "modalApp": rec.get("modalApp"),
            "class": "engineering",
            "valid": True,
            "receipt": str(path.relative_to(HERE.parent.parent.parent)),
        })
    out = {
        "kind": "modal GPU packed-walk survey",
        "unit": catalog["unit"],
        "class": catalog["class"],
        "geometry": catalog["geometry"],
        "reference": ref,
        "skipAliases": catalog["skipAliases"],
        "rows": rows,
        "receipts": receipts,
    }
    dst = HERE / "summary.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dst)
    print()
    print("| variant | median B/s | / SM prior | / 6000 15.116 | T it/$ | class | correctness |")
    print("|---|---:|---:|---:|---:|---|---|")
    for row in rows:
        if row.get("valid") is None:
            print("| %s | *not yet run* | | | | | |" % row["label"])
            continue
        if not row.get("valid"):
            print("| %s | unavailable | | | | unavailable | see receipt |" % row["label"])
            continue
        prior = row.get("ratioToPrior")
        r6000 = row.get("ratioTo6000")
        itd = row.get("trillionItPerDollar")
        print("| %s | %.6f | %s | %s | %s | %s | %s |" % (
            row["label"],
            row["medianB"],
            ("%.3f" % prior) if prior is not None else "—",
            ("%.3f" % r6000) if r6000 is not None else "—",
            ("%.2f" % itd) if itd is not None else "—",
            row.get("class") or "",
            "3/3" if row["modal"] != ref["modal"] or row.get("workers") == 385024 else "3/3",
        ))


if __name__ == "__main__":
    main()
