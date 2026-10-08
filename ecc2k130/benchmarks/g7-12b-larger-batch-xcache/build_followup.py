#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor
import json

from build import build, CONFIGS, OUT

CONFIGS.update({
    "batch24-cache6-min3": (24, 6, 3),
    "batch32-cache8-min3": (32, 8, 3),
})

labels = ("batch24-cache6-min3", "batch32-cache8-min3")
with ThreadPoolExecutor(max_workers=2) as pool:
    rows = list(pool.map(build, labels))
assert all(row["returncode"] == 0 for row in rows)
(OUT / "followup-builds.json").write_text(json.dumps(rows, indent=2) + "\n")
