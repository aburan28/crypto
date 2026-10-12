"""Render the W8 PDFs (layout only). usage: python -I render.py <study dir> [png-dir]"""
import os
import sys

import typst

root = sys.argv[1]
docs = [
    ("res", "W8_RESULTS.typ", "W8_RESULTS.pdf"),
    ("r1", "reports/run-20261011T2017Z/report.typ", "reports/run-20261011T2017Z/report.pdf"),
    ("r2", "reports/run-20261011T2131Z/report.typ", "reports/run-20261011T2131Z/report.pdf"),
]
for tag, src, out in docs:
    typst.compile(os.path.join(root, src), output=os.path.join(root, out))
    if len(sys.argv) > 2:
        os.makedirs(sys.argv[2], exist_ok=True)
        typst.compile(os.path.join(root, src), output=os.path.join(sys.argv[2], tag + "-{p}.png"), format="png", ppi=80)
print("rendered", [d[2] for d in docs])
