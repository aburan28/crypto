"""Render every Modal paired round; separate the two ephemeral workers."""

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROWS = list(csv.DictReader((HERE / "modal_paired.csv").open()))
COLORS = {"baseline": "#375a9e", "candidate": "#d27732"}
FIELDS = ("f2m113", "f2m131", "f2m163", "f2m233")
parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="670" viewBox="0 0 1000 670">',
    '<rect width="1000" height="670" fill="white"/>',
    '<style>text{font-family:Helvetica,Arial,sans-serif;fill:#202633}.title{font-size:20px;font-weight:bold}.label{font-size:13px}.tick{font-size:11px;fill:#555}.grid{stroke:#d9dde4;stroke-width:1}.axis{stroke:#687181;stroke-width:1}</style>',
    '<text class="title" x="30" y="32">ECC binary field multiply: two Modal x86-64 CPU runs</text>',
    '<text class="label" x="30" y="54">Elapsed ms for 500,000 products; every paired round shown. Runs used separate ephemeral workers; CPU model unavailable.</text>',
]

for index, field in enumerate(FIELDS):
    col, row = index % 2, index // 2
    left, top = 75 + col * 480, 105 + row * 260
    width, height = 395, 185
    selected = [r for r in ROWS if r["field"] == field]
    ns = [int(r["nanoseconds"]) for r in selected]
    low, high = min(ns) * 0.94, max(ns) * 1.06
    scale_y = lambda value: top + height * (high - value) / (high - low)
    scale_x = lambda run, rnd: left + width * (int(rnd) + (int(run) - 1) * 26) / 49
    parts.append(f'<text class="label" x="{left}" y="{top - 14}">{field}</text>')
    for value in (low, (low + high) / 2, high):
        y = scale_y(value)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left + width}" y2="{y:.2f}"/>')
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 4:.2f}" text-anchor="end">{value / 1e6:.1f}</text>')
    divider = left + width * 24.5 / 49
    parts.append(f'<line x1="{divider:.2f}" y1="{top}" x2="{divider:.2f}" y2="{top + height}" stroke="#a9adb7" stroke-dasharray="4 4"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top + height}" x2="{left + width}" y2="{top + height}"/>')
    for run in (1, 2):
        center = scale_x(run, 11.5)
        parts.append(f'<text class="tick" x="{center:.2f}" y="{top + height + 18}" text-anchor="middle">run {run}: rounds 0–23</text>')
        for arm in ("baseline", "candidate"):
            arm_rows = sorted((r for r in selected if int(r["run"]) == run and r["arm"] == arm), key=lambda r: int(r["round"]))
            assert len(arm_rows) == 24
            points = [(scale_x(run, r["round"]), scale_y(int(r["nanoseconds"]))) for r in arm_rows]
            coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
            parts.append(f'<polyline points="{coords}" fill="none" stroke="{COLORS[arm]}" stroke-width="1.4"/>')
            parts.extend(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.8" fill="{COLORS[arm]}"/>' for x, y in points)

parts.extend([
    '<line x1="335" y1="626" x2="365" y2="626" stroke="#375a9e" stroke-width="2"/><text class="label" x="372" y="630">Baseline</text>',
    '<line x1="470" y1="626" x2="500" y2="626" stroke="#d27732" stroke-width="2"/><text class="label" x="507" y="630">Candidate</text>',
    '<text class="tick" x="30" y="656">Source: modal_paired.csv | Timing includes public field multiply allocation, dispatch and reduction | descriptive L0 measurement</text>',
    '</svg>',
])
(HERE / "modal_paired.svg").write_text("\n".join(parts) + "\n")
