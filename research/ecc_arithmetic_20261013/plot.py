"""Lay out every recorded A/B round as an editable SVG; no filtering."""

import csv
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / "paired_mac_final.csv").open()))
parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="650" viewBox="0 0 1000 650">',
    '<rect width="1000" height="650" fill="white"/>',
    '<style>text{font-family:Helvetica,Arial,sans-serif;fill:#202633} .title{font-size:20px;font-weight:bold} .label{font-size:13px} .tick{font-size:11px;fill:#555} .grid{stroke:#d9dde4;stroke-width:1} .axis{stroke:#687181;stroke-width:1}</style>',
    '<text class="title" x="30" y="32">ECC binary field multiply: raw paired rounds on M4 Pro</text>',
    '<text class="label" x="30" y="54">Elapsed time for 500,000 products; log scale. All 24 pairs shown per field.</text>',
]

colors = {"baseline": "#375a9e", "candidate": "#d27732"}
for index, field in enumerate(["f2m113", "f2m131", "f2m163", "f2m233"]):
    col, row = index % 2, index // 2
    left = 75 + col * 480
    top = 95 + row * 255
    width, height = 395, 195
    parts.append(f'<text class="label" x="{left}" y="{top - 12}">{field}</text>')
    for value, label in [(10_000_000, "10 ms"), (100_000_000, "100 ms"), (1_000_000_000, "1000 ms")]:
        y = top + height * (1 - (math.log10(value) - 7) / (math.log10(2_000_000_000) - 7))
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left + width}" y2="{y:.2f}"/>')
        parts.append(f'<text class="tick" x="{left - 8}" y="{y + 4:.2f}" text-anchor="end">{label}</text>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top + height}" x2="{left + width}" y2="{top + height}"/>')
    for tick in [0, 6, 12, 18, 23]:
        x = left + width * tick / 23
        parts.append(f'<text class="tick" x="{x:.2f}" y="{top + height + 17}" text-anchor="middle">{tick}</text>')
    for arm in ["baseline", "candidate"]:
        selected = sorted((r for r in rows if r["field"] == field and r["arm"] == arm), key=lambda r: int(r["round"]))
        points = []
        for r in selected:
            x = left + width * int(r["round"]) / 23
            y = top + height * (1 - (math.log10(int(r["nanoseconds"])) - 7) / (math.log10(2_000_000_000) - 7))
            points.append((x, y))
        coords = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        parts.append(f'<polyline points="{coords}" fill="none" stroke="{colors[arm]}" stroke-width="1.4"/>')
        parts.extend(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.1" fill="{colors[arm]}"/>' for x, y in points)

parts.extend([
    '<line x1="370" y1="615" x2="400" y2="615" stroke="#375a9e" stroke-width="2"/><text class="label" x="407" y="619">Baseline</text>',
    '<line x1="490" y1="615" x2="520" y2="615" stroke="#d27732" stroke-width="2"/><text class="label" x="527" y="619">Candidate</text>',
    '<text class="tick" x="30" y="645">RUN-ECC-ARITH-20261013-MAC-FINAL | Source: paired_mac_final.csv | Round number on x-axis</text>',
    '</svg>',
])
(HERE / "paired_mac_final.svg").write_text("\n".join(parts) + "\n")
