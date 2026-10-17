"""Render the matched σ-table rate chart from the frozen headline receipt."""

import json
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "ecc2k130/benchmarks/sigma-table/headline/result.json"
data = json.loads(SOURCE.read_text())
assert data["correctness_preflight"] and data["aa_noise_gate"]

control = data["control_median_million_per_second"] / 1000
candidate = data["candidate_median_million_per_second"] / 1000
ratio = data["ab_median_ratio"]
ratios = data["ab_ratios"]
aa_drift = 100 * data["aa_max_absolute_drift"]

W, H = 1200, 600
X0, X1 = 305, 1075
LO, HI = 15.4, 16.0
fg, muted, blue, orange = "#172b40", "#52677a", "#2277a7", "#d77939"


def x(value):
    return X0 + (value - LO) * (X1 - X0) / (HI - LO)


def font(size, bold=False):
    base = "/System/Library/Fonts/Supplemental/Arial"
    return ImageFont.truetype(base + (" Bold.ttf" if bold else ".ttf"), size)


png = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(png)
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">', '<rect width="100%" height="100%" fill="white"/>']


def line(coords, color, width=2):
    draw.line(coords, fill=color, width=width)
    svg.append(f'<line x1="{coords[0]}" y1="{coords[1]}" x2="{coords[2]}" y2="{coords[3]}" stroke="{color}" stroke-width="{width}"/>')


def circle(cx, cy, r, color):
    draw.ellipse((cx-r, cy-r, cx+r, cy+r), fill=color)
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}"/>')


def text(px, py, label, size=24, color=fg, bold=False):
    draw.text((px, py), label, fill=color, font=font(size, bold))
    weight = "bold" if bold else "normal"
    svg.append(f'<text x="{px}" y="{py+size}" font-family="Arial, sans-serif" font-size="{size}" font-weight="{weight}" fill="{color}">{escape(label)}</text>')


text(65, 42, "Same-walk RTX PRO 6000 kernel rate", 37, fg, True)
text(65, 95, "Five paired A/B sessions; 385,024 workers; B16/T512/minBlocks1", 23, muted)

for tick in (15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0):
    xx = round(x(tick))
    line((xx, 180, xx, 430), "#e5ebef", 2)
    text(xx-24, 448, f"{tick:.1f}", 19, muted)

for name, value, yy, color in (("Fused control", control, 260, blue), ("Sigma-table", candidate, 355, orange)):
    text(65, yy-20, name, 24, fg, True)
    xx = round(x(value))
    line((X0, yy, xx, yy), color, 11)
    circle(xx, yy, 13, color)
    text(xx+23, yy-19, f"{value:.6f} B/s", 23, color, True)

text(390, 495, "Scalar updates per second (billions; zoomed axis)", 22, fg)
text(65, 545, f"Paired median ratio {ratio:.6f}; range {min(ratios):.6f}-{max(ratios):.6f}; A/A max drift {aa_drift:.3f}%", 19, muted)

svg.append("</svg>\n")
(HERE / "gpu-rate.svg").write_text("\n".join(svg))
png.save(HERE / "gpu-rate.png", dpi=(150, 150))
