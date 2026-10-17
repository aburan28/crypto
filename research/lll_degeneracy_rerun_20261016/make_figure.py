"""recovery.svg: share of probes that recovered the key, per curve, before
the Gram-Schmidt fixes (as the 2026-10-11 note reported them) and in the
2026-10-16 re-run (results.json).  Plots recorded outcomes only."""
import html
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
summary = json.loads((HERE / "results.json").read_text())["summary"]
# Pre-fix outcomes as printed in papers/lll-degeneracy-note (2026-10-11 PDF),
# Tables 1 and 3: (recovered, probes).
PRE = {"secp192k1": (1, 1), "secp224k1": (1, 1), "secp256k1": (0, 3), "P-256": (3, 3),
       "brainpoolP256r1": (1, 1), "P-384": (0, 1), "brainpoolP384r1": (0, 1), "P-521": (0, 1)}
CURVES = list(PRE)
post = {c: summary.get(f"{c} (f64)", {"recovered": 0, "probes": 0}) for c in CURVES}
hp = summary["P-521 (HP)"]
W, H, L, T, BW = 860, 380, 70, 60, 26
x0, plot_w, plot_h = L, W - L - 30, H - T - 90
out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Helvetica, Arial, sans-serif" font-size="12">',
       f'<rect width="{W}" height="{H}" fill="#fff"/>',
       '<text x="16" y="24" font-size="15" font-weight="bold">Boneh–Venkatesan HNP lattice: share of probes recovering the key</text>',
       '<text x="16" y="42" fill="#444">pre-fix: the 2026-10-11 note (unscaled f64 Gram–Schmidt); re-run: 2026-10-16, scaled f64 (P-521 also 2048-bit HP); labels are recovered/probes</text>']
for frac in (0, 0.5, 1):
    y = T + plot_h * (1 - frac)
    out.append(f'<line x1="{x0}" x2="{x0 + plot_w}" y1="{y}" y2="{y}" stroke="#ddd"/><text x="{x0 - 8}" y="{y + 4}" text-anchor="end">{frac:.1f}</text>')
step = plot_w / len(CURVES)
for i, c in enumerate(CURVES):
    cx = x0 + step * i + step / 2
    bars = [("#b3261e", PRE[c]), ("#2e7d32", (post[c]["recovered"], post[c]["probes"]))]
    if c == "P-521":
        bars.append(("#1565c0", (hp["recovered"], hp["probes"])))
    start = cx - BW * len(bars) / 2
    for j, (color, (r, n)) in enumerate(bars):
        frac = r / n if n else 0
        bx, h = start + j * BW, plot_h * frac
        out.append(f'<rect x="{bx:.1f}" y="{T + plot_h - h:.1f}" width="{BW - 3}" height="{max(h, 1.5):.1f}" fill="{color}"/>')
        out.append(f'<text x="{bx + (BW - 3) / 2:.1f}" y="{T + plot_h - h - 4:.1f}" text-anchor="middle" font-size="10">{r}/{n}</text>')
    out.append(f'<text x="{cx:.1f}" y="{T + plot_h + 18}" text-anchor="middle">{html.escape(c)}</text>')
ly = H - 30
for k, (color, label) in enumerate([("#b3261e", "pre-fix (as reported)"), ("#2e7d32", "re-run, scaled f64 GS"), ("#1565c0", "re-run, 2048-bit HP GS")]):
    out.append(f'<rect x="{x0 + k * 230}" y="{ly}" width="14" height="14" fill="{color}"/><text x="{x0 + k * 230 + 20}" y="{ly + 12}">{label}</text>')
out.append("</svg>")
(HERE / "recovery.svg").write_text("\n".join(out) + "\n")
print("recovery.svg")
