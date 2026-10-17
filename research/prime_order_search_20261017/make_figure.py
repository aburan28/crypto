"""Plan wall time per prime_search curve, main vs prime_order (log scale).
Reads timing-main.txt and timing-new.txt (`arm bits rep seconds`), writes
plan-time.svg and report.html next to this file. No dependencies."""
import math, statistics, pathlib
here = pathlib.Path(__file__).parent
def load(name):
    d = {}
    for line in (here / name).read_text().split("\n"):
        if line.strip():
            _, b, _, s = line.split()
            d.setdefault(int(b), []).append(float(s))
    return {b: statistics.median(v) for b, v in sorted(d.items())}
main, new = load("timing-main.txt"), load("timing-new.txt")
W, H, L, R, T, B = 640, 400, 70, 150, 30, 50
bits = sorted(set(main) | set(new))
x = lambda b: L + (b - bits[0]) / (bits[-1] - bits[0]) * (W - L - R)
lo, hi = math.log10(0.01), math.log10(1000)
y = lambda s: T + (hi - math.log10(s)) / (hi - lo) * (H - T - B)
out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="sans-serif" font-size="12">',
       f'<rect width="{W}" height="{H}" fill="#fff"/>']
for e in range(-2, 4):
    v = 10 ** e
    out.append(f'<line x1="{L}" x2="{W-R}" y1="{y(v):.1f}" y2="{y(v):.1f}" stroke="#ddd"/>')
    out.append(f'<text x="{L-8}" y="{y(v)+4:.1f}" text-anchor="end" fill="#555">{v:g} s</text>')
for b in bits:
    out.append(f'<text x="{x(b):.1f}" y="{H-B+18}" text-anchor="middle" fill="#555">{b}</text>')
out.append(f'<text x="{(L+W-R)/2}" y="{H-12}" text-anchor="middle" fill="#333">prime_search curve size (bits)</text>')
out.append(f'<text x="18" y="{(T+H-B)/2}" transform="rotate(-90 18 {(T+H-B)/2})" text-anchor="middle" fill="#333">ecbench plan wall time (median)</text>')
for data, colour, label in ((main, "#c0392b", "main: point_count, O(p)"), (new, "#2471a3", "prime_order, O(√p)")):
    pts = " ".join(f"{x(b):.1f},{y(s):.1f}" for b, s in data.items())
    out.append(f'<polyline points="{pts}" fill="none" stroke="{colour}" stroke-width="2"/>')
    for b, s in data.items():
        out.append(f'<circle cx="{x(b):.1f}" cy="{y(s):.1f}" r="3.5" fill="{colour}"/>')
    lb, ls = list(data.items())[-1]
    out.append(f'<text x="{x(lb)+8:.1f}" y="{y(ls)+4:.1f}" fill="{colour}">{label}</text>')
out.append("</svg>")
(here / "plan-time.svg").write_text("\n".join(out) + "\n")
rows = "".join(f"<tr><td>{b}</td><td>{main.get(b, float('nan')):.3f}</td><td>{new[b]:.3f}</td><td>{main[b]/new[b]:.0f}×</td></tr>" for b in bits if b in main)
(here / "report.html").write_text(f"""<!doctype html><meta charset="utf-8"><title>Prime-order curve search</title>
<style>body{{font-family:sans-serif;max-width:720px;margin:24px auto;color:#222}}td,th{{padding:3px 10px;text-align:right}}</style>
<h1>Prime-order curve search: O(p) → O(√p)</h1>
<p>Wall time of <code>ecbench plan</code> on a one-curve <code>prime_search</code> spec (seed 59297), median of 3 runs (one run for main at 26 bits). The plan is dominated by <code>find_prime_order_curve</code>.</p>
{(here / 'plan-time.svg').read_text()}
<table><tr><th>bits</th><th>main (s)</th><th>this change (s)</th><th>speedup</th></tr>{rows}</table>
<p>Correctness: the search draws the same random numbers and accepts the same curves. All 21 committed <code>prime_search</code> sessions on main re-plan to their recorded workload IDs (288 workloads, <code>replan-check.txt</code>); <code>prime_order</code> equals point_count + primality on 400 random curves, and the search equals the naive search at 8–16 bits for three seeds (unit test <code>prime_order_matches_the_naive_count</code>).</p>
<p>Host: Intel Xeon @ 2.10 GHz, 4 vCPU cloud container, shared, unpinned, L0 wall. Engineering change to setup time; no solve, cost unit or measured result changes.</p>
""")
