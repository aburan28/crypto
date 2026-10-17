"""Presentation metadata only: quote the native bounded comparison verbatim."""
import json
import re
from pathlib import Path

repo = Path(__file__).resolve().parents[3]
comparison_path = "research/full_volcano_graphwide_ic_20261017/execution/pilot-summary/two-comparison.json"
cmp = json.loads((repo / comparison_path).read_text())
ratio = cmp["ops"]["ratio_b_over_a"]
assert cmp["ops"]["bounded"] and cmp["wall"]["status"] == "descriptive"
path = repo / "docs/ic/progress-timeline.json"
doc = json.loads(path.read_text())
series = next(s for s in doc["series"] if s["id"] == "ec-cold")
point = {"date": "2026-10-17", "value": ratio, "reference": "matched",
         "unit": "quotient of incomplete cold GAE lower bounds; diagnostic",
         "label": "Volcano projected-base control/pilot; best bounded IC cell",
         "detail": "Projected size-16 MITM on icv1-fp57-t536870914-6cabb0c9, two public one-target attempts, one round, L0/debug. Both IC and rho have unpriced work. The ratio is a quotient of incomplete lower bounds, not a direction-of-cost or speedup claim. Full graph-wide completion and exact representations remain open.",
         "class": "accounting", "anchor": "#lab-volcano-ic",
         "source": "https://github.com/aburan28/crypto/blob/main/" + comparison_path}
assert not any(p["label"] == point["label"] for p in series["points"])
series["points"].append(point)
doc["updated"] = "2026-10-17"
path.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
page = repo / "docs/index-calculus-scoreboard.html"
text = page.read_text()
panel = '''<h3 id="lab-volcano-ic">Volcano subgroup IC: audited controls and bounded survey, 17 October 2026</h3>
<p>The projected-abscissa base charges cofactor projection and checks membership and recovered base logs. The root-vertex pilot has 48 attempts (30 verified, 14 exhausted, 4 timed out); the separate all-vertex survey's first 12 vertices/strata have 192 attempts (119 verified, 58 exhausted, 15 timeouts). All successful and exhausted deterministic reports replay. The queue covers 24,823 vertices and 51 strata, but queue coverage is not completed execution. Exact registration, qualified full-matrix timing, finite-field IC and the original A–E scope remain open. L0/debug counts only; every arm has unpriced native work. <a href="https://github.com/aburan28/crypto/blob/main/research/full_volcano_graphwide_ic_20261017/execution/README.md">Protocol, raw records, per-stratum min/max/mean and inspected reports</a>. No fastest-method or scaling claim; the best-of table and qualified frontier are unchanged.</p>
'''
needle = '<h3>Binary GHS curve panel: exact model and subgroup differences</h3>'
assert text.count(needle) == 1
text = text.replace(needle, panel + needle)
row = f'<tr><td>2026-10-17</td><td>Elliptic curve, whole method, cold, one target</td><td>Volcano projected size-16 MITM, two targets; <code>icv1-fp57-t536870914-6cabb0c9</code></td><td class="n">{ratio:.3f}&times; <small>quotient of incomplete GAE lower bounds; diagnostic, no cost-direction claim</small></td><td>matched; L0/debug</td><td><span class="chip">accounting</span></td><td><a href="#lab-volcano-ic">panel</a> &middot; <a href="{point["source"]}">native comparison</a></td></tr>\n'
marker = '<script type="application/json" id="progress-data">'
assert text.count(marker) == 1
before, rest = text.split(marker, 1)
# The last table before the timeline script is the plotted-point source table.
end = before.rfind('</tbody>')
assert end >= 0
before = before[:end] + row + before[end:]
old_json, after = rest.split('</script>', 1)
json.loads(old_json)
page.write_text(before + marker + json.dumps(doc, indent=2, ensure_ascii=False) + '</script>' + after)
