# Lays out the timing stage's medians (ANALYSIS.json, written by
# `pdp_rootfind analyze`) as Typst tables for report.typ.  Layout only.
import json, math, os
here = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(here, "..", "..", "ANALYSIS.json")))
def num(x, nd=1):
    return "—" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:,.{nd}f}"
out = ["#table(columns: 9, align: (left, right, right, right, right, right, right, right, right),",
       "  [arm], [$n$], [$ell$], [ns / target, median], [ns / target, min], [rounds pooled], [contended (excluded)], [mul / target], [lookups / target],"]
for r in sorted(A["binary"], key=lambda r: (r["n"] != 3 * r["l"], r["l"], r["arm"])):
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}],".format(r["arm"].replace("_", "\\_"), r["n"], r["l"], num(r.get("ns_median"), 1), num(r.get("ns_min"), 1), r.get("ns_runs", 0), r.get("ns_contended_runs", 0), num(r.get("mul"), 0) if r.get("mul") else "—", num(r.get("lookups"), 0) if r.get("lookups") else "—"))
out.append(")")
open(os.path.join(here, "table_binary.typ"), "w").write("\n".join(out) + "\n")
out = ["#table(columns: 8, align: (right, left, right, right, right, right, right, right),",
       "  [$m$], [arm], [$t$], [ns / full scan, median], [ns / full scan, min], [rounds pooled], [mul / target], [table build s],"]
for r in sorted(A["prime"], key=lambda r: (r["m"], r["arm"], r["t"])):
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}],".format(r["m"], r["arm"], r["t"], num(r.get("ns_median"), 1), num(r.get("ns_min"), 1), r.get("ns_runs", 0), num(r.get("mul"), 0), num(r.get("setup_secs"), 4) if r.get("setup_secs") else "—"))
out.append(")")
open(os.path.join(here, "table_prime.typ"), "w").write("\n".join(out) + "\n")
print("ok")
