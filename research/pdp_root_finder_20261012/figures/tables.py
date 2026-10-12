# Lays out ANALYSIS.json (written by `pdp_rootfind analyze`) as Typst tables
# for REPORT.typ.  Formats recorded and derived values only; computes nothing.
import json, math, os

here = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(here, "..", "ANALYSIS.json")))

def esc(s):
    return str(s).replace("#", "\\#").replace("$", "\\$").replace("_", "\\_")

def num(x, nd=1):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{x:,.{nd}f}"

def sci(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    return f"${m:.2f} dot 10^({e})$"

def mb(b):
    return "—" if not b else f"{b / 1e6:,.1f}"

out = []
# ── binary table ──
rows = sorted(A["binary"], key=lambda r: (r["n"] != 3 * r["l"], r["l"], r["arm"]))
out.append("#table(columns: 12, align: (left, right, right, right, right, right, right, right, right, right, right, right),")
out.append("  [arm], [$n$], [$ell$], [mul / target], [lookups], [candidates], [÷ $T_\"enum\"$], [÷ pair floor], [µs / target (median)], [÷ baseline wall], [rounds (contended)], [table MB],")
for r in rows:
    ns = r.get("ns_median")
    counted = bool(r.get("mul"))  # baseline and group-scalar carry no counters
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{} ({})], [{}],".format(
        esc(r["arm"]), r["n"], r["l"], num(r.get("mul"), 0) if counted else "—", num(r.get("lookups"), 0) if counted else "—", num(r.get("candidates"), 0) if counted else "—",
        sci(r.get("ratio_candidates_over_T_enum")) if (counted and r.get("ratio_candidates_over_T_enum") is not None) else "—",
        num(r.get("ratio_lookups_over_pair_floor"), 2) if r.get("ratio_lookups_over_pair_floor") else "—",
        num(ns / 1e3, 1) if ns else "—",
        num(r.get("ratio_ns_vs_baseline"), 1) if r.get("ratio_ns_vs_baseline") else "—",
        r.get("ns_runs", 0), r.get("ns_contended_runs", 0), mb(r.get("table_bytes"))))
out.append(")")
open(os.path.join(here, "table_binary.typ"), "w").write("\n".join(out) + "\n")

# ── prime ladder table ──
out = []
rows = sorted(A["prime"], key=lambda r: (r["m"], r["arm"], r["t"]))
out.append("#table(columns: 11, align: (right, left, right, right, right, right, right, right, right, right, right),")
out.append("  [$m$], [arm], [$t$], [$N$], [mul / target], [candidates], [÷ $T_\"enum\"$], [÷ pair floor], [µs / full scan (median)], [rounds (contended)], [table MB],")
for r in rows:
    ns = r.get("ns_median")
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{} ({})], [{}],".format(
        r["m"], esc(r["arm"]), r["t"], int(r["N"]), num(r.get("mul"), 0), num(r.get("candidates"), 0),
        sci(r.get("ratio_candidates_over_T_enum")) if r.get("ratio_candidates_over_T_enum") is not None else "—",
        num(r.get("ratio_lookups_over_pair_floor"), 2) if r.get("ratio_lookups_over_pair_floor") else "—",
        num(ns / 1e3, 1) if ns else "—", r.get("ns_runs", 0), r.get("ns_contended_runs", 0), mb(r.get("table_bytes"))))
out.append(")")
open(os.path.join(here, "table_prime.typ"), "w").write("\n".join(out) + "\n")

# ── slopes ──
out = []
out.append("#table(columns: 7, align: (left, left, left, right, right, right, right),")
out.append("  [track], [arm], [quantity], [bits per rung], [± 95 %], [$r^2$], [sizes],")
for s in A["slopes"]:
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{}],".format(esc(s["track"]), esc(s["arm"]), esc(s["quantity"]), num(s["slope_bits_per_rung"], 3), num(s.get("ci95_halfwidth"), 3), num(s["r2"], 4), s["rungs"]))
out.append(")")
open(os.path.join(here, "table_slopes.typ"), "w").write("\n".join(out) + "\n")

# ── F4 vs root-finders on the committed p₁ Kummer cells ──
out = []
out.append("#table(columns: 10, align: (right, right, right, right, right, right, right, left, right, right),")
out.append("  [$m$], [$t$], [$N$], [finished rows], [unfinished], [F4 multiply-adds], [F4 ms], [arm], [mul / target], [F4 multiply-adds ÷ mul],")
for r in sorted([r for r in A["f4_prime"] if r["kind"] == "kummer" and r["p"] == 2013265921], key=lambda r: (r["m"], r["t"], r["arm"])):
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}], [{}],".format(
        r["m"], r["t"], r["N"], r["rows"], r.get("unfinished_rows", 0), sci(r["f4_muladds"]), num(r["f4_ms"], 1), esc(r["arm"]), num(r["mul"], 0), sci(r["ratio_muladds_over_mul"])))
out.append(")")
open(os.path.join(here, "table_f4.typ"), "w").write("\n".join(out) + "\n")

# ── S·T² ──
out = []
out.append("#table(columns: 7, align: (right, right, right, right, right, right, right),")
out.append("  [$n$], [$ell$], [$S$ (keys)], [$T$ (group-op eq.)], [$log_2(S T^2)$], [$log_2(r\\/|\"Aut\"|)$], [gap (bits)],")
for r in sorted(A["st2_binary"], key=lambda r: (r["n"], r["l"])):
    out.append("  [{}], [{}], [{}], [{}], [{}], [{}], [{:+.1f}],".format(r["n"], r["l"], sci(r["S_keys"]), num(r["T_gae"], 1), num(r["log2_ST2"], 1), num(r["log2_r_over_aut"], 1), r["gap_bits"]))
out.append(")")
open(os.path.join(here, "table_st2.typ"), "w").write("\n".join(out) + "\n")
print("ok")
