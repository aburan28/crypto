#!/usr/bin/env python3
"""Turn `ic --out report.json bench --sweep ...` reports into the lead-ladder table.

  python3 -I summarise_runs.py runs/<stamp>/report_n*.json [--baseline-label SUBSTR]

For every report: one line per configuration with the medians over its repeats of the
structural metrics (LEAD_LADDER_METHODOLOGY_20261010.md section 1) and the per-phase
group-addition equivalents, then the ratio of every configuration to the baseline (the
first configuration of the sweep unless --baseline-label names another).  The second
copy of the baseline in an A/B/A sweep is the noise floor and is reported as such.
Operation counts only; wall seconds are printed in a separate column marked as a note.
"""
import argparse, json, statistics, sys

PHASES = [("fb", lambda r: r["factor_base"]["cost"]["gae"]),
          ("setup", lambda r: r["decomposition"]["setup"]["gae"]),
          ("decomp", lambda r: r["decomposition"]["cost"]["gae"]),
          ("la", lambda r: r["linear_algebra"]["cost"]["gae"]),
          ("verify", lambda r: r["verify"]["gae"]),
          ("total", lambda r: r["total_gae"])]

STRUCT = [("signed_pts", lambda r: r["factor_base"]["signed_points"]),
          ("columns", lambda r: r["factor_base"]["columns"]),
          ("pts/col", lambda r: r["factor_base"]["points_per_column"]),
          ("tried", lambda r: r["decomposition"]["targets_tried"]),
          ("found", lambda r: r["decomposition"]["relations_found"]),
          ("trials/rel", lambda r: r["decomposition"]["targets_tried"] / max(1, r["decomposition"]["relations_found"])),
          ("rows", lambda r: r["linear_algebra"]["rows"]),
          ("rank", lambda r: r["linear_algebra"]["rank"]),
          ("la_work", lambda r: r["linear_algebra"]["work"]),
          ("canon", lambda r: _native(r, "canonicalisations")),
          ("frob_maps", lambda r: _native(r, "frobenius_maps")),
          ("lookups", lambda r: _native(r, "lookups")),
          ("tbl_entries", lambda r: _native(r, "pair_table_entries")),
          ("fb_as_solves", lambda r: r["factor_base"]["cost"]["native"].get("as_solves", 0)),
          ("setup_adds", lambda r: r["decomposition"]["setup"]["group_ops"]["adds"]),
          ("decomp_adds", lambda r: r["decomposition"]["cost"]["group_ops"]["adds"]),
          ("decomp_doubles", lambda r: r["decomposition"]["cost"]["group_ops"]["doubles"]),
          ("decomp_smults", lambda r: r["decomposition"]["cost"]["group_ops"]["scalar_mults"]),
          ("la_row_ops", lambda r: r["linear_algebra"]["cost"]["native"].get("row_ops", 0)),
          ("S", lambda r: r["s"]),
          ("S/rho", lambda r: r.get("s_over_rho")),
          ("wall_s(note)", lambda r: r["wall_seconds"])]

def _native(r, key):
    tot = 0
    for part in (r["decomposition"]["setup"], r["decomposition"]["cost"], r["factor_base"]["cost"]):
        tot += part.get("native", {}).get(key, 0)
    return tot

def med(xs):
    xs = [x for x in xs if isinstance(x, (int, float))]
    return statistics.median(xs) if xs else None

def fmt(x):
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:.3g}" if abs(x) < 1e6 else f"{x:.3e}"
    return str(x)

def summarise(path, baseline_label):
    doc = json.load(open(path))
    rows = doc["rows"]
    by_label = {}
    order = []
    for r in rows:
        if r["label"] not in by_label:
            order.append(r["label"])
        by_label.setdefault(r["label"], []).append(r)
    rho = doc.get("rho_reference") or {}
    print(f"\n## {path}\n")
    print(f"instance `{doc['instance']}`  log2 r = {rows[0]['log2_r']:.1f}  status {doc['status']}  "
          f"rho reference: {rho.get('method','—')} A={rho.get('automorphisms','—')} mean S={fmt(rho.get('mean_s'))} "
          f"(plain A=1: S={fmt((doc.get('rho_reference_plain') or {}).get('mean_s'))})")
    pins = doc.get("calibration_pins") or {}
    print(f"calibration: pinned units {pins.get('pinned')}, measured on this host {pins.get('measured')} "
          f"(GAE of a measured unit is host-dependent; the counted columns are not)")
    skipped = doc.get("configurations_skipped") or []
    for s in skipped:
        print(f"skipped: {s['configuration']}: {s['why']}")
    # Which label is the baseline
    base = order[0]
    if baseline_label:
        for l in order:
            if baseline_label in l:
                base = l
                break
    # Medians per label
    table = {}
    for l in order:
        rs = by_label[l]
        table[l] = {
            "repeats": len(rs), "verified": sum(1 for r in rs if r["verified"]),
            "exhausted": sum(1 for r in rs if r["exhausted"]),
            **{k: med([f(r) for r in rs]) for k, f in STRUCT},
            **{"gae_" + k: med([f(r) for r in rs]) for k, f in PHASES},
        }
    cols = ["repeats", "verified", "exhausted"] + [k for k, _ in STRUCT] + ["gae_" + k for k, _ in PHASES]
    print("\n### medians over repeats\n")
    print("| configuration | " + " | ".join(cols) + " |")
    print("|---|" + "|".join("--:" for _ in cols) + "|")
    for l in order:
        print(f"| `{l}` | " + " | ".join(fmt(table[l][c]) for c in cols) + " |")
    # Ratios to baseline
    print(f"\n### ratio to baseline `{base}` (lead / baseline; the repeated baseline row is the noise floor)\n")
    rcols = ["columns", "pts/col", "trials/rel", "rows", "la_row_ops", "canon", "setup_adds", "decomp_adds"] + ["gae_" + k for k, _ in PHASES]
    print("| configuration | " + " | ".join(rcols) + " |")
    print("|---|" + "|".join("--:" for _ in rcols) + "|")
    b = table[base]
    for l in order:
        out = []
        for c in rcols:
            x, y = table[l][c], b[c]
            out.append(fmt(x / y) if (x is not None and y) else "—")
        print(f"| `{l}` | " + " | ".join(out) + " |")
    # Non-genericity quantity B*D^2/N with D = decomposition+LA+verify gae per target (online, excluding setup)
    print("\n### B·D²/N (B = signed points, D = online GAE per target = decomp + la + verify, N = r)\n")
    r_val = rows[0]["r"]
    for l in order:
        t = table[l]
        D = (t["gae_decomp"] or 0) + (t["gae_la"] or 0) + (t["gae_verify"] or 0)
        B = t["signed_pts"] or 0
        print(f"- `{l}`: B={fmt(B)} D={fmt(D)} → B·D²/N = {fmt(B * D * D / r_val)}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--baseline-label", default=None)
    a = ap.parse_args()
    for p in a.reports:
        summarise(p, a.baseline_label)

if __name__ == "__main__":
    main()
