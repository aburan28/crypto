"""Parse hyperelliptic_ic_vs_rho sweeps and compare baseline (A*) with candidate (B*).

Usage: python3 analyse.py DIR   (reads DIR/A*.txt, DIR/B*.txt, DIR/bench.jsonl)
Writes DIR/summary.json and prints the tables.
"""
import glob
import json
import os
import re
import statistics as st
import sys

ROW = re.compile(
    r"^\s*(\d+)\s+(baseline|ab-walk|fb-walk)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+"
    r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)(.*)$"
)
DETAIL = re.compile(
    r"c = (\d+), relation stage (\d+) ops \((\d+) precompute, ([\d.]+) ops/trial\) "
    r"\+ oracle (\d+) mul-mods \((\d+) equiv measured, (\d+) charged\) "
    r"\+ linear algebra (\d+) mul-mods \((\d+) equiv\), conv (\d+); smoothness ([\d.]+); "
    r"IC floor S = ([\d.]+); wall (\d+) ms IC vs (\d+) ms rho"
)


def parse(path):
    rows, genus, cur = [], None, None
    for line in open(path):
        g = re.match(r"--- genus (\d+) ---", line.strip())
        if g:
            genus = int(g.group(1))
            continue
        r = ROW.match(line)
        if r:
            cur = dict(
                genus=genus, p=int(r[1]), config=r[2], N=int(r[3]), m=int(r[4]),
                ic_ops=float(r[5]), rho_ops=float(r[6]), S_ic=float(r[7]),
                S_rho=float(r[8]), ratio=float(r[9]), S_walk=float(r[10]),
                ic_over_floor=float(r[11]), ic_over_h1=float(r[12]),
                verified="UNSOLVED" not in r[13],
            )
            rows.append(cur)
            continue
        d = DETAIL.search(line)
        if d and cur is not None:
            cur.update(
                c=int(d[1]), relation_ops=int(d[2]), precompute=int(d[3]),
                oracle_modmuls=int(d[5]), oracle_equiv=int(d[6]),
                la_modmuls=int(d[8]), la_equiv=int(d[9]), conv=int(d[10]),
                smoothness=float(d[11]), floor_S=float(d[12]),
                wall_ic_ms=int(d[13]), wall_rho_ms=int(d[14]),
            )
    return rows


def key(r):
    return (r["genus"], r["p"], r["config"])


def main(d):
    runs = {}
    for f in sorted(glob.glob(os.path.join(d, "[AB][0-9].txt"))):
        runs[os.path.basename(f)[:-4]] = {key(r): r for r in parse(f)}
    bench = {}
    if os.path.exists(os.path.join(d, "bench.jsonl")):
        for line in open(os.path.join(d, "bench.jsonl")):
            rec = json.loads(line)
            bench[rec["label"]] = rec["run"]
    A = sorted(k for k in runs if k[0] == "A")
    B = sorted(k for k in runs if k[0] == "B")
    keys = sorted(runs[A[0]])

    # Pinned output: relation stage, oracle mul-mods, N, m and verification
    # must be identical across every sweep, baseline and candidate alike.
    pinned = ("N", "m", "relation_ops", "precompute", "oracle_modmuls", "verified", "rho_ops")
    violations = []
    for k in keys:
        ref = runs[A[0]][k]
        for lab in A + B:
            r = runs[lab].get(k)
            if r is None:
                violations.append((k, lab, "missing"))
                continue
            for f in pinned:
                if r.get(f) != ref.get(f):
                    violations.append((k, lab, f, ref.get(f), r.get(f)))

    out = {"sweeps": {lab: {"wall_seconds": bench.get(lab, {}).get("wall_seconds"),
                            "contended": bench.get(lab, {}).get("contended"),
                            "other_cpu_seconds": bench.get(lab, {}).get("other_cpu_seconds")}
                      for lab in A + B},
           "pinned_fields": list(pinned), "pinned_violations": violations, "rows": []}
    print(f"sweeps A={A} B={B}; contended: "
          f"{[l for l in A + B if bench.get(l, {}).get('contended')]}")
    print(f"pinned-output violations: {len(violations)}")
    for v in violations[:20]:
        print("  ", v)

    hdr = (f"{'g':>2} {'p':>4} {'config':>8} {'N':>9} {'m':>4} | {'LA_A med':>8} {'A/A':>9} "
           f"{'LA_B med':>8} {'cut':>6} | {'S_ic A':>7} {'S_ic B':>7} {'ratio A':>7} {'ratio B':>7} "
           f"| {'mm A':>7} {'mm B':>7}")
    print(hdr)
    for k in keys:
        la_a = [runs[l][k]["la_equiv"] for l in A]
        la_b = [runs[l][k]["la_equiv"] for l in B]
        row = dict(
            genus=k[0], p=k[1], config=k[2], N=runs[A[0]][k]["N"], m=runs[A[0]][k]["m"],
            la_equiv_A=la_a, la_equiv_B=la_b,
            la_modmuls_A=[runs[l][k]["la_modmuls"] for l in A],
            la_modmuls_B=[runs[l][k]["la_modmuls"] for l in B],
            S_ic_A=[runs[l][k]["S_ic"] for l in A], S_ic_B=[runs[l][k]["S_ic"] for l in B],
            ratio_A=[runs[l][k]["ratio"] for l in A], ratio_B=[runs[l][k]["ratio"] for l in B],
            S_rho_A=[runs[l][k]["S_rho"] for l in A], S_rho_B=[runs[l][k]["S_rho"] for l in B],
            S_walk_A=[runs[l][k]["S_walk"] for l in A], S_walk_B=[runs[l][k]["S_walk"] for l in B],
            oracle_equiv_A=[runs[l][k]["oracle_equiv"] for l in A],
            oracle_equiv_B=[runs[l][k]["oracle_equiv"] for l in B],
            wall_ic_A=[runs[l][k]["wall_ic_ms"] for l in A],
            wall_ic_B=[runs[l][k]["wall_ic_ms"] for l in B],
            wall_rho_A=[runs[l][k]["wall_rho_ms"] for l in A],
            wall_rho_B=[runs[l][k]["wall_rho_ms"] for l in B],
            floor_S=runs[A[0]][k]["floor_S"],
            relation_ops=runs[A[0]][k]["relation_ops"],
            verified_all=all(runs[l][k]["verified"] for l in A + B),
        )
        ma, mb = st.median(la_a), st.median(la_b) if la_b else None
        # A/A spread: max/min of the baseline sweeps for this row.
        aa = (max(la_a) / max(min(la_a), 1)) if la_a else None
        row["la_A_median"], row["la_B_median"], row["la_AA_spread"] = ma, mb, aa
        row["la_cut"] = (ma / mb) if mb else None
        row["la_cut_min_over_max"] = (min(la_a) / max(la_b)) if la_b and max(la_b) else None
        out["rows"].append(row)
        print(f"{k[0]:>2} {k[1]:>4} {k[2]:>8} {row['N']:>9} {row['m']:>4} | {ma:>8.0f} "
              f"{aa:>8.2f}x {mb if mb is not None else float('nan'):>8.0f} "
              f"{row['la_cut'] or float('nan'):>5.1f}x | {st.median(row['S_ic_A']):>7.2f} "
              f"{st.median(row['S_ic_B']) if B else float('nan'):>7.2f} "
              f"{st.median(row['ratio_A']):>7.2f} "
              f"{st.median(row['ratio_B']) if B else float('nan'):>7.2f} | "
              f"{st.median(row['la_modmuls_A']):>7.0f} "
              f"{st.median(row['la_modmuls_B']) if B else float('nan'):>7.0f}")
    json.dump(out, open(os.path.join(d, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
