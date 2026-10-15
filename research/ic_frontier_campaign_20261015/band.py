"""Price each arm across the AGENTS.md §5a band and report where rankings hold.

Recorded values only: for a record, priced S(w) = (total_gae + w * U) / sqrt(r),
where U is the sum of the counts behind its unpriced counters (pilot/pricing.py)
and w is a weight in GAE per native operation.  A record with an uncounted
unpriced counter has no priced S above w = 0 and is reported as a lower bound.

    python3 band.py sessions/K13/records.jsonl [...]          # per-arm table
    python3 band.py --json out.json sessions/*/records.jsonl  # machine-readable

Ranking rule: arm A ranks below arm B on a curve only if mean S_A(w) < mean
S_B(w) at every w in WEIGHTS; otherwise the pair is "unresolved".
"""
import json
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "pilot"))
from pricing import unpriced_counts  # noqa: E402

WEIGHTS = (0.0, 0.1, 1.0)

# Native counters ecbench prices through a pinned calibration, with the unit
# that prices each (src/cryptanalysis/ecbench/methods.rs, PRICED_NATIVE), and
# the table those ratios come from (docs/ic/calibration.json, ns per unit with
# ns_per_add = 1, i.e. GAE per operation).
PRICED_NATIVE = {
    "sqrt_solves": "ns_per_sqrt", "as_solves": "ns_per_as_solve", "s4_pairs": "ns_per_s4_pair",
    "lookups": "ns_per_lookup", "row_ops": "ns_per_row_op", "legendre_symbols": "ns_per_legendre",
    "inversions": "ns_per_inversion", "frobenius_maps": "ns_per_frobenius",
    "canonicalisations": "ns_per_canon", "target_guard_probes": "ns_per_lookup",
}
_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
CALIBRATION = json.load(open(os.path.join(_ROOT, "docs", "ic", "calibration.json")))["instances"]


def calibrated_native(r):
    """(GAE the pinned calibration charged for native work, count of those operations)."""
    pinned = (r.get("detail") or {}).get("calibration_pinned_units") or []
    if not pinned:
        return 0.0, 0
    curve = r["workload"]["curve"]
    key = None
    if curve.get("family") == "koblitz":
        a = curve["construction"].split("(")[1].split(",")[0].strip()
        key = f"koblitz/K_{a} / GF(2^{curve['field_degree']})"
    table = CALIBRATION.get(key) if key else None
    if table is None:
        raise SystemExit(f"{r['record_id']}: pinned units {pinned} but no calibration row {key!r}")
    gae, ops = 0.0, 0
    for p in r["phases"]:
        for counter, unit in PRICED_NATIVE.items():
            n = p.get("native", {}).get(counter, 0)
            if n and unit in pinned:
                gae += n * table[unit]
                ops += n
    s = r.get("solver") or {}
    if s.get("ops") and "ns_per_word_xor" in pinned and "word XOR" in (s.get("op_unit") or ""):
        gae += s["ops"] * table["ns_per_word_xor"]
        ops += s["ops"]
    return gae, ops


def load(paths):
    rows = []
    for path in paths:
        for line in open(path):
            r = json.loads(line)
            if r.get("warmup"):
                continue
            c = r["cost"]
            curve = r["workload"]["curve"]
            sq = math.sqrt(int(curve["r"]))
            counts = unpriced_counts(r)
            missing = [k for k, v in counts.items() if v is None]
            u = sum(v for v in counts.values() if v is not None)
            # Take calibration-priced native work back out and price it with
            # everything else at the band weight, so every arm is priced alike.
            cal_gae, cal_ops = calibrated_native(r)
            base = c["total_gae"] - cal_gae
            u += cal_ops
            priced = {w: ((base + w * u) / sq if (w == 0.0 or not missing) else None)
                      for w in WEIGHTS}
            rows.append(dict(curve=curve["slug"], log2r=math.log2(int(curve["r"])), arm=r["arm"],
                             status=r["outcome"]["status"], target=r["workload"]["target_index"],
                             round=r["round"], s=c["s"], priced=priced, missing=missing,
                             unpriced_total=u, complete=not c["lower_bound"]))
    return rows


def summarise(rows):
    by = defaultdict(list)
    for x in rows:
        by[(x["curve"], x["arm"])].append(x)
    out = []
    for (curve, arm), xs in sorted(by.items(), key=lambda kv: (kv[1][0]["log2r"], kv[0][1])):
        ok = [x for x in xs if x["status"] == "verified"]
        means = {}
        for w in WEIGHTS:
            vals = [x["priced"][w] for x in ok]
            means[w] = (sum(vals) / len(vals)) if ok and all(v is not None for v in vals) else None
        out.append(dict(curve=curve, log2r=xs[0]["log2r"], arm=arm, n=len(xs), verified=len(ok),
                        mean=means, missing=sorted({m for x in xs for m in x["missing"]})))
    return out


def ranks(summary):
    """Pairs (a, b, verdict) per curve under the band rule."""
    by = defaultdict(list)
    for s in summary:
        by[s["curve"]].append(s)
    res = []
    for curve, arms in by.items():
        for i, a in enumerate(arms):
            for b in arms[i + 1:]:
                if a["verified"] < a["n"] or b["verified"] < b["n"]:
                    res.append((curve, a["arm"], b["arm"], "incomplete"))
                    continue
                ma, mb = a["mean"], b["mean"]
                if any(ma[w] is None or mb[w] is None for w in WEIGHTS):
                    res.append((curve, a["arm"], b["arm"], "unrankable (uncounted unpriced work)"))
                elif all(ma[w] < mb[w] for w in WEIGHTS):
                    res.append((curve, a["arm"], b["arm"], f"{a['arm']} < {b['arm']}"))
                elif all(ma[w] > mb[w] for w in WEIGHTS):
                    res.append((curve, a["arm"], b["arm"], f"{b['arm']} < {a['arm']}"))
                else:
                    res.append((curve, a["arm"], b["arm"], "unresolved inside the band"))
    return res


def main(argv):
    js = None
    if argv and argv[0] == "--json":
        js, argv = argv[1], argv[2:]
    rows = load(argv)
    summary = summarise(rows)
    if js:
        json.dump({"weights": WEIGHTS, "summary": summary, "ranks": ranks(summary)},
                  open(js, "w"), indent=1, default=str)
    print("| curve | log2 r | arm | verified/runs | S(w=0) | S(w=0.1) | S(w=1) | uncounted |")
    print("|:--|--:|:--|--:|--:|--:|--:|:--|")
    for s in summary:
        f = lambda v: "—" if v is None else f"{v:.3f}"
        print(f"| {s['curve']} | {s['log2r']:.2f} | {s['arm']} | {s['verified']}/{s['n']} | "
              f"{f(s['mean'][0.0])} | {f(s['mean'][0.1])} | {f(s['mean'][1.0])} | "
              f"{', '.join(s['missing']) or '—'} |")


if __name__ == "__main__":
    main(sys.argv[1:])
