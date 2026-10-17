"""Tabulates the recorded probe lines (reports/<run>/raw/*.jsonl): checks
per-target decisions agree across arms and stage units repeat across rounds,
and takes medians / ratios of recorded CPU times for display."""
import sys, json, glob, os, statistics, collections
run = sys.argv[1]
DECISION = ("hit", "indices", "exhausted", "reductions", "splits", "propagations",
            "infeasible", "geo_refute", "geo_witness", "partial_refute", "geo_adds")
dec = collections.defaultdict(set)     # (cell, engine, target) -> decision tuples over arms+rounds
units = collections.defaultdict(set)   # (cell, engine, arm) -> stage units over rounds
cpu = collections.defaultdict(list)    # (cell, engine, arm) -> cpu_ms per round
skips = {}
cells = []
for f in sorted(glob.glob(f"{run}/raw/*.jsonl")):
    lines = [json.loads(l) for l in open(f)]
    meta = lines[0]; cell, arm = meta["cell"], meta["arm"]
    if cell not in cells: cells.append(cell)
    for d in lines[1:]:
        if "target" in d:
            dec[(cell, d["arm"], d["target"])].add(tuple(json.dumps(d[k]) for k in DECISION))
        elif "total_ms" in d:
            units[(cell, d["arm"], arm)].add(d["word_xors"])
            cpu[(cell, d["arm"], arm)].append(d["cpu_ms"])
            skips[(cell, d["arm"], arm)] = (d.get("gate_first_nodes"), d.get("gate_first_specialisations"))
bad = [k for k, v in dec.items() if len(v) != 1]
nonrep = [k for k, v in units.items() if len(v) != 1]
print(f"targets checked {len(dec)}; decisions differing across arms/rounds: {len(bad)}")
print(f"arm totals whose stage unit varies across rounds: {len(nonrep)}")
order = ["K1_17", "K1_23", "K0_23", "K0_31", "K0_19", "K1_19", "K1_17d7", "K1_29", "K0_37"]
rows = []
for c in order:
    for e in ("f6p", "f4"):
        u = {a: next(iter(units[(c, e, a)])) for a in ("off", "gate", "rref", "gaterref")}
        m = {a: statistics.median(cpu[(c, e, a)]) for a in u}
        r = lambda a: (u["off"] / u[a], m["off"] / m[a],
                       min(cpu[(c, e, "off")]) / max(cpu[(c, e, a)]), max(cpu[(c, e, "off")]) / min(cpu[(c, e, a)]))
        rows.append((c, e, u, m, {a: r(a) for a in ("gate", "rref", "gaterref")}, skips[(c, e, "gate")]))
hdr = "cell\tengine\tunits_off\tunits_gate\tunits_rref\tunits_gaterref\tcpu_off\tcpu_gate\tcpu_rref\tcpu_gaterref\tU_off/gate\tU_off/rref\tU_off/gaterref\tC_off/gate\tC_off/gate_lo\tC_off/gate_hi\tC_off/rref\tC_off/gaterref\tgate_skip_nodes\tgate_skip_specs"
with open(f"{run}/table.tsv", "w") as f:
    f.write(hdr + "\n")
    for c, e, u, m, r, s in rows:
        f.write("\t".join([c, e] + [str(u[a]) for a in ("off", "gate", "rref", "gaterref")]
                + [f"{m[a]:.1f}" for a in ("off", "gate", "rref", "gaterref")]
                + [f"{r[a][0]:.4f}" for a in ("gate", "rref", "gaterref")]
                + [f"{r['gate'][1]:.3f}", f"{r['gate'][2]:.3f}", f"{r['gate'][3]:.3f}", f"{r['rref'][1]:.3f}", f"{r['gaterref'][1]:.3f}"]
                + [str(s[0]), str(s[1])]) + "\n")
for c, e, u, m, r, s in rows:
    print(f"{c:8s} {e:4s} U off/gate {r['gate'][0]:.3f} off/rref {r['rref'][0]:.3f} off/both {r['gaterref'][0]:.3f} | CPU off/gate {r['gate'][1]:.2f} [{r['gate'][2]:.2f}-{r['gate'][3]:.2f}] off/rref {r['rref'][1]:.2f} off/both {r['gaterref'][1]:.2f}")
