"""Tabulates the recorded probe lines (reports/<run>/raw/*.jsonl): checks
per-target decisions agree across arms and stage units repeat across rounds,
and takes medians / ratios of recorded CPU times against the `scatter` arm
(the rewrite a host without BMI2 ran before) for display."""
import sys, json, glob, statistics, collections
run = sys.argv[1]
ARMS = ("bmi2", "runs", "lazy", "scatter")
DECISION = ("hit", "indices", "exhausted", "reductions", "splits", "propagations",
            "infeasible", "geo_refute", "geo_witness", "partial_refute", "geo_adds")
dec = collections.defaultdict(set)
units = collections.defaultdict(set)
cpu = collections.defaultdict(list)
for f in sorted(glob.glob(f"{run}/raw/*.jsonl")):
    lines = [json.loads(l) for l in open(f)]
    cell, arm = lines[0]["cell"], lines[0]["arm"]
    for d in lines[1:]:
        if "target" in d:
            dec[(cell, d["arm"], d["target"])].add(tuple(json.dumps(d[k]) for k in DECISION))
        elif "total_ms" in d:
            units[(cell, d["arm"], arm)].add(d["word_xors"])
            cpu[(cell, d["arm"], arm)].append(d["cpu_ms"])
print(f"targets checked {len(dec)}; decisions differing across arms/rounds: {sum(len(v) != 1 for v in dec.values())}")
print(f"arm totals whose stage unit varies across rounds: {sum(len(v) != 1 for v in units.values())}")
order = ["K1_17", "K1_23", "K0_23", "K0_31", "K0_19", "K1_19", "K1_17d7", "K1_29", "K0_37"]
with open(f"{run}/table.tsv", "w") as f:
    f.write("cell\tengine\t" + "\t".join(f"units_{a}" for a in ARMS) + "\t" + "\t".join(f"cpu_{a}" for a in ARMS)
            + "\tC_scatter/runs\tC_scatter/runs_lo\tC_scatter/runs_hi\tC_lazy/runs\tC_bmi2/runs\tC_scatter/bmi2\n")
    for c in order:
        for e in ("f4", "f6p"):
            u = {a: next(iter(units[(c, e, a)])) for a in ARMS}
            m = {a: statistics.median(cpu[(c, e, a)]) for a in ARMS}
            sr = m["scatter"] / m["runs"]
            lo = min(cpu[(c, e, "scatter")]) / max(cpu[(c, e, "runs")])
            hi = max(cpu[(c, e, "scatter")]) / min(cpu[(c, e, "runs")])
            row = [c, e] + [str(u[a]) for a in ARMS] + [f"{m[a]:.1f}" for a in ARMS] + [
                f"{sr:.3f}", f"{lo:.3f}", f"{hi:.3f}", f"{m['lazy'] / m['runs']:.3f}",
                f"{m['bmi2'] / m['runs']:.3f}", f"{m['scatter'] / m['bmi2']:.3f}"]
            f.write("\t".join(row) + "\n")
            print(f"{c:8s} {e:4s} units runs==bmi2 {u['runs']==u['bmi2']} lazy==scatter {u['lazy']==u['scatter']} | CPU scatter/runs {sr:.2f} [{lo:.2f}-{hi:.2f}] lazy/runs {m['lazy']/m['runs']:.2f} bmi2/runs {m['bmi2']/m['runs']:.2f}")
