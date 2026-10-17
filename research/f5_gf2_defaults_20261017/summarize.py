"""Tabulates reports/<run>/raw/*.jsonl: checks every arm returns the same
rows / bases / decisions and word counts as `on`, and takes medians over
rounds of the recorded times and their ratios against `on` for display."""
import sys, json, glob, statistics, collections
run = sys.argv[1]
f5 = collections.defaultdict(list)    # (form, threads, case, arm) -> wall_ms per round
f5key = collections.defaultdict(set)  # (form, threads, case) -> exactness keys over arms+rounds
f4 = collections.defaultdict(list); f4key = collections.defaultdict(set)
inh = collections.defaultdict(list); inhkey = collections.defaultdict(set)
DEC = ("hit", "indices", "exhausted", "reductions", "splits", "propagations", "infeasible", "word_xors")
for path in sorted(glob.glob(f"{run}/raw/*.jsonl")):
    lines = [json.loads(l) for l in open(path)]
    m = lines[0]
    for d in lines[1:]:
        if m["tool"] == "f5":
            k = (m["form"], m["threads"], d["case"])
            f5[k + (m["arm"],)].append(d["wall_ms"])
            f5key[k].add((d["rows_fp"], d["row_space_fp"], d["output_terms"], d["reduce_word_ops"], d["rank"]))
        elif m["tool"] == "f4":
            k = (d["case"],)
            f4[k + (m["arm"],)].append(d["wall_ms"])
            f4key[k].add((d["basis_fp"], d["word_xors"]))
        elif m["tool"] == "inh":
            if "target" in d:
                inhkey[(m["form"], d["target"])].add(tuple(json.dumps(d[x]) for x in DEC))
            elif "total_ms" in d:
                inh[(m["form"], m["arm"])].append(d["cpu_ms"])
bad = [k for k, v in list(f5key.items()) + list(f4key.items()) + list(inhkey.items()) if len(v) != 1]
print(f"exactness keys checked: f5 {len(f5key)}, f4 {len(f4key)}, inherited {len(inhkey)}; differing across arms/rounds: {len(bad)}")
med = lambda v: statistics.median(v)
rows = []
for (form, th, case, arm), v in sorted(f5.items()):
    if arm == "on":
        continue
    on = f5[(form, th, case, "on")]
    rows.append(("f5", f"form{form}", th, case, arm, med(on), med(v), med(v) / med(on), min(v) / max(on), max(v) / min(on)))
for (case, arm), v in sorted(f4.items()):
    if arm != "on":
        on = f4[(case, "on")]
        rows.append(("f4", "-", 1, case, arm, med(on), med(v), med(v) / med(on), min(v) / max(on), max(v) / min(on)))
for (cell, arm), v in sorted(inh.items()):
    if arm != "on":
        on = inh[(cell, "on")]
        rows.append(("inh", "-", 1, cell, arm, med(on), med(v), med(v) / med(on), min(v) / max(on), max(v) / min(on)))
with open(f"{run}/table.tsv", "w") as f:
    f.write("tool\tform\tthreads\tcase\tarm\ton_ms\tarm_ms\tarm/on\tlo\thi\n")
    for r in rows:
        f.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.4f}" for x in r) + "\n")
for r in rows:
    if r[0] != "f5" or r[3].startswith("f5_n2") or r[4] == "off":
        print(f"{r[0]:3s} {r[1]:5s} t{r[2]} {r[3]:16s} {r[4]:8s} on {r[5]:9.2f} arm {r[6]:9.2f} arm/on {r[7]:.3f} [{r[8]:.3f}-{r[9]:.3f}]")
