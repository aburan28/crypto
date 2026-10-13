# Packs scripts/*.magma into calculator-sized submissions (bundles/NN.magma,
# at most 49,000 bytes each, the calculator's limit being 50,000), easiest
# first by our own solve time.  Light scripts (ours under 5 s) share bundles;
# a heavy one gets a bundle of its own with the other target of its cell,
# so a 60 s cut-off on one hard cell cannot cost another.  Every script
# prints its own RESULT line as it finishes.  Text layout only; nothing is
# computed here.
import json, os, glob

here = os.path.dirname(os.path.abspath(__file__))
ours = {}
for f in ("ours_k4.jsonl", "ours_pkm.jsonl", "ours_d4.jsonl"):
    for l in open(os.path.join(here, f)):
        r = json.loads(l)
        ours[r["name"]] = r.get("ours_total_ms", r.get("ours_ms"))
bool_ms = {}
for l in open(os.path.join(here, "ours_bool.jsonl")):
    r = json.loads(l)
    if r.get("arm") == "f4" and "ms" in r:
        bool_ms[f"bool_n17_root{r['target']}"] = r["ms"]
ours.update(bool_ms)
scripts = sorted(glob.glob(os.path.join(here, "scripts", "*.magma")),
                 key=lambda p: (ours.get(os.path.basename(p)[:-6], 1e12), p))
LIMIT = 49000
os.makedirs(os.path.join(here, "bundles"), exist_ok=True)
for old in glob.glob(os.path.join(here, "bundles", "*.magma")):
    os.remove(old)
def cell(name):
    # the two targets of a PKM cell share a bundle
    return name.rsplit("_x", 1)[0] if name.startswith("pkm_") else name

bundles, cur, size = [], [], 0
heavy = {}
for p in scripts:
    name = os.path.basename(p)[:-6]
    body = open(p).read()
    # Heavy cells, and the D₄ system (added after the first 14 bundles were
    # published, so it gets its own and leaves those byte for byte as they were).
    if ours.get(name, 1e12) >= 5000 or name.startswith("d4_"):
        heavy.setdefault(cell(name), []).append((name, body))
        continue
    if cur and size + len(body) > LIMIT:
        bundles.append(cur); cur, size = [], 0
    cur.append((name, body)); size += len(body)
if cur:
    bundles.append(cur)
for c, members in heavy.items():
    if not c.startswith("d4_"):
        bundles.append(members)
for c, members in heavy.items():
    if c.startswith("d4_"):
        bundles.append(members)
index = []
for i, b in enumerate(bundles, 1):
    head = f"// bundle {i:02d} of {len(bundles)}: " + ", ".join(n for n, _ in b) + "\n"
    text = head + "\n".join(body for _, body in b)
    assert len(text) <= 50000, (i, len(text))
    open(os.path.join(here, "bundles", f"{i:02d}.magma"), "w").write(text)
    index.append({"bundle": i, "bytes": len(text), "scripts": [n for n, _ in b]})
json.dump(index, open(os.path.join(here, "bundles", "index.json"), "w"), indent=1)
for e in index:
    print(f"{e['bundle']:02d}  {e['bytes']:6d} B  {len(e['scripts'])} scripts: {', '.join(e['scripts'])}")
