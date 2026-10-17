"""Presentation/export of recorded native values only; no solver or aggregation.

Usage: ROOT_RESULT PREFIX_RESULT APPLICABILITY FRESH_REPORT LOCAL_CERT STRATUM_STATS
"""
import csv
import html
import json
import sys
from pathlib import Path
import typst

root_path, prefix_path, inventory_path, out, cert_path, stratum_path = map(Path, sys.argv[1:])
root = json.loads(root_path.read_text())
prefix = json.loads(prefix_path.read_text())
inventory = json.loads(inventory_path.read_text())
local = json.loads(cert_path.read_text())
stratum_stats = json.loads(stratum_path.read_text())
assert prefix["prefix"]["sha256"] == local["db_sha256"] == stratum_stats["db_sha256"]
out.mkdir(parents=True, exist_ok=False)
stats = prefix["summary"]
(out / "plot-data.json").write_text(json.dumps({"prefix": stats, "root": root["summary"], "applicability": inventory,
    "local_certificate": {k:v for k,v in local.items() if k != "rows"}, "stratum_stats": stratum_stats}, indent=2) + "\n")

def flatten(value, path="", result=None):
    result = {} if result is None else result
    if isinstance(value, dict):
        for key, item in value.items():
            flatten(item, f"{path}/{key}" if path else key, result)
    else:
        result[path] = json.dumps(value) if isinstance(value, list) else value
    return result

flat = []
for observation in prefix["rows"]:
    certified = local["rows"][observation["graph_vertex"]]
    flat.append(flatten(dict(observation, local_conductor=certified["conductor"], local_conductor_status=certified["status"])))
fields = sorted({key for row in flat for key in row})
with (out / "prefix-observations.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(flat)

svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="380" viewBox="0 0 1000 380">',
       '<rect width="1000" height="380" fill="white"/><g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="28" font-size="22" font-weight="bold">Original C prefix: finite-field work per verified target</text>',
       f'<text x="20" y="52" font-size="14">n={stats["verified"]}; all stored prefix vertices; two public target searches per model. Not the full C class.</text>',
       '<text x="20" y="74" font-size="13">Each row has its own scale. Orange dot: native arithmetic mean. Blue line: observed min–max. No confidence interval.</text>']
labels = {"queries": "Destination queries", "candidate_pair_tests": "Candidate pair products", "relation_count": "Collected relations", "matrix_entry_ops": "Matrix entry operations"}
for i, (key, label) in enumerate(labels.items()):
    m = stats["metrics"][key]
    y = 120 + i * 65
    svg.append(f'<text x="20" y="{y}" font-size="15">{label}</text>')
    if m["n"]:
        scale = 430 / max(1, m["max"])
        svg += [f'<line x1="{330+m["min"]*scale}" x2="{330+m["max"]*scale}" y1="{y-5}" y2="{y-5}" stroke="#256f9a" stroke-width="5"/>',
                f'<circle cx="{330+m["mean"]*scale}" cy="{y-5}" r="6" fill="#bd6824"/>',
                f'<text x="790" y="{y}" font-size="14">{m["min"]} / {m["mean"]:.3f} / {m["max"]}</text>']
    else:
        svg.append(f'<text x="340" y="{y}" font-size="14">No verified observations</text>')
svg += ['</g></svg>']
(out / "prefix-work.svg").write_text("\n".join(svg) + "\n")

groups = [g for g in stratum_stats["groups"] if g["conductor"] is not None and g["status"] == "verified"]
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="650" viewBox="0 0 1100 650">',
       '<rect width="1100" height="650" fill="white"/><g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="28" font-size="22" font-weight="bold">Original C: locally certified conductor strata</text>',
       f'<text x="20" y="52" font-size="14">{local["strata_exact_observed"]} observed strata; {stratum_stats["records_exact_conductor"]} target records with exact local conductor. Not all 165 C strata.</text>',
       '<text x="20" y="74" font-size="13">Orange: native mean. Blue: observed min–max. n is target records; two searches per stored vertex. Ordinal conductor axis.</text>']
for panel, (metric, label) in enumerate([("relations_per_query", "Collected equations per destination query"), ("matrix_entry_ops", "Matrix entry operations per target")]):
    top = 115 + panel * 260
    axis_max = max([1] + [g["metrics"][metric]["max"] for g in groups])
    svg.append(f'<text x="35" y="{top}" font-size="15" font-weight="bold">{label}</text>')
    for tick in [0, axis_max / 2, axis_max]:
        y = top + 140 - tick * 115 / axis_max
        svg += [f'<line x1="70" x2="1065" y1="{y}" y2="{y}" stroke="#dce5ea"/>', f'<text x="25" y="{y+4}" font-size="12">{tick:g}</text>']
    for index, g in enumerate(groups):
        x = 80 + index * 975 / max(1, len(groups)-1)
        m = g["metrics"][metric]
        svg += [f'<line x1="{x}" x2="{x}" y1="{top+140-m["max"]*115/axis_max}" y2="{top+140-m["min"]*115/axis_max}" stroke="#256f9a" stroke-width="3"/>',
                f'<circle cx="{x}" cy="{top+140-m["mean"]*115/axis_max}" r="3" fill="#bd6824"/>',
                f'<text transform="translate({x-2},{top+157}) rotate(50)" font-size="10">{g["conductor"]}</text>',
                f'<text x="{x}" y="{top+230}" text-anchor="middle" font-size="10">n={g["records"]}</text>']
svg += ['</g></svg>']
(out / "local-strata.svg").write_text("\n".join(svg) + "\n")

pipeline = '''<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="180" viewBox="0 0 1000 180">
<rect width="1000" height="180" fill="white"/>
<g font-family="sans-serif" fill="#18364c" font-size="15">
<text x="20" y="25" font-size="20" font-weight="bold">Executed diagnostic route: original C, largest prime subgroup r=7</text>
<path d="M220 90H270M490 90H540M760 90H810" stroke="#256f9a" stroke-width="3"/>
<g fill="#edf4f8" stroke="#256f9a"><rect x="20" y="53" width="200" height="76" rx="5"/><rect x="270" y="53" width="220" height="76" rx="5"/><rect x="540" y="53" width="220" height="76" rx="5"/><rect x="810" y="53" width="170" height="76" rx="5"/></g>
<text x="38" y="80">EC generator / target</text><text x="38" y="107">cofactor projection</text>
<text x="288" y="80">Reduced Tate / Miller</text><text x="288" y="107">shifted divisor, Fp</text>
<text x="558" y="80">Projected prime base</text><text x="558" y="107">two-summand relations</text>
<text x="826" y="80">Gauss mod 7</text><text x="826" y="107">Fp + EC checks</text>
<text x="20" y="162">Two factor-base logs and one target log are unknowns. No rho/BSGS fallback. This is not NFS or integer-smoothness IC.</text>
</g></svg>'''
(out / "route.svg").write_text(pipeline)

def row(cells):
    return ", ".join(f"[{c}]" for c in cells) + ",\n"

app = "".join(row([r["case"], r["largest_prime_subgroup"], r["embedding_degree"], r["dense_u64_coefficient_bytes_per_element"]]) for r in inventory["rows"])
coverage = "".join(row([label, doc["summary"]["requested"], doc["summary"]["verified"], doc["summary"]["failed"]]) for label, doc in [("Original source model", root), ("Stored prefix", prefix)])
m = stats["end_to_end_ns"]
phase_names = {"setup": "Setup", "pairing": "Pairing", "factor_base_relations_matrix_and_ff_verification": "FF solve/checks", "ec_verification": "EC check"}
phase_rows = ""
for key, value in stats["phase_operations"].items():
    phase, primitive = key.split("/", 1)
    phase_rows += row([phase_names[phase], primitive.replace("_", "-"), value["n"], value["min"], f'{value["mean"]:.3f}' if value["mean"] is not None else "missing", value["max"]])
yield_rows = "".join(row([key.replace("_", "-"), v["n"], v["min"], f'{v["mean"]:.6f}' if v["mean"] is not None else "missing", v["max"]]) for key, v in stats["metrics"].items() if key in ("relations_per_query", "query_hit_yield"))

source = r'''#set page(paper: "a4", margin: 13mm)
#set text(size: 9pt)
= Original C: finite-field reduction diagnostic
Source model icv1-fp61-t2-fefb1f9c. Public synthetic inputs, 17 October 2026.
This report is descriptive; no admitted runtime or scaling claim.

The original C model has p=1637958981148213249,
a=571607911220311176, b=29329524900995129 and
N=1637958981148213248. Native complete factorization gives
N=2^26 × 3^20 × 7. Its largest prime subgroup is 7 and embedding degree is 1.
The 61-bit field is not a 61-bit prime-order DLP.

#image("route.svg", width: 100%)

The shifted-divisor Miller implementation computes standard reduced Tate
evaluations. Finite tests check all six nonzero control/C scalar inputs and
all nontrivial destination generators/targets; these are implementation tests,
not a new bilinearity theorem or independent implementation proof. Recovered
factor-base logs, the destination answer and the original EC equation are
checked. The solver receives neither planted scalars nor factor-base logs.

== Original-input arithmetic inventory (not five pairing runs)
#table(columns: (.4fr,1.4fr,1.6fr,1.8fr), inset: 4pt,
 [Case], [Largest prime r], [Exact embedding degree], [Dense u64 bytes / field element],
''' + app + r''')
Dense coefficient sizes describe that specific layout only. They are not
universal lower bounds, impossibility proofs or estimates for other algorithms.
The A/B/D/E finite-field routes have not been executed by this diagnostic.

== Coverage and retained limitations
#table(columns: (2fr,1fr,1fr,1fr), inset: 4pt,
 [Population], [Requested targets], [Verified], [Failed],
''' + coverage + r''')
The prefix contains 1,869 stored vertices, 1,024 expanded vertices and 4,903
replayed maps. Its frontier is not exhausted; global stratum populations remain
NOT-ESTABLISHED. The separate conservative local certificate is described on
page 4. Original C's expected class has 1,934,901,248 vertices.
The source-model and prefix populations stay separate, including any repeated
target. Two public searches may select the same subgroup point; they are not
independent random samples. Other prime-power components are not measured.

#pagebreak()
= Recorded prefix work and elapsed observations
#image("prefix-work.svg", width: 100%)

Relation yield has two explicit meanings: number of collected equations per
destination query (several per query are allowed), and fraction of queries
with any relation. Identity queries are skipped separately. Relations satisfy
the checked destination equation before matrix insertion. Full rank is three.
#table(columns: (2fr,.5fr,1fr,1fr,1fr), inset: 4pt,
 [Metric], [n], [Min], [Mean], [Max],
''' + yield_rows + r''')

Native per-target end-to-end elapsed interval: setup/order factorization,
point search/projection, pairing, cold factor base, collection/matrix and
destination/EC verification. Input file loading, batch serialization and
output I/O are outside that per-target interval; complete process resource
logs are separate. All observations are unqualified; there is no external
reference or fastest-method claim.
#table(columns: (.5fr,1fr,1fr,1fr), inset: 4pt,
 [n], [Min ns], [Mean ns], [Max ns],
''' + row([m["n"], m["min"], f'{m["mean"]:.3f}' if m["mean"] is not None else "missing", m["max"]]) + r''')

The archive retains both executions, exact binary, source, protocol, fixed
prefix DB and input hashes, native summaries and per-target data. The native
comparison removes only elapsed time fields, then compares every remaining
value, including source/model identities, outcomes, relations and operation
ledgers. It is same-binary reproduction, not independent implementation replay.
See root-replay-match.json and prefix-replay-match.json for actual checks.

Unpriced work stays explicit: order factorization internals; square-root
internals (calls counted); allocation/comparisons; matrix inversions (matrix
entry operations counted). Explicit Fp inversions are expanded into their
multiplication/square counts: pricing those inversions again would double count.
Native EC field counts are separate from the explicit Miller/destination work.
There is no common-cost conversion or complete end-to-end GAE claim.

The three complete additional classes, full A--E graphs, all requested strata,
full benchmark matrix, exact representation registration and large-extension
finite-field backends remain separate obligations. The original graph explorer
and fastest-method scoreboards are unchanged: this diagnostic establishes no
new complete graph, general algebraic identity or admitted speed comparison.

#pagebreak()
= Per-phase operation distributions (native summary)
#set text(size: 6.7pt)
#table(columns: (1fr,2.25fr,.5fr,.8fr,.9fr,.8fr), inset: 2pt,
 [Phase], [Primitive ledger], [n], [Min], [Mean], [Max],
''' + phase_rows + r''')
#set text(size: 8pt)
Source: #raw("execution/finite-field-k1-evidence"),
#raw("execution/finite-field-k1-summary.json"),
#raw("execution/FINITE_FIELD_K1_PROTOCOL.md") and the delivery manifest. Counts cover
present verified observations only; failed rows remain in raw data. All raw
target and primitive fields are also exported in prefix-observations.csv.
'''
source += r'''
#pagebreak()
= Local conductor bounds: partial coverage of original C
#image("local-strata.svg", width: 100%)
''' + f'''
The local certificate labels {local["vertices_exact"]} vertices across
{local["strata_exact_observed"]} conductor strata, from
{local["exact_populations"][0]["conductor"]:,} to
{local["exact_populations"][-1]["conductor"]:,}. The source model's exact local
conductor is {local["rows"][0]["conductor"]:,}. The other
{local["vertices_unknown"]} vertices retain unknown conductor;
{stratum_stats["records_unknown_conductor"]} target records remain in that separate bucket.
''' + r'''
This uses the existing ordinary-volcano floor-distance correspondence. Complete
kernel multiplicities identify known floors. Distance to known floors bounds
height above; distance to known floors or unexpanded frontier bounds it below.
Only equality of finite bounds produces an exact local height/conductor.
Unexpanded nodes are possible floors for the bound, never asserted floors.
The frozen prefix's global stratum receipt and database are unchanged.

The second labeling implementation passed all checks against the three
complete graphs at visibility cutoffs 0, 1, 17, 1024 and all. These are finite
calibrations using the same standard volcano theorem and isogeny arithmetic,
not a formal proof or independent arithmetic implementation. The full-class
population and all 165 strata remain unestablished. Local annotation is study
preparation and is not used by the DLP solver or charged to its elapsed interval.

Native per-stratum summaries keep all phases, statuses, min/max/arithmetic
means and present-value n. An unknown conductor is a null bucket, not conductor
zero or an omitted observation. No causal effect of conductor on cost is
inferred from these differing points and sample sizes.
'''
(out / "report.typ").write_text(source)
typst.compile(str(out / "report.typ"), output=str(out / "report.pdf"))
for name, width, height in [("prefix-work", 1000, 380), ("route", 1000, 180), ("local-strata",1100,650)]:
    (out / f"{name}.typ").write_text(f'#set page(width: {width}pt, height: {height}pt, margin: 0pt)\n#image("{name}.svg", width: 100%)\n')
    typst.compile(str(out / f"{name}.typ"), output=str(out / f"{name}.pdf"))
print(out / "report.pdf")
