"""Presentation of native QA summaries and profile-match receipts only."""
import json
import sys
from pathlib import Path
import typst

summary_path, control_path, pilot_path, out = map(Path, sys.argv[1:])
doc = json.loads(summary_path.read_text())
control = json.loads(control_path.read_text())
pilot = json.loads(pilot_path.read_text())
observations = [json.loads(line) for line in summary_path.with_name("observations.jsonl").read_text().splitlines()]
assert control["ok"] and pilot["ok"]
out.mkdir(parents=True, exist_ok=False)
(out / "plot-data.json").write_text(json.dumps({"summary":doc,"control_match":control,"pilot_match":pilot},indent=2)+"\n")
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="330" viewBox="0 0 1000 330"><rect width="1000" height="330" fill="white"/><g font-family="sans-serif" fill="#18364c">',
       '<text x="25" y="30" font-size="22" font-weight="bold">Repaired runner: audited control and three-root pilot</text>',
       '<text x="25" y="56" font-size="14">60 attempts, two public targets per model, one round, zero warmups. Counts only; no admitted runtime comparison.</text>']
for i, (status,color) in enumerate([("verified","#256f9a"),("exhausted","#bd6824"),("timeout","#84949b")]):
    n = doc["outcome_totals"][status]
    y = 100 + i*55
    svg += [f'<text x="25" y="{y+19}" font-size="17">{status}</text>',f'<rect x="170" y="{y}" width="{n*16}" height="28" fill="{color}"/>',f'<text x="{183+n*16}" y="{y+20}" font-size="16">{n}</text>']
svg += ['<text x="25" y="295" font-size="15">Completed pairs match baseline: control 12, pilot 46. Two pilot pairs remain censored, not equal-count comparisons.</text>', '</g></svg>']
(out / "qa-outcomes.svg").write_text("\n".join(svg)+"\n")
rows = ""
for r in doc["groups"]:
    if r["conductor_0_means_all"] != 0:
        continue
    m = r["metrics"].get("cold_gae_lower_bound")
    values = [str(m["n"]),f'{m["min"]:,.0f}',f'{m["mean"]:,.2f}',f'{m["max"]:,.0f}'] if m else ["0","missing","missing","missing"]
    cells = [r["graph"],r["arm"],r["status"],str(r["records"])] + values
    rows += ", ".join(f'[{v}]' for v in cells) + ",\n"
source = r'''#set page(paper:"a4",margin:13mm)
#set text(size:9pt)
= Watchdog repair: control and pilot evidence
17 October 2026. New campaign data, separate from historical dev/release runs.
#image("qa-outcomes.svg",width:100%)

The earlier release survey halted after 70 audited vertices on a failed
standard audit. One late child message carried a checked scalar while the
parent labeled the attempt timeout. Its failed 16-record capsule is retained;
no historical record or auditor rule was changed.

The repaired runner makes watchdog censoring authoritative. Timed-out attempts
have no accepted scalar, success check, final cost/counters or solve timing.
Any parsed late message remains in the sealed diagnostic error. Solver
arithmetic, public targets, factor bases and budgets are unchanged.

Fresh small control: 12 attempts. Fresh three-root pilot: 48 attempts. Every
session passes the standard audit; successful and exhausted deterministic
reports reproduce. Native comparison to the earlier release baseline checks
method identities, exact outcomes, scalars, GAE, phases, field operations,
counters and factor bases. It finds 12 matching completed control pairs and
46 matching completed pilot pairs. Two censored pilot pairs have no completed
reports to compare. Elapsed values are excluded. Same-implementation replay
and profile matching are not independent implementation validation.

The frozen runner hash is
#text(size:7pt)[35d89073f4d666b094c213a1d31b0c819b63600123e067c111af3c947d7609a9].
The driver hash is
#text(size:7pt)[b8b489c7225c64ed492ffce39286c995f5ca820960956cf706c0d3110f5bdd2d].
Clean source 2bf7700b7 contains the repair. The native build source/protocol
and input/binary hashes are frozen in each campaign. The separate full serial
release library check passes 4,236 tests, zero failures, 117 existing ignored.
The previously failing parallel Joux--Vitse case passes serially; its parallel
cause remains unresolved. No solver budget or exact-answer assertion changed.

Host: Xeon Platinum 8573C, Linux x86-64, four-CPU quota, 16 GiB memory limit.
Release optimization; wall observations remain unqualified. No external
comparison, fastest-method claim or speedup is made. The methods retain actual
largest prime orders 61751, 19768889568671 and 4253 in the pilot; the small
control has order 7 over p=113. No smaller subgroup replaced the 45-bit case.

Passing these gates permits a separately frozen 24,823-vertex bounded survey,
not the full eight-target/five-round/million-trial matrix. Planning is not
execution. Its actual launch and current coverage require the separate launch
receipt and live ledger. Original A--E completion and large-extension
finite-field IC remain open. Source protocols, raw capsules and native
comparison receipts are in the accompanying archive and study execution tree.

#pagebreak()
= Native per-arm GAE distributions
#set text(size:6.7pt)
#table(columns:(.65fr,1.5fr,.9fr,.35fr,.35fr,.85fr,.9fr,.85fr),inset:2pt,
 [Graph], [Arm], [Outcome], [Attempts], [Cost n], [Min GAE], [Mean GAE], [Max GAE],
''' + rows + r''')
#set text(size:8pt)
GAE with unpriced work is a lower bound. Exhaustion is an attempted-budget
cost, not a successful DLP cost. Timeout counters/costs are missing, not zero.
Each outcome stays separate, with present-value n. Raw phase/native counters,
memory and unqualified elapsed observations remain in the native summary and
CSV. Repeated baseline and candidate targets are not pooled.
'''
source += '\n== Censored observations (one row per attempt)\n#set text(size:7pt)\n#table(columns:(.5fr,.5fr,1.5fr,.7fr,1fr,1fr),inset:3pt,\n[Graph], [Arm], [Record id], [Final GAE], [Raw process s], [Peak RSS KiB],\n'
for r in observations:
    if r["status"] == "timeout":
        cells = [r["graph"],r["arm"],r["record_id"],"missing",r["metrics"].get("process_seconds","missing"),r["metrics"].get("peak_rss_kib","missing")]
        source += ", ".join(f'[{v}]' for v in cells) + ",\n"
source += ')\n#set text(size:8pt)\nThese process observations are not completed solve times. The final-cost table above contains present reports only; censored costs stay missing.\n'
(out / "report.typ").write_text(source)
typst.compile(str(out / "report.typ"),output=str(out / "report.pdf"))
(out / "figure.typ").write_text('#set page(width:1000pt,height:330pt,margin:0pt)\n#image("qa-outcomes.svg",width:100%)\n')
typst.compile(str(out / "figure.typ"),output=str(out / "figure.pdf"))
