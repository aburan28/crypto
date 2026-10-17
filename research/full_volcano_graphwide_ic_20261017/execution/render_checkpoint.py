"""Presentation only: native summary and integrity inventory; no new results."""
import html
import json
import sys
from pathlib import Path
import typst

summary_dir, checkpoint, out = map(Path, sys.argv[1:4])
doc = json.loads((summary_dir / "summary.json").read_text())
inv = json.loads((checkpoint / "inventory.json").read_text())
meta = json.loads(Path(sys.argv[4]).read_text()) if len(sys.argv) > 4 else None
out.mkdir(exist_ok=False, parents=True)
(out / "plot-data.json").write_text(json.dumps({"summary": doc, "inventory": inv}, indent=2) + "\n")
classes = [("two", "icv1-fp57-t536870914-6cabb0c9"),
           ("three", "icv1-fp57-tm536870941-e6d0c766"),
           ("mixed", "icv1-fp59-t1073741954-4dd0b1ef")]
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="740" viewBox="0 0 1100 740">',
       '<rect width="1100" height="740" fill="white"/><g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="28" font-size="22" font-weight="bold">All 51 conductor strata: measured size-64 MITM yields</text>',
       '<text x="20" y="51" font-size="13">Blue: relations per source trial. Orange: relations per actual oracle call. Points: native mean; bars: observed min–max.</text>',
       '<text x="20" y="71" font-size="12">Conductor order is ordinal on each x-axis; n is measured target attempts, not independent vertices or a confidence interval.</text>']
tables = {}
for ci, (graph, slug) in enumerate(classes):
    base = 112 + 202 * ci
    rows = [r for r in doc["groups"] if r["graph"] == graph and r["conductor_0_means_all"] > 0 and r["arm"] == "projected-64-mitm"]
    svg.append(f'<text x="65" y="{base-15}" font-size="14" font-weight="bold">{slug}</text>')
    for tick in [0, .5, 1]:
        y = base + 103 - tick * 96
        svg += [f'<line x1="65" x2="1060" y1="{y}" y2="{y}" stroke="#dce5ea"/>',
                f'<text x="29" y="{y+4}" font-size="12">{tick:.1f}</text>']
    prior = [None, None]
    table = []
    for i, row in enumerate(rows):
        x = 81 + i * 960 / max(1, len(rows) - 1)
        f = row["conductor_0_means_all"]
        for si, (key, color) in enumerate([("trial_relation_yield", "#256f9a"), ("oracle_call_relation_yield", "#bd6824")]):
            m = row["metrics"][key]
            y = base + 103 - m["mean"] * 96
            dx = -1.5 if si == 0 else 1.5
            if prior[si]:
                px, py = prior[si]
                svg.append(f'<line x1="{px}" y1="{py}" x2="{x+dx}" y2="{y}" stroke="{color}" stroke-width="1.5"/>')
            svg += [f'<line x1="{x+dx}" x2="{x+dx}" y1="{base+103-m["max"]*96}" y2="{base+103-m["min"]*96}" stroke="{color}" stroke-width="2"/>',
                    f'<circle cx="{x+dx}" cy="{y}" r="3" fill="{color}"/>']
            prior[si] = (x+dx, y)
        svg += [f'<text transform="translate({x}, {base+123}) rotate(35)" font-size="11">f={f}</text>',
                f'<text x="{x}" y="{base+159}" text-anchor="middle" font-size="10">n={row["records"]}</text>']
        cost = row["metrics"]["cold_gae_lower_bound"]
        trial = row["metrics"]["trial_relation_yield"]
        actual = row["metrics"]["oracle_call_relation_yield"]
        table.append(f'[{f}], [{row["status"]}], [{row["records"]}], [{cost["min"]:,.0f}], [{cost["mean"]:,.1f}], [{cost["max"]:,.0f}], [{trial["mean"]:.6f}], [{actual["mean"]:.6f}],')
    tables[graph] = "\n".join(table)
svg += ['<text x="20" y="729" font-size="11">Run 20261017T0337Z: 99 audited vertices / 24,823 planned. All strata covered by samples; entire-graph distributions remain incomplete.</text>', '</g></svg>']
if meta:
    svg[-2] = (f'<text x="20" y="729" font-size="11">Run {html.escape(meta["run_id"])}: '
               f'{inv["vertices_audited"]} audited vertices / {inv["vertices_planned"]:,} planned. '
               'All strata sampled; entire-graph distributions incomplete. One failed capsule excluded.</text>')
(out / "all-strata-yield.svg").write_text("\n".join(svg) + "\n")
coverage = "\n".join(f'[{r["graph"]}], [{r["vertices_audited"]}], [{r["vertices_planned"]}], [{r["strata_audited"]}], [{r["strata_planned"]}], [{r["max_audited_conductor"]}],' for r in inv["graphs"])
source = r'''#set page(paper: "a4", margin: 13mm)
#set text(size: 8.3pt)
#set heading(numbering: none)
= Volcano IC: all strata, 99-vertex checkpoint
Run 20261017T0337Z. Counted diagnostic; no admitted runtime or speedup claim.

#image("all-strata-yield.svg", width: 100%)

#table(columns: (1fr,1fr,1fr,1fr,1fr,1fr), inset: 3pt,
 [Graph key], [Audited vertices], [Planned vertices], [Audited strata], [Planned strata], [Largest End(E) conductor],
''' + coverage + r''')

*Observed:* 1,584 measured attempts: 1,065 verified recoveries, 422 exhausted
budgets, 97 watchdog timeouts. All 51 conductor strata have measured vertices;
only 99 of 24,823 vertices have been audited. Stratum coverage is not complete
graph coverage. Maximal sampled conductors are 8192, 6561 and 1728.

*Integrity:* every standard audit receipt passes; successful/exhausted
deterministic reports replayed when collected. The stopped checkpoint inventory
rechecks every receipt/file hash, all 1,584 record seals, exact replay populations
and ledger/queue identities. It is a receipt-binding check, not a fresh or
independent solver replay. Exact inventory and raw frozen files accompany this report.

#pagebreak()
= Yield definitions and complete stratum tables
*Correction to earlier wording:* the earlier reports called 301,715 source
trials “attempted decompositions.” Trials include duplicate queries skipped
before the oracle. Their recorded values and historical PDFs are unchanged.
The new summary retains both denominators: 2,223,204 source trials and 2,075,145
actual decomposition calls, producing 57,176 relations across 990 IC reports.
For these guarded classic-walk arms, calls equal trials minus repeated queries
skipped. Identity skips occur before trial increment. All 594 two-summand MITM
records independently agree with their lookup counters. A zero/unknown call
denominator has no yield; timeout counters are never assumed zero.

The legacy field relation-yield remains a per-trial alias. Summary schema v2
adds explicit trial-relation-yield and oracle-call-relation-yield, numerator,
denominator, sample counts and definitions. Global totals across distinct
subgroups/bases are descriptive accounting, not a common success probability.

The tables show only size-64 projected MITM, with every conductor stratum.
GAE values are lower bounds with unpriced native work. Exhausted rows are
attempt costs, not successful DLP costs. Means use present observations and
the explicit n; multiple vertices per stratum are sampled, not a census.
Native CSV/JSON contain every arm, phase, status, range, matrix rank/dependency,
memory and raw elapsed observation. Exact generator representations remain
explicit and unregistered where noted; no new curve-name alias is inferred.
'''
for index, (graph, slug) in enumerate(classes):
    if index == 2:
        source += "\n#pagebreak()\n"
    source += f'\n== {slug}\n'
    source += '#set text(size: 7pt)\n#table(columns: (.6fr,.9fr,.4fr,.9fr,.9fr,.9fr,1fr,1fr), inset: 2.5pt,\n[Conductor], [Status], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],\n' + tables[graph] + ')\n#set text(size: 8.3pt)\n'
source += r'''
*Protocol and provenance:* projected bases 16/64 with MITM/subtraction,
unchanged raw-base MITM, rho.negation, textbook BSGS and rho A/A. Two public
targets per vertex, one alternating round, zero warmups, target seed 20261017,
order seed 1, 4,096 source trials and a 3-second child watchdog. Actual largest
prime subgroups are retained, including the 45-bit subgroup. Cold solver setup,
projection and verification are charged. All original measured sessions remain
byte-frozen. Source commits 4eecf6b53 and ee54ad603; complete stopped campaign,
native source, queue, ledger, receipts and frozen runner are in the archive.

*Host and limits:* Intel Xeon Platinum 8573C, Linux x86-64, four-CPU quota,
16 GiB RAM. Dev-profile ecbench hash
#text(size: 6.5pt)[1869ca2e7be271794df8128c08a6dca0201038e2c81a3e0ce6434f7ec0d8379d].
All wall observations are L0/debug and unqualified. No builds/test suites ran
during the resumed measured block; it was stopped before the release build.
Raw process time includes orchestration; no optimized speedup is inferred.

The full million-trial/eight-target/five-round/four-size matrix, original
A--E full graph request, all-vertex completion, exact representation registry,
transported-target benchmarks, composite PH, actual finite-field IC and
independent implementation replay remain open. Exact embedding degrees
30,875; 3,953,777,913,734; and 2,126 exceed the implemented MOV teaching backend.
This is an implementation gap, not proof that mathematical reductions do
not exist. A separately named release-profile survey must pass matched count
checks and audits before launch. Canonical fastest-method claims are unchanged.
'''
if meta:
    # Presentation-only substitutions from the native inventory and summary.
    # No result is recomputed here; old published report inputs remain usable.
    source = source.replace("99-vertex checkpoint", f'{inv["vertices_audited"]}-vertex release checkpoint')
    source = source.replace("20261017T0337Z", meta["run_id"])
    source = source.replace("1,584", f'{doc["measured_records"]:,}').replace("1,065", f'{doc["outcome_totals"]["verified"]:,}')
    source = source.replace("422 exhausted", f'{doc["outcome_totals"]["exhausted"]:,} exhausted')
    source = source.replace("97 watchdog", f'{doc["outcome_totals"]["timeout"]:,} watchdog')
    source = source.replace("99 of 24,823", f'{inv["vertices_audited"]} of {inv["vertices_planned"]:,}')
    source = source.replace("every standard audit receipt passes", "every included standard audit receipt passes")
    source = source.replace("Exact inventory and raw frozen files accompany this report.",
        "The failed three-v10-attempt-001 capsule (16 raw records) and failed audit are retained separately in the archive; they are excluded from these accepted summaries. Its record 10 combines timeout status with a checked scalar. No old record was relabeled. A watchdog recording repair is tested separately.")
    source = source.replace("2,223,204", f'{doc["relation_trials"]:,}').replace("2,075,145", f'{doc["decomposition_calls"]:,}')
    source = source.replace("57,176", f'{doc["relations_found"]:,}').replace("990 IC reports", f'{doc["relation_observations"]} IC reports')
    source = source.replace("All 594", f'All {doc["mitm_call_lookup_checks"]}')
    source = source.replace("Source commits 4eecf6b53 and ee54ad603", f'Source commit {meta["source_commit"]}')
    start = source.index("Dev-profile ecbench hash")
    end = source.index("The full million-trial", start)
    source = source[:start] + f'''Release-profile ecbench hash
#text(size: 6.5pt)[{doc["binary_sha256"]}].
All wall observations are unqualified. No optimized runtime comparison is
admitted. This worker halted on its failed audit after 432 elapsed seconds;
that invocation includes orchestration and replay, not only solver work.
The historical dev checkpoint remains a separate dataset. Neither its repeated
targets nor its timings are pooled here.

''' + source[end:]
    source = source.replace("A separately named release-profile survey must pass matched count\nchecks and audits before launch.", "A repaired runner must pass fresh controls and completed-report matching before a separate new campaign launches.")
(out / "report.typ").write_text(source)
typst.compile(str(out / "report.typ"), output=str(out / "report.pdf"))
(out / "figure.typ").write_text('#set page(width: 1100pt, height: 740pt, margin: 0pt)\n#image("all-strata-yield.svg", width: 100%)\n')
typst.compile(str(out / "figure.typ"), output=str(out / "figure.pdf"))
print(out / "report.pdf")
