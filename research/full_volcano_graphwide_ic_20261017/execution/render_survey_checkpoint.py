"""Presentation only: plot native checkpoint statistics without recomputing them.

Usage: NATIVE_SUMMARY_DIR CHECKPOINT_DIR FRESH_REPORT RUN_ID
"""
import json
import sys
from pathlib import Path
import typst

summary_dir, checkpoint, out = map(Path, sys.argv[1:4])
run_id = sys.argv[4]
doc = json.loads((summary_dir / "summary.json").read_text())
inv = json.loads((checkpoint / "inventory.json").read_text())
assert inv["ok"] and inv["outcome_totals"] == doc["outcome_totals"]
out.mkdir(parents=True, exist_ok=False)
classes = [("two", "icv1-fp57-t536870914-6cabb0c9"),
           ("three", "icv1-fp57-tm536870941-e6d0c766"),
           ("mixed", "icv1-fp59-t1073741954-4dd0b1ef")]
plot_rows = [r for r in doc["groups"] if r["conductor_0_means_all"] > 0 and
             r["arm"] == "projected-64-mitm" and "oracle_call_relation_yield" in r["metrics"]]
(out / "plot-data.json").write_text(json.dumps({"run_id":run_id,"groups":plot_rows,
    "graphs":inv["graphs"],"outcome_totals":inv["outcome_totals"],"states":inv["states"]},indent=2)+"\n")
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="740" viewBox="0 0 1100 740">',
       '<rect width="1100" height="740" fill="white"/><g font-family="sans-serif" fill="#18364c">',
       '<text x="20" y="28" font-size="22" font-weight="bold">Repaired survey: all 51 stratum samples, projected size-64 MITM</text>',
       '<text x="20" y="51" font-size="13">Blue: equations/source trial. Orange: equations/actual oracle call. Native means and observed min–max; not confidence intervals.</text>',
       '<text x="20" y="71" font-size="12">Ordinal conductor axis; each point is a stratum/outcome group. n counts target attempts, not independent vertices.</text>']
tables = {}
for ci,(graph,slug) in enumerate(classes):
    base=112+202*ci
    rows=[r for r in plot_rows if r["graph"]==graph]
    svg.append(f'<text x="65" y="{base-15}" font-size="14" font-weight="bold">{slug}</text>')
    for tick in [0,.5,1]:
        y=base+103-tick*96
        svg += [f'<line x1="65" x2="1060" y1="{y}" y2="{y}" stroke="#dce5ea"/>',f'<text x="29" y="{y+4}" font-size="12">{tick:.1f}</text>']
    table=[]
    for i,r in enumerate(rows):
        x=81+i*960/max(1,len(rows)-1)
        for key,color,dx in [("trial_relation_yield","#256f9a",-1.5),("oracle_call_relation_yield","#bd6824",1.5)]:
            m=r["metrics"][key]
            svg += [f'<line x1="{x+dx}" x2="{x+dx}" y1="{base+103-m["max"]*96}" y2="{base+103-m["min"]*96}" stroke="{color}" stroke-width="2"/>',
                    f'<circle cx="{x+dx}" cy="{base+103-m["mean"]*96}" r="3" fill="{color}"/>']
        f=r["conductor_0_means_all"]
        svg += [f'<text transform="translate({x},{base+123}) rotate(35)" font-size="11">f={f}</text>',f'<text x="{x}" y="{base+159}" text-anchor="middle" font-size="10">n={r["records"]}</text>']
        cost=r["metrics"]["cold_gae_lower_bound"]
        cells = [f, r["status"], r["records"], f'{cost["min"]:,.0f}',
                 f'{cost["mean"]:,.1f}', f'{cost["max"]:,.0f}',
                 f'{r["metrics"]["trial_relation_yield"]["mean"]:.6f}',
                 f'{r["metrics"]["oracle_call_relation_yield"]["mean"]:.6f}']
        table.append(', '.join(f'[{v}]' for v in cells)+',')
    tables[graph]='\n'.join(table)
svg += [f'<text x="20" y="729" font-size="11">Checkpoint {run_id}: {inv["vertices_audited"]} / {inv["vertices_planned"]} vertices. All strata sampled; full-vertex distributions incomplete.</text>','</g></svg>']
(out/"all-strata-yield.svg").write_text('\n'.join(svg)+'\n')
coverage='\n'.join(', '.join(f'[{g[k]}]' for k in ['graph','vertices_audited','vertices_planned','strata_audited','strata_planned'])+',' for g in inv['graphs'])
o=doc['outcome_totals']
source=r'''#set page(paper:"a4",margin:13mm)
#set text(size:8.5pt)
= Repaired volcano survey checkpoint
''' + f'''Run {run_id}. {inv['vertices_audited']:,} audited vertices, {doc['measured_records']:,} accepted attempts.
#image("all-strata-yield.svg",width:100%)
#table(columns:(1fr,1fr,1fr,1fr,1fr),inset:3pt,
 [Graph], [Audited vertices], [Planned vertices], [Audited strata], [Planned strata],
{coverage})

Observed outcomes: {o.get('verified',0):,} verified recoveries,
{o.get('exhausted',0):,} exhausted attempts, {o.get('timeout',0):,} timeouts.
Native source trials {doc['relation_trials']:,}; actual decomposition calls
{doc['decomposition_calls']:,}; equations {doc['relations_found']:,}.
''' + r'''
These pooled counts span different subgroups/bases and are accounting totals,
not a common success probability. Multiple outcomes and all missing/censored
observations remain distinct. Stratum samples are not complete populations.

The environment returned without the previous worker. One in-flight attempt
was preserved and marked interrupted by the unchanged driver's zero-job stop;
it is excluded from accepted means. Exact termination cause was not observed.
No sealed session, binary, queue or mathematical budget was changed.

The native inventory checks record seals, receipt/file hashes, ledger/queue
identities and exact existing successful/exhausted replay coverage. This is
receipt binding, not fresh or independent solver replay. The entire stopped
campaign, including the unaccepted interruption, accompanies the inventory.

#pagebreak()
= Stratum distributions: recorded size-64 MITM
Trial yield includes repeated queries skipped before the oracle; actual-call
yield excludes them. Identity skips occur before trial increment. Zero or
unknown denominators have no yield. Timeout final counters are unknown, not
zero. GAE with unpriced work is a lower bound; exhaustion costs an attempt,
not a successful DLP. Native summary/CSV retain all other arms, phases and
present-value n. The tables below use only recorded size-64 MITM reports.
'''
for index,(graph,slug) in enumerate(classes):
    if index==2:source+='\n#pagebreak()\n'
    source+=f'\n== {slug}\n#set text(size:7pt)\n#table(columns:(.65fr,.9fr,.45fr,.9fr,.9fr,.9fr,1fr,1fr),inset:2.5pt,\n[Conductor], [Outcome], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],\n'+tables[graph]+')\n#set text(size:8.5pt)\n'
source+=r'''
Protocol: actual largest prime subgroups; two public targets/model, one
alternating round, zero warmups, bases 16/64, raw-base control, MITM/subtraction,
rho.negation, BSGS and rho A/A. 4,096 source trials, 3-second child watchdog.
Cold setup/projection/verification charged; unpriced work remains explicit.
The full eight-target/five-round/million-trial matrix remains unexecuted.
''' + f'''
Frozen release runner {doc['binary_sha256']}. Source 2bf7700b7, unchanged
queue and driver. Earlier dev99, stopped release70, control/pilot and original
C finite-field diagnostic are separate datasets; repeated workloads and
different binaries are never pooled here. No performance improvement follows
from this checkpoint's larger sample count.
''' + r'''
Wall observations are unqualified and descriptive. Every capsule retains its
host/environment observation; any resumed runtime epoch must remain identifiable.
No admitted runtime, external fastest-method, scaling or sub-rho claim is made.
Exact representation registration, transported targets, composite PH,
independent arithmetic replay, original A--E full graphs and large-extension
finite-field backends remain open. Original C's r=7 k=1 diagnostic is separate.
Complete graph topology is for the declared fields and finite prime degrees,
not every other prime or infinitely many compositions.
'''
(out/'report.typ').write_text(source)
typst.compile(str(out/'report.typ'),output=str(out/'report.pdf'))
(out/'figure.typ').write_text('#set page(width:1100pt,height:740pt,margin:0pt)\n#image("all-strata-yield.svg",width:100%)\n')
typst.compile(str(out/'figure.typ'),output=str(out/'figure.pdf'))
