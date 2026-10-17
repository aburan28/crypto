"""Presentation only. Native summaries provide all counts, means and ranges.

Usage: CAMPAIGN_ROOT RUN_NAME FRESH_REPORT_DIR
"""
import json
import datetime
import sys
from pathlib import Path
import typst

root, name, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
doc = json.loads((root/f'derived-{name}/summary.json').read_text())
meta = json.loads((root/name/'meta.json').read_text())
iso = json.loads((root/f'{name}-isolation.jsonl').read_text())
out.mkdir(parents=True, exist_ok=False)
started=datetime.datetime.fromtimestamp(iso['run']['before']['unix_time'],datetime.timezone.utc).isoformat()
methods = [('ec.ph_bsgs.full_point_v1', 'EC PH / BSGS', '#24698a'),
           ('ff.k1_ph_relations.full_point_v1', 'Tate + FF PH relations', '#c16b24')]
panels = [('ec_adds_and_doubles', 'Curve additions and doublings: call counts, including scalar multiplication'),
          ('explicit_fp_muls_and_squares', 'Explicit Fp multiplications and squarings: curve-field work remains separate'),
          ('ec_field_inversions_unexpanded', 'Curve-field inversions: calls counted, internals unexpanded')]
rows = [r for r in doc['strata'] if r['status']=='verified']
conductors = sorted({r['conductor'] for r in rows}, key=lambda f:(f is None,f or 0))
lookup = {(r['method'],r['conductor']):r for r in rows}
(out/'plot-data.json').write_text(json.dumps({'source_records_sha256':doc['records_sha256'],
    'native_rows':rows,'overall':doc['overall'],'run':name},indent=2)+'\n')
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="760" viewBox="0 0 1100 760">',
     '<rect width="1100" height="760" fill="white"/><g font-family="sans-serif" fill="#17354a">',
     f'<text x="24" y="28" font-size="22" font-weight="bold">Original C: full selected-point-order diagnostic — {name}</text>',
     '<text x="24" y="53" font-size="14">Native per-target means and observed min–max, not confidence intervals. Blue: EC control. Orange: finite-field route.</text>']
for pi,(key,title) in enumerate(panels):
    top=96+pi*193
    values=[r['metrics'][key] for r in rows]
    maximum=max((v['max'] for v in values),default=1) or 1
    svg.append(f'<text x="67" y="{top-9}" font-size="15" font-weight="bold">{title}</text>')
    for fraction in [0,.5,1]:
        y=top+131-fraction*112
        svg += [f'<line x1="68" x2="1070" y1="{y}" y2="{y}" stroke="#dbe4eb"/>',
                f'<text x="60" y="{y+4}" text-anchor="end" font-size="13">{maximum*fraction:,.0f}</text>']
    for i,f in enumerate(conductors):
        x=550 if len(conductors)==1 else 88+i*958/(len(conductors)-1)
        for method,label,color in methods:
            m=lookup[(method,f)]['metrics'][key]
            dx=-3 if method.startswith('ec.') else 3
            svg += [f'<line x1="{x+dx}" x2="{x+dx}" y1="{top+131-m["max"]/maximum*112}" y2="{top+131-m["min"]/maximum*112}" stroke="{color}" stroke-width="2"/>',
                    f'<circle cx="{x+dx}" cy="{top+131-m["mean"]/maximum*112}" r="3.5" fill="{color}"/>']
        if pi==2:
            label='?' if f is None else f'{f/1e6:.2f}'
            n=lookup[(methods[0][0],f)]['metrics'][key]['n']
            svg += [f'<text transform="translate({x},{top+151}) rotate(55)" font-size="13">{label}</text>',
                    f'<text x="{x}" y="{top+209}" text-anchor="middle" font-size="12">{n}</text>']
svg += ['<text x="24" y="731" font-size="14">Conductor f in millions, ordinal positions; ? keeps unknown labels. Bottom numbers: target observations per arm.</text>',
        '<text x="24" y="752" font-size="13">Point orders vary and are certified; these samples do not isolate a causal conductor effect or establish a priced total.</text>','</g></svg>']
(out/'operations.svg').write_text('\n'.join(svg)+'\n')
flow='''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="220" viewBox="0 0 1100 220"><rect width="1100" height="220" fill="white"/><defs><marker id="a" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><path d="M0,0 L10,4 L0,8" fill="#386078"/></marker></defs><g fill="#edf3f7" stroke="#386078" stroke-width="2"><rect x="20" y="20" width="220" height="170" rx="10"/><rect x="300" y="20" width="250" height="70" rx="10"/><rect x="610" y="20" width="250" height="70" rx="10"/><rect x="300" y="120" width="560" height="70" rx="10"/><rect x="920" y="20" width="160" height="170" rx="10"/><path d="M240 55 H295 M240 155 H295 M550 55 H605 M860 55 H915 M860 155 H915" marker-end="url(#a)"/></g><g font-family="sans-serif" fill="#17354a" text-anchor="middle" font-size="17"><text x="130" y="74">Public G, Q</text><text x="130" y="103">Exact point order M</text><text x="130" y="132">Same target in both arms</text><text x="425" y="48">Reduced Tate / Miller</text><text x="425" y="74">Primitive image checked</text><text x="735" y="48">FF PH + factor-base relations</text><text x="735" y="74">Small-prime digits + CRT</text><text x="580" y="148">Direct EC PH / bounded BSGS digits + CRT</text><text x="580" y="175">Same known smooth factorization; no unreduced rho baseline</text><text x="1000" y="79">Recovered x</text><text x="1000" y="108">EC point check</text><text x="1000" y="137">Truth check</text></g></svg>'''
(out/'pipeline.svg').write_text(flow)
global_rows=[]
for r in doc['overall']:
    for metric,title in panels:
        m=r['metrics'][metric]
        label=next(v[1] for v in methods if v[0]==r['method'])
        short={'ec_adds_and_doubles':'Curve add + double','explicit_fp_muls_and_squares':'Explicit Fp mul + square','ec_field_inversions_unexpanded':'Curve-field inverse'}[metric]
        global_rows.append(', '.join(f'[{x}]' for x in [label,short,m['n'],f'{m["min"]:,.0f}',f'{m["mean"]:,.2f}',f'{m["max"]:,.0f}'])+',')
replay='replay' in name
population='original source model only' if name.startswith('source') else 'every retained C-prefix model'
contended=str(iso['run']['contended']).lower()
totals=doc['totals'];ff='ff.k1_ph_relations.full_point_v1/'
source=f'''#set page(paper:"a4",margin:13mm)
#set text(size:9pt)
= Original C: complete selected point orders
Run {name}. Started {started}.
Host: {iso['host']['cpu_model']}, {iso['host']['machine']}, kernel {iso['host']['kernel']}.
Population: {population}. {'Separate replay; never pooled with primary observations.' if replay else 'Primary observations; replays are kept separately.'}

#image("pipeline.svg",width:100%)

Native accounting: {doc['records']:,} method records; {doc['paired_verified']:,} matched verified targets.
Point orders M range from {doc['point_order_min']} to {doc['point_order_max']}.
Both arms recover modulo the complete verified M. A selected point need not
generate the full curve group; maximality is not asserted.

The field is p=1637958981148213249, and N=p-1=2^26*3^20*7.
Every selected point order is smooth over 2,3,7. Its large numeric value is
not a large-prime security level. Both arms exploit the known factorization.
Solver arguments contain the public curve, G, Q and M, never scalar truth.
The runner verifies the recovered point and the separate public known answer.

The finite-field route varies the first argument of reduced Tate/Miller,
checks exact prime-power image order, then performs multiplicative PH.
Initially unknown projected small-prime base logs are recovered by checked
one-/two-factor relations and modular Gaussian elimination. No generic-log
fallback is used. This is a smooth-order diagnostic, not NFS/integer-smoothness
index calculus, a large-k backend or a new asymptotic result.

Recorded finite-field setup equations/queries:
{totals.get(ff+'base_relations_relations',0):,} / {totals.get(ff+'base_relations_queries',0):,}.
Target digit equations/queries:
{totals.get(ff+'target_relations_relations',0):,} / {totals.get(ff+'target_relations_queries',0):,}.
Identity and dependent setup equations remain counted; per-component ranks,
matrix work and pair tests are retained. Tiny prime digit groups explain the
dense relation coverage; these yields do not describe a large-order oracle.

Fixture construction, exact order finding and known-target generation are
recorded separately, outside both solver intervals. Cold solver work includes
input-order checks, projections, method setup, auxiliary searches, pairing,
relations/matrix work, digit recovery, CRT and final EC point verification.
Raw end-to-end observations remain in the native CSV, with min/max/mean and n.
No admitted runtime comparison is made: one round, no A/A, PSI unavailable.
CPU 4 was reserved; wrapper observed contended={contended}. Replays do not
contribute to primary means. The full-prefix replay is explicitly contended.

The immutable C prefix contains 1,869 vertices and 4,903 replayed maps;
1,024 vertices are expanded. Local labels certify 517 vertices / 24 observed
strata and keep 1,352 unknown. Original full C remains 1,934,901,248 vertices
and 165 strata. This does not complete original A--E, the full eight-target /
five-round matrix, transported-target comparisons, representation registration
or independent arithmetic replay.

Raw records, fixtures, scalar truth, native CSV/statistics, isolation journals,
exact binaries and both source/replay populations accompany the delivery.
Native implementation is source-pinned at d954e80d9; linked-library source
2bf7700b7. Source hash #text(size:7pt)[#raw("{meta['source_sha256']}")].

#pagebreak()
= Recorded primitive work across local strata
#image("operations.svg",width:100%)
#set text(size:7pt)
#table(columns:(1.2fr,1.8fr,.5fr,.9fr,1fr,.9fr),inset:3pt,
[Arm], [Metric], [n], [Min], [Mean], [Max],
{chr(10).join(global_rows)})
#set text(size:8pt)
These are native descriptive summaries. Curve additions/doublings, curve-field
inversions and explicit Fp multiplication/square counts are different units;
they are not a priced total. EC inversion and square-root internals remain
unexpanded. Explicit Fp inversions already contribute to the multiplication /
square ledger and must not be priced twice. Allocation, comparisons, scalar
integer work and Gaussian inverse details remain unpriced with counters where
available. All per-phase counters and present-value n are in the native CSV.
Observed ranges are not confidence intervals. Unknown strata remain separate.
Source and replay plots may match because every non-clock field reproduces;
they are not extra independent samples or an independent arithmetic proof.
'''
(out/'report.typ').write_text(source)
typst.compile(str(out/'report.typ'),output=str(out/'report.pdf'))
for figure,height in [('operations',760),('pipeline',220)]:
    (out/f'{figure}.typ').write_text(f'#set page(width:1100pt,height:{height}pt,margin:0pt)\n#image("{figure}.svg",width:100%)\n')
    typst.compile(str(out/f'{figure}.typ'),output=str(out/f'{figure}.pdf'))
