#!/usr/bin/env python3
"""Render frozen measurements into prose and the canonical static scoreboard."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
first=json.loads((HERE/'run-01/comparison.json').read_text())
confirm=json.loads((HERE/'confirmation-01/comparison.json').read_text())
assert first['complete'] and confirm['complete']
variants=json.loads((HERE/'contract.json').read_text())['variants']
lookup={(r['variant'],r['kind'],r['metric']):r for r in first['aggregates']}
clook={(r['variant'],r['metric']):r for r in confirm['aggregates']}
names={'legacy':'Legacy accounting','reference':'Corrected reference','setup':'Borrowed setup','fixed':'Fixed slots','horner':'Grouped Horner','combined':'Combined (retained)'}
def cell(r):return f"{r['ratio']:.6f} [{r['paired_95'][0]:.6f}, {r['paired_95'][1]:.6f}]"
columns=[('micro','setup_s'),('micro','specialization_s'),('stage','process_wall_s'),('full_dlp','process_wall_s')]
lines=['| Variant | Symbolic setup | Specialization | Residual corpus | Cold full-DLP corpus |',
       '|---|---:|---:|---:|---:|']
for v in variants:lines.append('| '+names[v]+' | '+' | '.join(cell(lookup[v,k,m]) for k,m in columns)+' |')
confirmation=['| Variant | Whole-process wall time | Child user+system CPU time |', '|---|---:|---:|']
for v in ['legacy','reference','horner','combined']:
 confirmation.append('| '+names[v]+' | '+cell(clook[v,'wall_s'])+' | '+cell(clook[v,'cpu_s'])+' |')
raw=[json.loads(x) for x in (HERE/'confirmation-01/raw.jsonl').read_text().splitlines()]
checks=sum(json.loads((HERE/'confirmation-01'/x['result']).read_text())[0]['gaudry']['cross_checked'] for x in raw)
micro_ref=json.loads((HERE/'run-01/micro-p67-seed11-r0-reference.json').read_text())
micro_new=json.loads((HERE/'run-01/micro-p67-seed11-r0-combined.json').read_text())
den=micro_ref['targets']*micro_ref['evaluation_repetitions']
bm=micro_ref['specialization_fp_muls']//den;cm=micro_new['specialization_fp_muls']//den
report=f'''# Results: each summation-polynomial path pursued

**Retained:** borrowed symbolic expansions plus once-per-curve fixed-slot grouping
and Horner specialization, with corrected operation accounting. The independent
fixed-accumulation candidate is validated and preserved, but its slower per-target
path is superseded by grouped Horner. The prior row/cache allocation optimization
remains in place.

The three local results relative to the corrected reference are:

* Borrowing cached expansions lowers symbolic setup time by about 1.0%.
* Fixed-slot accumulation lowers specialization time by about 49.4% without
  changing the field-operation count.
* Grouped Horner lowers specialization time by about 55.4%, and lowers actual
  specialization multiplications. The combined setup/Horner candidate retains
  those evaluation gains and lowers setup time by about 1.2%.

These are **local CPU measurements**, not a new end-to-end ECDLP speedup.
The independent full-DLP confirmation below establishes neither an improvement
nor a regression. Its CPU-time control explains why the apparent initial
full-DLP gains should not be promoted to an overall speedup claim.

## Original frozen six-way ablation

One unit in every numerical cell: **candidate time / corrected-reference time**,
with input-cluster bootstrap 95% intervals in brackets. All variants perform the
same verified workload within a column. Micro timings separate symbolic setup
from specialization; residual and full-DLP columns include entire child processes.
The reference is explicitly 1.0. Absolute times and every pair remain in
[the comparison](run-01/comparison.json).

{chr(10).join(lines)}

The original 342 runs completed with no failures. This comprises 144 micro runs
(589,824 timed specializations), 108 stage runs (6,480 verified residuals), and
90 verified cold full-DLP runs (16,200 independent oracle cross-checks).
Every specialized-polynomial digest agrees, as do complete decomposition sets,
relation statistics, rho results and all counters apart from exactly reconciled
specialization charges. All 21 prior completed stage/full-DLP outputs replay
exactly in the legacy control, ignoring only timing fields. No outliers are
removed. See [raw records](run-01/raw.jsonl) and [provenance](run-01/provenance.json).

## Why the initial full-DLP timing is not the verdict

The accounting-only legacy control also appeared faster than the corrected
reference, even though it executes the same polynomial arithmetic. Individual
paired runs showed descheduling outliers: for example, seed 11 repetition 0 took
0.6083 seconds for reference and 0.4395 for combined, whereas the other two pairs
were about 0.427–0.429 seconds each. This is an unresolved timing confound, not
proof that algebra became substantially cheaper.

Before selecting a default, [a second contract](confirmation-contract.json)
froze five repetitions on the five previous curves plus two fresh curves,
seeds 907 and 1009. It reran legacy, reference, Horner and combined and retained
whole-process wall time **and child CPU time separately**. All 140 runs recovered
the correct scalar and passed {checks:,} independent oracle checks; exact outputs
and specialization-cost deltas agree. The original evidence is retained unchanged.

One unit below: time / corrected-reference time, with input-cluster 95% intervals.

{chr(10).join(confirmation)}

All non-reference confirmation intervals include 1. The retained combined
candidate's CPU-time ratio is {clook['combined','cpu_s']['ratio']:.6f}; this supports
neither a full-DLP speedup nor a regression. The large, repeatable local
specialization benefit is retained, with no claim that it materially accelerates
complete scalar recovery. See [confirmation raw records](confirmation-01/raw.jsonl),
[comparison](confirmation-01/comparison.json), and [provenance](confirmation-01/provenance.json).

## Accounting and limits

At p=67, seed 11, specialization actually costs {bm} Fp multiplications per
reference call and {cm} with Horner: {(1-cm/bm)*100:.2f}% fewer products in this
phase. Symbolic setup uses 41,360 products in both variants. The historical
estimate was 1,770 per specialization, which overcharged the reference by 428.
That correction is separate from the additional 143 products saved by Horner.
The runner checks the exact per-curve delta using the measured microbenchmark
charges and number of solver calls, for every stage and full-DLP report.

This small phase's reduction is diluted by Macaulay elimination, normal forms,
root extraction and relation collection. Full common-operation totals, normalized
S, cost/rho and cost/floor remain unmeasured. Existing reported solver totals are
partial metrics; they include the corrected measured specialization charge, but
do not newly calibrate cold setup and runtime bookkeeping. Classification is
**engineering**, plus an explicitly separated **accounting correction**. The
summation-polynomial degree, relation-yield ceiling and asymptotic boundary do
not change. No extrapolated crossover or GPU/RDMA improvement is claimed.

The validation corpus is small and the host is shared. Confidence intervals are
clustered by curve/seed, with repetitions kept together. The confirmation has
seven input clusters, not 35 independent curves. Full-DLP times include independent
MITM and rho auditing. These limits preclude a general performance guarantee.

## Implementation and reproducibility

All 15 Gaudry unit tests pass for the retained implementation, including direct
termwise specialization checks over every element of F_(7³), sampled larger-field
inputs, zero/base-field targets, symbolic S4 identities, border retries, independent
MITM equivalence and dense/sparse/large-prime scalar recovery. [Test log](gaudry-tests.log).
All source variants and their hashes are preserved under [sources/](sources/).
[selection.json](selection.json) records the default and its source hash.

The current shared checkout also passes `cargo check --lib`. Build commands,
workload definitions and timing boundaries are in [README.md](README.md).
`reference.patch` records the accounting change; the remaining patches compare
each optimization to that corrected reference. The standalone fixed-slot path
remains reproducible even though grouped Horner is the selected implementation.
'''
(HERE/'RESULTS.md').write_text(report)
url='https://github.com/aburan28/crypto/blob/main/research/summation_polynomials_20260914/'
html='''  <div class="panel" id="summation-polynomials-20260914">
    <div class="panel-head"><h2>Summation-polynomial setup and specialization</h2>
      <p>Class <span class="chip">engineering</span>, with a separate accounting correction.
        Six ablations preserve exact polynomials and solver outputs. Retained: borrowed
        symbolic expansions and grouped Horner. Fixed-slot specialization remains a
        validated experimental alternative. About 55% less specialization time;
        no established full-DLP speedup or regression.</p>
      <p><a href="performance-gains.html">Explore interactive performance charts</a> ·
        <a href="performance-gains/summary.pdf">Download chart PDF</a></p></div>
    <div class="table-scroll"><table>
      <caption>Dimensionless time / corrected reference; input-cluster 95% intervals.
        <a href="'''+url+'''run-01/comparison.json">Frozen ablation</a>;
        <a href="'''+url+'''RESULTS.md">results and scope</a>.</caption>
      <thead><tr><th>Variant</th><th>Symbolic setup</th><th>Specialization</th><th>Residual corpus</th><th>Cold full DLP</th></tr></thead><tbody>
'''
for v in variants:
 html+='        <tr><td>'+names[v]+'</td>'+''.join('<td>'+cell(lookup[v,k,m])+'</td>' for k,m in columns)+'</tr>\n'
html+='''      </tbody></table></div>
    <div class="panel-head"><p>The initial full-DLP timing is confounded by shared-host
      outliers in the accounting-only control. A separate 140-run confirmation adds
      two fresh curves and child CPU time. All scalar recoveries and exact cost deltas
      pass; its intervals include 1. Original measurements remain preserved.</p></div>
    <div class="table-scroll"><table><caption>Confirmation time / corrected reference,
      input-cluster 95% intervals. <a href="'''+url+'''confirmation-01/comparison.json">Frozen confirmation</a>.</caption>
      <thead><tr><th>Variant</th><th>Whole-process wall time</th><th>Child CPU time</th></tr></thead><tbody>
'''
for v in ['legacy','reference','horner','combined']:
 html+='        <tr><td>'+names[v]+'</td><td>'+cell(clook[v,'wall_s'])+'</td><td>'+cell(clook[v,'cpu_s'])+'</td></tr>\n'
html+='''      </tbody></table></div>
    <div class="panel-head"><p>On the p=67 seed-11 example, specialization uses 1,342
      counted Fp products in the corrected reference and 1,199 with Horner. The
      historical 1,770 estimate is an accounting issue, not an optimization saving.
      Full common-operation S, cost/rho and cost/floor remain unmeasured. Polynomial
      degree, relation yield, asymptotic exponent and generic-group boundary are unchanged.</p></div>
  </div>

'''
page=REPO/'docs/index-calculus-scoreboard.html';content=page.read_text();anchor='  <div class="panel" id="gaudry-allocation-20260914">'
assert content.count(anchor)==1
marker='  <div class="panel" id="summation-polynomials-20260914">'
if marker in content:
    start=content.index(marker)
    end=content.index(anchor,start)
    content=content[:start]+content[end:]
page.write_text(content.replace(anchor,html+anchor))
print('Rendered frozen results and appended the canonical scoreboard panel')
