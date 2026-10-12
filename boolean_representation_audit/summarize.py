"""Render tables from the saved measurements, without rerunning timings."""
import json
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'runs/20260929/results.json').read_text())
lines=['# Boolean row representation audit: exactness passes, performance is conditional','',
'2026-09-29 UTC. Standalone synthetic algebra, 6 and 8 polynomial variables. These are not field degrees.', '',
'The packed backend preserves every certified row and algebraic-work counter of the recovered sparse backend under each matching schedule. The receipt audit reconciles 504 fresh-process samples, 1,416 independently verified computations, and 48 fresh output recomputations. Six new test methods include 60 random systems and the 12 frozen systems; the original nine tests also pass.', '',
'Classification: representation engineering / stage diagnostic. No novel algebraic algorithm, asymptotic improvement, blocked-workload result, ECC experiment or end-to-end cryptographic improvement is established. Full ideal closure still has an exponential output ceiling. The m=83 confidence gate is unperformed.', '',
'## Decision', '',
'The exhaustive schedule meets the predeclared numerical four-holdout criterion; frontier scheduling fails it because the chain and cycle controls stay near parity. This is not a blanket packed-backend performance claim. Baseline-only A/A controls expose substantial host noise (one tiny-case pair reaches 27.88x), and process/NUMA controls were incomplete. Treat timing improvements as provisional observations on this virtualized host, not validated hardware rates. No post-hoc tuning or reruns selected for favorable timing were performed.', '',
'The useful structural observation is that sparse input does not imply sparse intermediates: the new sparse quadratic system starts with at most six terms per equation and reaches 50 during frontier elimination. Both new random quadratic systems have a complete equation-level interaction graph (min-fill width 7), unlike the width-1 chain and width-2 cycle. These are diagnostic examples, not a statistical relationship or optimal-width proof for arbitrary systems.', '',
'## Paired representation comparison', '',
'Each ratio is packed time / sparse time under the same schedule; lower is less time. Intervals are exploratory 95% percentile bootstrap intervals on seven paired worker means, without multiple-comparison correction. Compute includes setup, conversion, products, queueing, elimination and certificate creation. The separate total adds independent certificate verification; neither includes process startup/imports, JSON or fingerprinting.', '',
'| Case | Schedule | Compute ratio [95% interval] | Compute + verify ratio | A/A min–max | Python allocation ratio |',
'| --- | --- | ---: | ---: | ---: | ---: |']
for c in r['cases']:
 aa=[p['b_over_a'] for p in c['aa']]
 for schedule,s in c['schedules'].items():
  q=s['compute_s_packed_over_sparse']; lo,hi=q['bootstrap_95pct']; total=s['compute_plus_verify_s_packed_over_sparse']['median']; arms=s['arms']
  lines.append(f"| {c['name']} | {schedule} | {q['median']:.3f} [{lo:.3f}, {hi:.3f}] | {total:.3f} | {min(aa):.3f}–{max(aa):.3f} | {arms['packed']['peak_compute_python_bytes']/arms['sparse']['peak_compute_python_bytes']:.3f} |")
lines+=['','A/A controls use the sparse frontier schedule only. Broad A/A ranges block fine-grained timing conclusions. Fresh untraced worker RSS values span 14,848–15,360 KiB and include interpreter/import/verifier memory; they do not establish a process-memory gain. Python allocations are separate traced compute-only runs. Packed payload bytes are minimal bit payload, while sparse payload counts 32-bit indices; neither is actual process memory.', '',
'## Absolute stage times and deterministic work', '',
'Every row below is verified. Counters are algebraic row counts, not calibrated machine operations. The exact output and work counts match between backends of the same schedule; representation cost alone changes.', '',
'| Case | Variant | Median compute ms | Median compute + verify ms | Submitted | XOR rows | Rank | Stored terms | Peak row terms |',
'| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
for c in r['cases']:
 for schedule,s in c['schedules'].items():
  for backend,a in s['arms'].items():
   st=a['stats']; lines.append(f"| {c['name']} | {backend}/{schedule} | {1000*a['compute_s']['median']:.3f} | {1000*a['compute_plus_verify_s']['median']:.3f} | {st['submitted']} | {st['xor_steps']} | {st['rank']} | {st['stored_terms']} | {st['peak_row_terms']} |")
lines+=['','The inherited peak-row counter observes post-XOR rows and stored pivots, not every incoming row; retain the maximum input-row size separately. It is sufficient to witness the reported 6-to-50 growth but must not be described as a universal all-transient peak.', '',
'## Host, provenance and deviations', '',
'- CPU: Intel Xeon Platinum 8370C, virtualized Linux x86-64, Python 3.12. CPU 0 affinity succeeded in all samples; one virtual NUMA node exposed. Affinity is not exclusive core reservation.',
'- Node-0 memory binding returned Invalid argument. Physical NUMA placement and DDR generation are unknown. The runner requested a process summary with ps -eo; that command failed with fatal library error, lookup self, leaving an empty list. No process-isolation claim is supported. Load averages are saved before/after.',
'- Five A/A pairs per case precede seven alternating AB/BA rounds per schedule. Each untraced worker performs three fresh computations. Four separate traced workers per case measure Python allocations. Raw failures/timeouts: none. Host noise and memory-binding failure are retained, not hidden.',
'- Source, inputs and protocol were committed before measurement: aa56cbb51ff53d8df64d65bdebd8aaefeb594d51. PRE_RUN_SHA256SUMS verifies every frozen source/input and the preserved previous artifact. The source run has not been retuned.',
'- A post-run receipt-auditor import-order error was fixed before its successful execution; measured sources and inputs were unchanged. This was audit plumbing, not a rerun or a numerical-data correction.',
'- The prior artifact remains unchanged in prior/. This follow-up adds a backend by loading the identical closure engine with different row/space/reducer classes. The packed backend is a bounded reference, not a scalable solver.', '',
'## Reproduce and audit', '',
'Run from boolean_representation_audit with Python 3.12 (standard library only):', '',
'```sh',
'sha256sum -c PRE_RUN_SHA256SUMS',
'python3 -m unittest -v test_backends',
'python3 audit_results.py',
'python3 measure.py --output runs/my-new-run',
'```', '',
'The last command refuses an existing directory. The delivered result is runs/20260929/results.json. The 504 individual receipts are compressed losslessly as samples.jsonl.gz; each JSON line has filename and data fields. Inspect with python3 -m gzip -d on a copy, or gzip.open in Python. SAMPLE_MANIFEST.json records bytes and SHA-256. Individual sample paths in results.json refer to the filename fields in that archive. AUDIT_RECEIPT.json and TEST_RECEIPT.txt retain verification output. The CI job replays correctness and audits evidence; it does not impose timing thresholds on shared runners.', '',
'## Next bounded question', '',
'- [x] Preserve and replay the original negative control.',
'- [x] Compare exact sparse and packed representations with identical scheduling and outputs.',
'- [x] Record intermediate fill-in, interaction graph, allocations, process RSS and raw timings.',
'- [x] Preserve failed performance gates and host-control limitations.',
'- [ ] Compare a compact support-sharing representation with these two references on a newly frozen, small synthetic suite. Include low-width and fully coupled inputs and charge construction/conversion. Do not raise the current variable limit or integrate an application.',
'- [ ] Test preservation of narrow variable interactions before considering a chordal method. A min-fill diagnostic is not chordal elimination or a completeness proof.',
'- [ ] Obtain a quieter, verifiably controlled host before treating the observed timing reductions as hardware results.', '',
'Prior art remains relevant: [PolyBoRi](https://polybori.sourceforge.net/features.html) uses shared ZDD representations for Boolean polynomials; [Cifuentes–Parrilo](https://arxiv.org/abs/1604.02618) develops chordal networks for structured polynomial ideals. Their features/abstract were reviewed in this follow-up. These experiments do not implement either system or establish novelty.', '']
(ROOT/'README.md').write_text('\n'.join(lines))
summary={'scope':'standalone synthetic algebra stage diagnostic','source_commit':'aa56cbb51ff53d8df64d65bdebd8aaefeb594d51','cases':12,'worker_samples':504,'verified_computations':1416,'gate':{},'end_to_end_speedup':None,'rows':[]}
for strategy in ('exhaustive','frontier'):
 summary['gate'][strategy]=all(c['schedules'][strategy]['compute_s_packed_over_sparse']['median']<.9 and c['schedules'][strategy]['compute_s_packed_over_sparse']['bootstrap_95pct'][1]<1 for c in r['cases'] if c['role']=='holdout')
for c in r['cases']:
 if c['role']=='holdout':
  s=c['schedules']['frontier']; summary['rows'].append({'case':c['name'],'width_upper_bound':c['structure']['min_fill_width_upper_bound'],'max_input_terms':c['structure']['max_input_row_terms'],'peak_reduced_terms':s['arms']['sparse']['stats']['peak_row_terms'],'ratio':s['compute_s_packed_over_sparse']})
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
