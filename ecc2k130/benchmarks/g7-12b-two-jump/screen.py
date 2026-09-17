#!/usr/bin/env python3
import json,math,pathlib,statistics
import runbench
ROOT=pathlib.Path(__file__).resolve().parents[2];OUT=pathlib.Path(__file__).resolve().parent
collision=json.loads((OUT/'collision-summary.json').read_text())
ratio=collision['two_jump_over_eight_jump_mean_iterations']
variants={'selected':ROOT/'build/ecc2k130-local-packed','two-jump-indexed':ROOT/'build/g7-two-jump-indexed','two-jump-immediate':ROOT/'build/g7-two-jump-immediate'}
orders=[['selected','two-jump-indexed','two-jump-immediate'],['two-jump-immediate','two-jump-indexed','selected'],['two-jump-indexed','selected','two-jump-immediate']]
rows=[]
for rep,order in enumerate(orders):
 for name in order: rows.append(runbench.bench(variants[name],f'timing-r{rep}-{name}',run_id=49271))
rates={n:[] for n in variants}
for row in rows: rates[row['label'].split('-',2)[2]].append(row['million_updates_per_second']/1000)
summary=[];base=rates['selected'];tcrit=4.3026527297
for name,v in rates.items():
 ratios=[v[i]/base[i] for i in range(3)];logs=[math.log(x) for x in ratios];mean=statistics.mean(logs);se=statistics.stdev(logs)/math.sqrt(3) if name!='selected' else 0
 work_ratio=1 if name=='selected' else ratio
 summary.append({'variant':name,'median_billion_updates_per_second':statistics.median(v),'rate_over_12b':statistics.median(v)/12,'collision_work_ratio':work_ratio,'collision_adjusted_billion_iterations_per_second':statistics.median(v)/work_ratio,'paired_geomean_speedup':math.exp(mean),'paired_speedup_95_ci':[math.exp(mean-tcrit*se),math.exp(mean+tcrit*se)],'rates':v,'decision':'reference' if name=='selected' else ('screen qualifier' if math.exp(mean-tcrit*se)>1 else 'unconfirmed or regression')})
out={'passed':all(r['passed'] for r in rows),'unit':'billion_complete_scalar_updates_per_second','updates_per_sample':34359738368,'timing_samples':len(rows),'rows':summary,'raw_labels':[r['label'] for r in rows],'goal_billion_updates_per_second':12,'goal_met':max(statistics.median(v) for v in rates.values())>=12,'generic_work_boundary':'sqrt(n/262)','generic_work_ratio':{'selected':1,'two_jump':ratio},'full_dlp_S':None}
(OUT/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
