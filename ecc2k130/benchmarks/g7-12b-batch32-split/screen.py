#!/usr/bin/env python3
import json,math,pathlib,statistics
import runbench
ROOT=pathlib.Path(__file__).resolve().parents[2];OUT=pathlib.Path(__file__).resolve().parent
variants={
 'selected':(ROOT/'build/ecc2k130-local-packed',524286,16),
 'batch24':(ROOT/'build/g7-split-batch24',349524,24),
 'batch32':(ROOT/'build/g7-split-batch32',262143,32),
}
orders=[['selected','batch24','batch32'],['batch32','batch24','selected'],['batch24','selected','batch32']]
rows=[]
for rep,order in enumerate(orders):
 for name in order:
  binary,workers,batch=variants[name]
  rows.append(runbench.bench(binary,f'timing-r{rep}-{name}',workers=workers,batch=batch,run_id=58401))
rates={name:[] for name in variants}
for row in rows:rates[row['label'].split('-',2)[2]].append(row['million_updates_per_second']/1000)
base=rates['selected'];tcrit=4.3026527297;summary=[]
for name,values in rates.items():
 paired=[values[i]/base[i] for i in range(3)];logs=[math.log(x) for x in paired];mean=statistics.mean(logs);se=statistics.stdev(logs)/math.sqrt(3) if name!='selected' else 0;lo=math.exp(mean-tcrit*se);hi=math.exp(mean+tcrit*se)
 summary.append({'variant':name,'median_billion_updates_per_second':statistics.median(values),'rate_over_12b':statistics.median(values)/12,'paired_geomean_speedup':math.exp(mean),'paired_speedup_95_ci':[lo,hi],'rates':values,'decision':'reference' if name=='selected' else ('screen qualifier' if lo>1 else 'unconfirmed or regression')})
out={'passed':all(row['passed'] for row in rows),'unit':'billion_complete_scalar_updates_per_second','scalar_walks_per_arm':8388576,'updates_per_sample':34359607296,'timing_samples':len(rows),'rows':summary,'raw_labels':[row['label'] for row in rows],'goal_billion_updates_per_second':12,'goal_met':max(statistics.median(v) for v in rates.values())>=12,'generic_work_boundary':'sqrt(n/262)','generic_work_ratio':1,'full_dlp_S':None}
(OUT/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
