#!/usr/bin/env python3
import json, math, statistics
from pathlib import Path
OUT=Path(__file__).resolve().parent; TCRIT=4.302652729911275
labels=("block","owner")
rates={label:[json.loads((OUT/f"owner-r{rep}-{label}.json").read_text())["million_updates_per_second"]/1000 for rep in range(3)] for label in labels}
base=rates["block"]; rows=[]
for label in labels:
 logs=[math.log(rates[label][i]/base[i]) for i in range(3)]; mean=statistics.mean(logs); error=0 if label=="block" else TCRIT*statistics.stdev(logs)/math.sqrt(3)
 rows.append({"variant":label,"rates_billion_per_second":rates[label],"median_billion_per_second":statistics.median(rates[label]),"paired_geomean_vs_block":math.exp(mean),"paired_95_ci":[math.exp(mean-error),math.exp(mean+error)]})
result={"status":"complete-map-owner-local-rare-phase","checkpoint_mismatches":{"owner":0},"rows":rows,"decision":"reject owner-local rare phase; retain block-local event queue"}
(OUT/'owner-comparison.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
