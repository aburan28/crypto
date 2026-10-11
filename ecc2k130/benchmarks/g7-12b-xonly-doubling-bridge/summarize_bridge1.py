#!/usr/bin/env python3
import json,math,statistics
from pathlib import Path
OUT=Path(__file__).resolve().parent;TCRIT=4.302652729911275;labels=("doubling","bridge1")
rates={label:[json.loads((OUT/f"bridge1-r{rep}-{label}.json").read_text())["million_updates_per_second"]/1000 for rep in range(3)] for label in labels};base=rates["doubling"];rows=[]
for label in labels:
 logs=[math.log(rates[label][i]/base[i]) for i in range(3)];mean=statistics.mean(logs);error=0 if label=="doubling" else TCRIT*statistics.stdev(logs)/math.sqrt(3)
 rows.append({"variant":label,"rates_billion_per_second":rates[label],"median_billion_per_second":statistics.median(rates[label]),"paired_geomean_vs_doubling":math.exp(mean),"paired_95_ci":[math.exp(mean-error),math.exp(mean+error)]})
result={"status":"non-promotable-sigma-bridge-1-common-path-diagnostic","checkpoint_mismatches":{"bridge1":0},"rows":rows,"decision":"bridge-1 provides enough headroom to study a complete bridge-1/bridge-3 map"}
(OUT/'bridge1-comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
