#!/usr/bin/env python3
import json,math,statistics
from pathlib import Path
OUT=Path(__file__).resolve().parent;TCRIT=4.302652729911275;labels=("bridge13","pipeline")
rates={label:[json.loads((OUT/f"pipeline-r{rep}-{label}.json").read_text())["million_updates_per_second"]/1000 for rep in range(3)] for label in labels};base=rates["bridge13"];rows=[]
for label in labels:
 logs=[math.log(rates[label][i]/base[i]) for i in range(3)];mean=statistics.mean(logs);error=0 if label=="bridge13" else TCRIT*statistics.stdev(logs)/math.sqrt(3)
 rows.append({"variant":label,"rates_billion_per_second":rates[label],"median_billion_per_second":statistics.median(rates[label]),"paired_geomean_vs_bridge13":math.exp(mean),"paired_95_ci":[math.exp(mean-error),math.exp(mean+error)]})
result={"status":"complete-map-pipelined-classification-screen","checkpoint_mismatches":{"pipeline":0},"collision_work_ratio_vs_selected":171.018/169.8805,"rows":rows,"decision":"retain pipelined bridge-1/bridge-3 as the fastest verified complete map; continue toward 12 B/s"}
(OUT/'bridge13-pipeline-comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
