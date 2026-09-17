#!/usr/bin/env python3
from runbench import ROOT, bench
arms=[("reset",ROOT/"build/g7-xonly-sparse-bridge13-pipeline-reset",524288,16),("pipeline",ROOT/"build/g7-xonly-sparse-bridge13-pipeline",524288,16)]
for rep in range(3):
    order=arms if rep%2==0 else list(reversed(arms))
    for label,binary,workers,batch in order:
        bench(binary,f"reset-r{rep}-{label}",workers=workers,batch=batch,run_id=62573)
