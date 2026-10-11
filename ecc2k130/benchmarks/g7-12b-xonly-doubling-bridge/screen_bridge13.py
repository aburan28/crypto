#!/usr/bin/env python3
from runbench import ROOT, bench
arms=[("bridge13",ROOT/"build/g7-xonly-sparse-bridge13",524288,16),("block",ROOT/"build/g7-xonly-sparse-bridge-block",524288,16)]
for rep in range(3):
    order=arms if rep%2==0 else list(reversed(arms))
    for label,binary,workers,batch in order:
        bench(binary,f"bridge13-r{rep}-{label}",workers=workers,batch=batch,run_id=62573)
