#!/usr/bin/env python3
from runbench import ROOT, bench
binary=ROOT/"build/g7-xonly-sparse-bridge13-pipeline-mod72"
for rep in range(3):
    bench(binary,f"mod72-confirm-r{rep}",workers=524288,steps=1024,launches=8,batch=16,run_id=62574)
