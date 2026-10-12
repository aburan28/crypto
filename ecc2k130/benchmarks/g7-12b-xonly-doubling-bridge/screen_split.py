#!/usr/bin/env python3
from runbench import bench,ROOT
arms=[('base',ROOT/'build/g7-xonly-doubling',524288,16),('split',ROOT/'build/g7-xonly-sparse-bridge-split',524288,16)]
for rep in range(3):
 order=arms if rep%2==0 else list(reversed(arms))
 for label,binary,workers,batch in order:bench(binary,f'split-r{rep}-{label}',workers=workers,batch=batch,run_id=62573)
