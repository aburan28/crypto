#!/usr/bin/env python3
from runbench import bench,ROOT
arms=[('b16',ROOT/'build/g7-xonly-doubling',524288,16),('b24',ROOT/'build/g7-xonly-doubling-b24',349525,24),('b32',ROOT/'build/g7-xonly-doubling-b32',262144,32)]
for rep in range(3):
 order=arms if rep%2==0 else list(reversed(arms))
 for label,binary,workers,batch in order:bench(binary,f'geometry-r{rep}-{label}',workers=workers,batch=batch,run_id=62573)
