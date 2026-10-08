#!/usr/bin/env python3
from runbench import bench,ROOT
selected=ROOT/'build/ecc2k130-local-packed';candidate=ROOT/'build/g7-xonly-doubling'
for rep in range(3):
 order=[('selected',selected),('doubling-only',candidate)]
 if rep&1: order.reverse()
 for label,binary in order: bench(binary,f'diagnostic-r{rep}-{label}',run_id=62473)
