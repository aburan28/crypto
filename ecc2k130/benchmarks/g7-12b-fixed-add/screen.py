#!/usr/bin/env python3
from runbench import bench,ROOT
selected=ROOT/'build/ecc2k130-local-packed';candidate=ROOT/'build/g7-fixed-add'
rows=[]
for rep in range(3):
 order=[('selected',selected),('fixed-add',candidate)]
 if rep&1: order.reverse()
 for label,binary in order: rows.append(bench(binary,f'screen-r{rep}-{label}',run_id=62173))
