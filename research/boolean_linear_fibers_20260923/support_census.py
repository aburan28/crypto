#!/usr/bin/env python3
"""Bounded discovery census of equations unaffected by selected variables."""
import argparse
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('selection_reference',HERE/'selection_reference.py')
selection=importlib.util.module_from_spec(spec);spec.loader.exec_module(selection)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def census(polys,n,cap=100000):
    graph=selection.interaction_graph(polys,n);dependencies=[0]*n
    for e,p in enumerate(polys):
        parity=set()
        for m in p:
            if m in parity:parity.remove(m)
            else:parity.add(m)
        support=0
        for m in parity:support|=m
        for j in range(n):
            if support&(1<<j):dependencies[j]|=1<<e
    all_equations=(1<<len(polys))-1;best={};visited=0
    def visit(mask,available,touched):
        nonlocal visited
        if visited==cap:return False
        visited+=1;k=mask.bit_count();untouched=all_equations&~touched
        value={'selected_mask':mask,'untouched_mask':untouched,'untouched_equations':untouched.bit_count()}
        if k not in best or (value['untouched_equations'],-mask)>(best[k]['untouched_equations'],-best[k]['selected_mask']):best[k]=value
        while available:
            bit=available&-available;available^=bit;j=bit.bit_length()-1
            if not visit(mask|bit,available&~graph[j],touched|dependencies[j]):return False
        return True
    complete=visit(0,(1<<n)-1,0)
    return {'complete':complete,'visited_sets':visited,'cap':cap,'best_by_size':best,'dependencies':dependencies}
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    out=parser.parse_args().out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=HERE/'resource_probe_02'
    for name,digest in json.loads((source/'manifest.json').read_text())['files'].items():assert sha(source/name)==digest
    for name in ['support_census.py','selection_reference.py']:shutil.copyfile(HERE/name,out/name)
    cells=[]
    for n in [16,24]:
        for seed in [17,937]:
            for family in ['planted','cross_planted','unplanted']:
                path=source/f'n{n}-{seed}-{family}.jsonl'
                with path.open() as f:fixture=json.loads(next(f))
                got=census(fixture['polys'],n)
                graph=selection.interaction_graph(fixture['polys'],n);mask,_=selection.maximum_independent(graph)
                touched=0
                for i in range(n):
                    if mask&(1<<i):touched|=got['dependencies'][i]
                cells.append({'n':n,'seed':seed,'family':family,'polys':fixture['polys'],'source_sha256':sha(path),
                    'current_selected_mask':mask,'current_untouched_equations':len(fixture['polys'])-touched.bit_count(),
                    **got,'solver_cost':None,'speed_ratio':None})
    result={'scope':'Discovery structure only. No shortcut solver or timing claim. Capped enumeration does not certify an optimum.',
        'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_manifest_sha256':sha(source/'manifest.json'),'cells':cells}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps({'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}},indent=2)+'\n')
    print(json.dumps([{'n':c['n'],'seed':c['seed'],'family':c['family'],'current_untouched':c['current_untouched_equations'],
        'best_untouched_by_size':{k:v['untouched_equations'] for k,v in c['best_by_size'].items()},'complete':c['complete'],'visited':c['visited_sets']} for c in cells]))

if __name__=='__main__':main()
