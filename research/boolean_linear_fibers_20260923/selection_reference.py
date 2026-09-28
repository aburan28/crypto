#!/usr/bin/env python3
"""Exact bounded discovery census of variable sets that leave linear fibers."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def interaction_graph(polys,n):
    assert 0<=n<=24
    graph=[0]*n
    for p in polys:
        parity=set()
        for m in p:
            assert 0<=m<1<<n and m.bit_count()<=2
            if m in parity:parity.remove(m)
            else:parity.add(m)
        for m in parity:
            if m.bit_count()==2:
                i=(m&-m).bit_length()-1;j=(m^(1<<i)).bit_length()-1
                graph[i]|=1<<j;graph[j]|=1<<i
    return graph
def independent(graph,mask):
    return all(not (mask&neighbors) for i,neighbors in enumerate(graph) if mask&(1<<i))
def choose(a,b):
    return a if (a.bit_count(),-a)>=(b.bit_count(),-b) else b
def maximum_independent(graph):
    n=len(graph);assert n<=24
    left=n//2;right=n-left
    best=[0]*(1<<right)
    for mask in range(1,1<<right):
        bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
        best[mask]=choose(best[rest],bit|best[rest&~(graph[left+i]>>left)])
    answer=0;left_count=0
    for mask in range(1<<left):
        if not independent(graph,mask):continue
        left_count+=1;allowed=(1<<right)-1
        for i in range(left):
            if mask&(1<<i):allowed &= ~(graph[i]>>left)
        answer=choose(answer,mask|(best[allowed]<<left))
    assert independent(graph,answer)
    return answer,{'right_subset_states':len(best),'left_subsets':1<<left,'independent_left_subsets':left_count}
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    out=parser.parse_args().out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=HERE/'resource_probe_02'
    for name,digest in json.loads((source/'manifest.json').read_text())['files'].items():assert sha(source/name)==digest
    shutil.copyfile(Path(__file__),out/'independent_set.py')
    rows=[]
    for n in [16,24]:
        for seed in [17,937]:
            for family in ['planted','cross_planted','unplanted']:
                path=source/f'n{n}-{seed}-{family}.jsonl'
                with path.open() as f:fixture=json.loads(next(f))
                graph=interaction_graph(fixture['polys'],n);mask,work=maximum_independent(graph)
                rows.append({'n':n,'seed':seed,'family':family,'input_sha256':sha(path),'polys':fixture['polys'],
                    'graph':graph,'edge_count':sum(g.bit_count() for g in graph)//2,'independent_mask':mask,
                    'independent_variables':[i for i in range(n) if mask&(1<<i)],'maximum_size':mask.bit_count(),
                    'outside_assignments':1<<(n-mask.bit_count()),'assignments_per_linear_fiber':1<<mask.bit_count(),
                    'exact_search_work':work,'complete':True,'solver_cost':None,'speed_ratio':None})
    value={'schema_version':1,'scope':'Discovery structure only; no new solver or timing claim. All costs and speed ratios are unmeasured.',
        'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_manifest_sha256':sha(source/'manifest.json'),'cells':rows}
    (out/'results.json').write_text(json.dumps(value,indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps({'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}},indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['n','seed','family','maximum_size','edge_count']} for r in rows]))

if __name__=='__main__':main()
