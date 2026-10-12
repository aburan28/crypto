import datetime,hashlib,json,platform,sys,time
from pathlib import Path
sys.path.insert(0,'/opt/homebrew/opt/cryptominisat/lib/python3.14/site-packages')
import pycryptosat
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'resource_probe_04'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
metadata={'scope':'Discovery-only native SAT feasibility probe. Python wrapper, not a same-binary performance comparison; no promotion claim. Solver context, encoding, solve, model checking and destruction are timed; interpreter/import startup is excluded.','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'executable':sys.executable,'machine':platform.machine(),'binding_path':pycryptosat.__file__,'binding_sha256':sha(Path(pycryptosat.__file__)),'library_sha256':sha(Path('/opt/homebrew/opt/cryptominisat/lib/libcryptominisat5.5.14.dylib')),'threads':1,'time_limit_seconds':5.0,'conflict_limit':200000,'source_fixtures':{f:sha(SOURCE/(f+'.jsonl')) for f in ['planted','cross_planted','unplanted']},'complete':False}
(HERE/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
records=[]
for family in ['planted','cross_planted','unplanted']:
    fixture=json.loads((SOURCE/(family+'.jsonl')).read_text().splitlines()[0]);polys=fixture['polys'];n=fixture['n']
    start=time.perf_counter_ns();solver=pycryptosat.Solver(threads=1,verbose=0,time_limit=5.0,confl_limit=200000)
    nonlinear=sorted({m for p in polys for m in p if m.bit_count()>=2})
    auxiliary={m:n+1+i for i,m in enumerate(nonlinear)};clauses=0
    for m,z in auxiliary.items():
        variables=[i+1 for i in range(n) if m&(1<<i)]
        for x in variables:solver.add_clause([x,-z]);clauses+=1
        solver.add_clause([-x for x in variables]+[z]);clauses+=1
    for p in polys:
        variables=[auxiliary[m] if m.bit_count()>=2 else m.bit_length() for m in p if m]
        solver.add_xor_clause(variables,0 in p)
    encoded=time.perf_counter_ns();status,assignment=solver.solve();solved=time.perf_counter_ns()
    model=None if status is not True else sum(1<<i for i in range(n) if i+1<len(assignment) and assignment[i+1] is True)
    verified=False
    if status is True:
        assert all(sum((model&m)==m for m in p)%2==0 for p in polys);verified=True
    elif status is False:
        assert fixture['search_reference']=='UNSAT';verified=True
    del solver
    total=time.perf_counter_ns()-start
    row={'n':n,'seed':fixture['seed'],'family':family,'status':'SAT' if status is True else 'UNSAT' if status is False else 'UNKNOWN','model':model,'verified':verified,'variables':n+len(nonlinear),'cnf_clauses':clauses,'xor_clauses':len(polys),'encoding_ns':encoded-start,'solver_ns':solved-encoded,'observed_total_ns':total,'completion_ns':total if verified else None}
    records.append(row);print(json.dumps(row),flush=True)
(HERE/'results.json').write_text(json.dumps(records,indent=2)+'\n')
metadata.update(complete=True,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(HERE/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
(HERE/'manifest.json').write_text(json.dumps({'files':{p.name:sha(p) for p in sorted(HERE.iterdir()) if p.is_file()},'scope':metadata['scope']},indent=2)+'\n')
