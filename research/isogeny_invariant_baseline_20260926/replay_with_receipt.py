"""Run, independently verify, and retain execution evidence, including failures."""
import argparse,copy,datetime,hashlib,json,os,platform,subprocess,sys,time,traceback
from pathlib import Path
from verify_independent import verify

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def need(ok,message):
    if not ok:raise RuntimeError(message)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--require-clean',action='store_true');args=ap.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    receipt={'schema_version':1,'status':'running','started_utc':utc(),'commands':[],
             'evidence_class':'executed_and_independently_recomputed_finite_case',
             'formal_proof':False,'source_files':{},'input_files':{},'artifacts':{}}
    try:
        need(sys.flags.optimize==0,'Python optimization must be disabled: benchmark uses assertions')
        git=lambda *a:subprocess.check_output(['git','-C',str(ROOT),*a],text=True).strip()
        receipt['source_commit']=git('rev-parse','HEAD')
        receipt['worktree_dirty']=bool(git('status','--porcelain','--untracked-files=all'))
        if args.require_clean:need(not receipt['worktree_dirty'],'Worktree is dirty')
        receipt['environment']={'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),'worker_threads':1}
        receipt['ci']={k:os.environ.get(k) for k in ['GITHUB_REPOSITORY','GITHUB_SHA','GITHUB_RUN_ID','GITHUB_RUN_ATTEMPT','GITHUB_WORKFLOW']}
        if os.environ.get('GITHUB_RUN_ID'):
            receipt['ci_run_url']='https://github.com/'+os.environ['GITHUB_REPOSITORY']+'/actions/runs/'+os.environ['GITHUB_RUN_ID']
        source_paths=[HERE/n for n in ['benchmark.py','summarize.py','verify_independent.py','replay_with_receipt.py']]+[ROOT/'.github/workflows/isogeny-invariant-baseline.yml']
        inputs=[HERE/'results/run-001.json.gz',HERE/'results/summary-001.json']
        receipt['source_files']={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
        receipt['input_files']={str(p.relative_to(ROOT)):sha(p) for p in inputs}
        env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None)
        commands=[['benchmark.py','--output',str(out/'fresh.json')],
                  ['summarize.py',str(out/'fresh.json'),'--output',str(out/'summary.json')],
                  ['summarize.py',str(inputs[0]),'--output',str(out/'archive-summary.json')],
                  ['verify_independent.py',str(out/'fresh.json'),'--output',str(out/'independent.json')]]
        for i,c in enumerate(commands):
            argv=[sys.executable,str(HERE/c[0]),*c[1:]];start=time.perf_counter()
            with (out/f'command-{i}.stdout').open('w') as stdout,(out/f'command-{i}.stderr').open('w') as stderr:
                p=subprocess.run(argv,cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
            receipt['commands'].append({'argv':argv,'exit_code':p.returncode,'wall_seconds':time.perf_counter()-start,'stdout':f'command-{i}.stdout','stderr':f'command-{i}.stderr'})
            need(p.returncode==0,'Command failed: '+c[0])
        fresh=json.loads((out/'summary.json').read_text());archive=json.loads((out/'archive-summary.json').read_text());frozen=json.loads(inputs[1].read_text())
        need(fresh==archive==frozen,'Fresh/archive/frozen summaries differ')
        # Negative controls must fail in the independent checker.
        raw=json.loads((out/'fresh.json').read_text());controls=[]
        for mutation in ['pair_count','curve_point','orbit_count']:
            bad=copy.deepcopy(raw)
            if mutation=='pair_count':bad['tiny']['fraction']['counts_in_point_order'][0]+=1
            elif mutation=='curve_point':bad['tiny']['points'][1][1]=(bad['tiny']['points'][1][1]+1)%81
            else:bad['tiny']['fraction']['field_orbits']['1']+=1
            try:verify(bad)
            except ValueError as exc:controls.append({'mutation':mutation,'rejected':True,'reason':str(exc)})
            else:raise RuntimeError('Independent checker accepted corrupt evidence: '+mutation)
        receipt['negative_controls']=controls
        receipt['independent_verification']=json.loads((out/'independent.json').read_text())
        receipt['status']='verified'
        receipt['limits']=['Receipt and hashes are execution/provenance evidence, not formal mathematical proof.',
            'Independent checker covers F81 incidence; F13^7 remains checked by original implementation.',
            'No solver-speed or full-DLP-cost result. Same author wrote both implementations; external review remains valuable.']
    except Exception as exc:
        receipt['status']='failed';receipt['error']=str(exc)
        (out/'failure.txt').write_text(traceback.format_exc())
    finally:
        receipt['finished_utc']=utc()
        receipt['artifacts']={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(out.iterdir()) if p.is_file()}
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({'status':receipt['status'],'receipt':str(out/'receipt.json'),'source_commit':receipt.get('source_commit')}))
    return 0 if receipt['status']=='verified' else 1
if __name__=='__main__':sys.exit(main())
