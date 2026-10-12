#!/usr/bin/env python3
"""Bounded discovery-only comparison; never overwrites an evidence directory."""
import argparse
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bounded_runner', HERE / 'run.py')
runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    source=list(HERE.glob('*.rs'))+[HERE/'run.py',Path(__file__)]
    for path in source:shutil.copyfile(path,out/path.name)
    retained=json.loads((HERE.parent/'boolean_initial_restriction_20260923/protocol.json').read_text())['variants']
    additions=['fiber_rows','fiber_columns','fiber_simd','fiber_zero_simd']
    protocol={'scope':'Discovery only; generated public Boolean systems, no holdouts or promotion.',
        'variables':[16,24],'seeds':[17,937],'families':['planted','cross_planted','unplanted'],
        'repetitions':3,'variants':retained+additions,'retained':retained,'node_cap':200000,'worker_seconds':300}
    runner.dump(out/'protocol.json',protocol)
    metadata={'complete':False,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'git_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True).strip(),
        'source_hashes':{p.name:runner.sha(out/p.name) for p in source},
        'rustc':subprocess.check_output(['rustc','--version','--verbose'],text=True),
        'host':{'platform':platform.platform(),'machine':platform.machine()}}
    runner.dump(out/'metadata.json',metadata)
    receipts=[];summary=[]
    with tempfile.TemporaryDirectory(prefix='boolean-fiber-probe-') as temp:
        tmp=Path(temp)
        for name,flags in [('worker',[]),('tests',['--test'])]:
            command=['rustc','--edition','2021','-O',*flags,str(out/'worker.rs'),'-o',str(tmp/name)]
            built=subprocess.run(command,capture_output=True,timeout=90)
            (out/(name+'-compile.txt')).write_bytes(built.stdout+built.stderr)
            assert built.returncode==0,'Compilation failed; retained source and diagnostics'
            metadata[name+'_sha256']=runner.sha(tmp/name)
        so,se=tmp/'test.stdout',tmp/'test.stderr'
        checked=runner.child([str(tmp/'tests')],so,se,120)
        checked.update(stdout=so.read_text(),stderr=se.read_text());runner.dump(out/'test_receipt.json',checked)
        assert checked['exit_code']==0 and not checked['timed_out']
        for n in protocol['variables']:
            for seed in protocol['seeds']:
                for family in protocol['families']:
                    cell=f'n{n}-{seed}-{family}'
                    so,se=out/(cell+'.jsonl'),out/(cell+'.stderr')
                    command=[str(tmp/'worker'),str(n),str(seed),family,str(protocol['repetitions']),str(protocol['node_cap']),'with-final']
                    receipt=runner.child(command,so,se,protocol['worker_seconds'])
                    receipt['cell']=cell;receipts.append(receipt);runner.dump(out/'receipts.json',receipts)
                    assert receipt['exit_code']==0 and not receipt['timed_out'],cell
                    rows=[json.loads(line) for line in so.read_text().splitlines()]
                    fixture=rows[0];samples=rows[1:]
                    assert len(samples)==len(protocol['variants'])*protocol['repetitions']
                    assert all(s['verified'] and s['outcome']!='UNKNOWN' for s in samples)
                    for sample in samples:
                        if sample['outcome']=='SAT':
                            assert all(sum((sample['model']&m)==m for m in p)%2==0 for p in fixture['polys'])
                        else:assert fixture['search_reference']=='UNSAT'
                    costs={arm:statistics.median(s['total_ns'] for s in samples if s['variant']==arm) for arm in protocol['variants']}
                    work={arm:{k:next(s for s in samples if s['variant']==arm)[k] for k in ['fiber_dimension','fiber_prefixes','fiber_queries','fiber_zero_rejected','fiber_linear_rejected']} for arm in additions}
                    ratios={arm:statistics.median(min(s['total_ns'] for s in samples if s['rep']==rep and s['variant'] in retained)/next(s['total_ns'] for s in samples if s['rep']==rep and s['variant']==arm) for rep in range(protocol['repetitions'])) for arm in additions}
                    summary.append({'cell':cell,'n':n,'seed':seed,'family':family,'cost_ns':costs,'retained_ratio':ratios,'new_work':work,'verified':True})
                    print(cell,{k:round(ratios[k],3) for k in ['fiber_rows','fiber_columns','fiber_simd','fiber_zero_simd']},flush=True)
    metadata.update(complete=True,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    runner.dump(out/'metadata.json',metadata)
    runner.dump(out/'summary.json',{'scope':protocol['scope'],'observations':len(summary)*len(protocol['variants'])*protocol['repetitions'],'cells':summary,'promotion':None})
    runner.dump(out/'manifest.json',{'files':{str(p.relative_to(out)):runner.sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}})

if __name__=='__main__':main()
