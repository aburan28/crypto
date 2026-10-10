#!/usr/bin/env python3
"""Freeze, build, and run the half-word protocol under repository CPU isolation."""
import argparse
import contextlib
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
from types import SimpleNamespace
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
REPO=HERE.parent.parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()


def command_receipt(command,stdout,stderr,timeout,env):
    started=time.monotonic_ns()
    with stdout.open('wb') as out,stderr.open('wb') as err:
        proc=subprocess.Popen(command,stdout=out,stderr=err,env=env,start_new_session=True)
        timed_out=False
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out=True
            # SIGINT lets the isolation controller run its affinity-restoration
            # finally blocks. Then reap/kill the whole group, including a worker
            # that might otherwise survive its timed-out controller.
            with contextlib.suppress(ProcessLookupError):os.killpg(proc.pid,signal.SIGINT)
            try:proc.wait(timeout=5)
            except subprocess.TimeoutExpired:pass
            with contextlib.suppress(ProcessLookupError):os.killpg(proc.pid,signal.SIGKILL)
            proc.wait()
    return {'command':command,'exit_code':proc.returncode,'timed_out':timed_out,
        'process_wall_ns':time.monotonic_ns()-started,'wall_scope':'Operational wrapper duration; not a solver performance sample.',
        'stdout_sha256':sha(stdout),'stderr_sha256':sha(stderr)}


def choose_cpus():
    allowed=set(os.sched_getaffinity(0))
    for cpu in sorted(allowed,reverse=True):
        path=Path(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list')
        text=path.read_text().strip() if path.exists() else str(cpu)
        siblings=set()
        for part in text.split(','):
            if '-' in part:
                a,b=map(int,part.split('-'));siblings.update(range(a,b+1))
            else:siblings.add(int(part))
        if siblings<=allowed and siblings!=allowed:return ','.join(map(str,sorted(siblings)))
    raise RuntimeError('No reservable core leaves another CPU available for housekeeping')


def qualified_conditions(path,cpus):
    if not path.exists():return False,None
    rows=[json.loads(line) for line in path.read_text().splitlines()]
    if len(rows)!=1:return False,None
    record=rows[0];run=record.get('run',{})
    good=(record.get('schema')=='isolated-bench/1' and record.get('mode')=='run'
          and record.get('reserved_cpus')==cpus and run.get('cpus')==cpus
          and run.get('exit_status')==0 and run.get('contended') is False
          and isinstance(record.get('host'),dict) and isinstance(record.get('preflight'),dict))
    try:
        pre=record['preflight']
        good=good and run['wall_seconds']>0 and 0<=run['other_cpu_seconds']<=.10*run['wall_seconds']
        good=good and pre['settle']['seconds']==2 and 0<=pre['settle']['other_cpu_seconds']<=.20
        good=good and all(isinstance(pre['conditions'][key],dict) and 0<=pre['conditions'][key]['some']['avg10']<=5 for key in ['psi_cpu','psi_memory'])
        good=good and all(isinstance(run[phase][key],dict) for phase in ['before','after'] for key in ['psi_cpu','psi_memory'])
        good=good and record['host']['logical_cpus']>len(cpus) and bool(record['host']['machine']) and bool(record['host']['kernel'])
    except (KeyError,TypeError):good=False
    return good,record


def wait_for_quiet(isolation,path,max_samples=30):
    spec=importlib.util.spec_from_file_location('frozen_isolation',isolation)
    tool=importlib.util.module_from_spec(spec);spec.loader.exec_module(tool)
    observations=[]
    args=SimpleNamespace(settle=2.0,max_other_cpu=.10,max_psi=5.0)
    for attempt in range(max_samples):
        record={'attempt':attempt+1,'recorded_utc':now()}
        try:
            record.update(accepted=True,preflight=tool.preflight(args,{os.getpid()}))
        except SystemExit as error:
            record.update(accepted=False,reason=str(error))
        observations.append(record);dump(path,observations)
        if record['accepted']:return
    raise RuntimeError(f'Host did not become quiet in {max_samples} preflight samples; no timed worker was launched')


def split_paired(path,cell,out,receipt):
    data=path.read_bytes();lines=data.splitlines(keepends=True)
    markers=[];offset=0
    for line in lines:
        row=json.loads(line)
        if row['type']=='fixture':markers.append((offset,row['mode']))
        offset+=len(line)
    assert len(markers)==2 and markers[0]==(0,'aa') and markers[1][1]=='ab'
    boundary=markers[1][0];stages={}
    for mode,start,end in [('aa',0,boundary),('ab',boundary,len(data))]:
        suffix='-aa' if mode=='aa' else ''
        target=out/(cell+suffix+'.jsonl');target.write_bytes(data[start:end])
        error=out/(cell+suffix+'.stderr');error.write_bytes(b'')
        stages[mode]={**receipt,'stdout_sha256':sha(target),'stderr_sha256':sha(error),'resource_mode':'paired',
            'derived_from':{'file':path.name,'sha256':sha(path),'byte_start':start,'byte_count':end-start}}
    assert (out/(cell+'-aa.jsonl')).read_bytes()+(out/(cell+'.jsonl')).read_bytes()==data
    return stages


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=['discovery','full'],required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--discovery',type=Path)
    parser.add_argument('--cpus')
    args=parser.parse_args()
    if platform.system()!='Linux' or not hasattr(os,'sched_getaffinity') or not Path('/proc/pressure/cpu').exists():
        parser.error('Qualified timings require Linux CPU affinity and /proc pressure records; no unisolated fallback is permitted.')
    cpus_text=args.cpus or choose_cpus()
    cpus=[]
    for part in cpus_text.split(','):
        if '-' in part:
            a,b=map(int,part.split('-'));cpus.extend(range(a,b+1))
        else:cpus.append(int(part))
    cpus=sorted(set(cpus))
    protocol_path=HERE/('protocol_'+args.phase+'.json')
    protocol=json.loads(protocol_path.read_text());assert protocol['schema_version'] in [2,3]
    discovery_binding=None
    if args.phase=='full':
        assert args.discovery is not None,'Full comparison requires sealed qualified discovery'
        discovery=args.discovery.resolve();manifest=json.loads((discovery/'manifest.json').read_text())
        for name,digest in manifest['files'].items():assert sha(discovery/name)==digest,name
        meta=json.loads((discovery/'metadata.json').read_text())
        assert meta['complete'] and meta['qualified'] and meta['phase']=='discovery'
        assert not (discovery/'EXECUTION_STATUS.json').exists(),'A failed execution or analysis is not an admitted discovery source'
        result=json.loads((discovery/'results.json').read_text())
        assert result['qualified'] and result['all_complete'] and result['phase']=='discovery'
        for name,digest in meta['source_hashes'].items():
            if name.endswith('.rs'):assert sha(HERE/name)==digest,'Changed timed source: '+name
        discovery_binding={'path':str(discovery),'manifest_sha256':sha(discovery/'manifest.json'),
            'rust_source_hashes':{n:h for n,h in meta['source_hashes'].items() if n.endswith('.rs')}}
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    names=sorted(p.name for p in HERE.glob('*.rs'))+['REFERENCE_PROTOCOL.json','REFERENCE_FIXTURES.json','SOURCE_LINEAGE.json','analyze.py','reference.py','isolated_run.py']
    for name in names:shutil.copyfile(HERE/name,out/name)
    shutil.copyfile(protocol_path,out/'protocol.json');names.append('protocol.json')
    shutil.copyfile(REPO/'tools/isolated_bench.py',out/'isolated_bench.py');names.append('isolated_bench.py')
    env=dict(os.environ,RAYON_NUM_THREADS='1',OMP_NUM_THREADS='1')
    rustc=shutil.which('rustc');assert rustc
    meta={'schema_version':protocol['schema_version'],'phase':args.phase,'complete':False,'qualified':False,'started_utc':now(),
        'source_hashes':{name:sha(out/name) for name in names},'cpus':cpus,'threads':1,
        'rustc':subprocess.check_output([rustc,'--version','--verbose'],text=True),
        'git_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        'git_status':subprocess.check_output(['git','status','--short'],cwd=REPO,text=True),
        'host':{'system':platform.system(),'machine':platform.machine(),'platform':platform.platform(),'logical_cpus':os.cpu_count()},
        'scope':protocol['scope'],'discovery_binding':discovery_binding}
    cpu_info={}
    allowed_keys={'model name','Processor','Features','flags','CPU implementer','CPU architecture','CPU variant','CPU part','CPU revision','Hardware'}
    for line in Path('/proc/cpuinfo').read_text().splitlines():
        if ':' in line:
            key,value=line.split(':',1);key=key.strip()
            if key in allowed_keys and key not in cpu_info:cpu_info[key]=value.strip()
    meta['host']['cpu_info']=cpu_info
    memory=next(line for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemTotal:'))
    meta['host']['memory_total_kib']=int(memory.split()[1])
    dump(out/'metadata.json',meta);receipts=[]
    isolation=out/'isolated_bench.py'
    try:
        # Keep binaries in the evidence bundle; source and executable digests both
        # describe the actual baseline/candidate image used for the comparison.
        for name,flags in [('baseline',[]),('worker',[]),('tests',['--test'])]:
            entry='baseline_measure.rs' if name=='baseline' else 'measure.rs'
            cmd=[sys.executable,str(isolation),'busy','--',rustc,'--edition','2021','-O',*flags,str(out/entry),'-o',str(out/name)]
            check=command_receipt(cmd,out/(name+'-compile.stdout'),out/(name+'-compile.stderr'),180,env)
            dump(out/(name+'-compile-receipt.json'),check)
            assert check['exit_code']==0 and not check['timed_out'],'Build failure retained'
            meta[name+'_sha256']=sha(out/name)
        cmd=[sys.executable,str(isolation),'busy','--',str(out/'tests'),'--test-threads=1']
        check=command_receipt(cmd,out/'test.stdout',out/'test.stderr',180,env)
        check.update(stdout=(out/'test.stdout').read_text(),stderr=(out/'test.stderr').read_text())
        dump(out/'test_receipt.json',check)
        assert check['exit_code']==0 and not check['timed_out']
        capabilities=subprocess.check_output([str(out/'worker'),'--capabilities'],text=True,env=env)
        (out/'capabilities.json').write_text(capabilities);meta['capabilities']=json.loads(capabilities)
        dump(out/'metadata.json',meta)
        started=time.monotonic()
        quiet_samples=protocol['quiet_wait']['max_samples']
        wait_for_quiet(isolation,out/'startup_quiet.json',quiet_samples)
        cells=[(n,split,seed,family) for n in protocol['variables'] for split in protocol['splits'] for seed in protocol[split+'_seeds'] for family in protocol['families']]
        for index,(n,split,seed,family) in enumerate(cells,1):
            assert time.monotonic()-started<protocol['limits']['campaign_seconds'],'Campaign cap; partial evidence retained'
            cell=f'n{n}-{split}-{seed}-{family}'
            # PSI is a trailing average: a quiet campaign startup does not
            # qualify every later fixture. Wait before launching any sample,
            # retain every observation, then let the locked run check again.
            wait_for_quiet(isolation,out/(cell+'-quiet.json'),quiet_samples)
            assert time.monotonic()-started<protocol['limits']['campaign_seconds'],'Campaign cap after quiet wait; partial evidence retained'
            order_seed=(protocol['order_seed']^(n<<48)^(seed<<8)^protocol['families'].index(family))&((1<<64)-1)
            record={'cell':cell,'n':n,'split':split,'seed':seed,'family':family,'order_seed':order_seed,'stages':{}}
            reps=protocol.get('repetitions_by_n',{}).get(str(n),protocol['repetitions'])
            modes=['paired'] if protocol['schema_version']>=3 else ['aa','ab']
            for mode in modes:
                condition=out/f'{cell}-{mode}-conditions.jsonl'
                worker=[str(out/'worker'),str(n),str(seed),family,str(reps),str(protocol['limits']['nodes']),str(order_seed),mode]
                cmd=[sys.executable,str(isolation),'run','--cpus',cpus_text,'--out',str(condition),'--label',cell+'/'+mode,
                     '--settle','2','--max-other-cpu','0.10','--max-psi','5.0','--',*worker]
                suffix='-'+mode if mode in ['aa','paired'] else ''
                receipt=command_receipt(cmd,out/(cell+suffix+'.jsonl'),out/(cell+suffix+'.stderr'),protocol['limits']['worker_seconds'],env)
                good,conditions=qualified_conditions(condition,cpus)
                receipt.update(worker_command=worker,conditions_file=condition.name,
                    conditions_sha256=sha(condition) if condition.exists() else None,
                    qualified=good and receipt['exit_code']==0 and not receipt['timed_out'])
                record['stages'][mode]=receipt
                if conditions is not None:receipt['worker_peak_rss_bytes']=conditions.get('run',{}).get('max_rss_kib',0)*1024
                dump(out/'receipts.json',receipts+[record])
                assert receipt['qualified'],'Isolation refused, contended, or failed: '+cell+'/'+mode
                if mode=='paired':
                    assert (out/(cell+suffix+'.stderr')).read_bytes()==b''
                    record['paired']=receipt
                    record['stages']=split_paired(out/(cell+suffix+'.jsonl'),cell,out,receipt)
            receipts.append(record);dump(out/'receipts.json',receipts)
            print('completed',cell,index,'of',len(cells),flush=True)
        meta.update(complete=True,qualified=True,ended_utc=now(),campaign_seconds=time.monotonic()-started)
        dump(out/'metadata.json',meta)
        subprocess.run([sys.executable,str(out/'analyze.py'),str(out)],check=True,env=env)
    except BaseException as error:
        dump(out/'EXECUTION_STATUS.json',{'status':'FAILED_OR_UNQUALIFIED','exception_type':type(error).__name__,'message':str(error),
            'recorded_utc':now(),'completed_cells':len(receipts)})
        raise
    finally:
        dump(out/'manifest.json',{'scope':protocol['scope'],'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'}})


if __name__=='__main__':main()
