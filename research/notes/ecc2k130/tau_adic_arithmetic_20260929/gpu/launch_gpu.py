"""Bounded local/container entrypoint for synthetic arithmetic device studies."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
from receipt_io import reserve,save,read_partial

HERE=Path(__file__).resolve().parent
TOOLS=HERE.parents[4]/'tools'


def select_cpus(allowed,siblings):
    for cpu in sorted(allowed):
        group=siblings(cpu)
        if group and group<allowed:return ','.join(map(str,sorted(group)))
    raise RuntimeError('No complete SMT sibling group with another CPU left free; timing isolation unavailable')


def timing_eligible(isolation):
    return (isolation.get('schema')=='isolated-bench/1'
            and isolation.get('mode')=='run'
            and isolation.get('run',{}).get('exit_status')==0
            and isolation.get('run',{}).get('contended') is False
            and isolation.get('left_on_reserved',{}).get('user_threads')==[])


def command_for(mode,study,benchmark,conditions,cpus=None):
    command=[sys.executable,str(HERE/'run_gpu.py'),'--mode',mode,'--study',study,'--output',str(benchmark)]
    if mode=='benchmark':
        if not cpus:raise ValueError('Timing requires reserved CPUs')
        command=[sys.executable,str(HERE/'isolated_entry.py'),'run','--cpus',cpus,
                 '--out',str(conditions),'--label','tau-known-scalar:'+study,'--',*command]
    return command


def bounded(command,env,seconds=540):
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
                             env=env,start_new_session=True)
    def stop_group():
        try:os.killpg(process.pid,signal.SIGTERM)
        except ProcessLookupError:pass
        try:stdout,stderr=process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            try:os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            stdout,stderr=process.communicate()
        return stdout,stderr
    timed_out=False
    try:
        stdout,stderr=process.communicate(timeout=seconds)
    except subprocess.TimeoutExpired:
        timed_out=True;stdout,stderr=stop_group()
    except BaseException:
        stop_group()
        raise
    return {'returncode':process.returncode,'stdout':stdout,'stderr':stderr,'timed_out':timed_out}


def launch(output,mode='smoke',study='original'):
    if mode not in ('smoke','benchmark') or study not in ('original','square-unroll'):
        raise ValueError('Unknown mode or study')
    result={'status':'started','mode':mode,'study':study,'gpu_executed':None,
            'timing_eligible':False,'timeout_seconds':540,'termination_grace_seconds':10}
    reserve(output,result)
    try:
        with tempfile.TemporaryDirectory(prefix='tau-device-') as folder:
            benchmark=Path(folder)/'benchmark.json';conditions=Path(folder)/'conditions.jsonl';cpus=None
            if mode=='benchmark':
                sys.path.insert(0,str(TOOLS))
                import isolated_bench
                cpus=select_cpus(set(os.sched_getaffinity(0)),isolated_bench.smt_siblings)
                result['isolation_tool_sha256']=hashlib.sha256((TOOLS/'isolated_bench.py').read_bytes()).hexdigest()
            result['reserved_cpus']=cpus
            command=command_for(mode,study,benchmark,conditions,cpus)
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
            save(output,result)
            result.update(bounded(command,env))
            result['benchmark']=read_partial(benchmark)
            result['gpu_executed']=result['benchmark'].get('gpu_executed')
            result['isolation']=read_partial(conditions) if mode=='benchmark' else None
            verified=result['returncode']==0 and result['benchmark'].get('status')=='passed'
            result['timing_eligible']=(mode=='benchmark' and not result['timed_out']
                                       and verified and timing_eligible(result['isolation']))
            if result['timed_out']:result['status']='process_timeout'
            elif not verified:result['status']='failed'
            elif mode=='benchmark' and not result['timing_eligible']:result['status']='isolation_rejected'
            else:result['status']='passed'
    except Exception as error:
        result.update(status='failed',error_type=type(error).__name__,error=str(error))
    save(output,result)
    print(json.dumps({k:result[k] for k in ('status','mode','study','gpu_executed','timing_eligible')}))
    return 0 if result['status']=='passed' else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('smoke','benchmark'),default='smoke')
    parser.add_argument('--study',choices=('original','square-unroll'),default='original')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    raise SystemExit(launch(args.output,args.mode,args.study))
