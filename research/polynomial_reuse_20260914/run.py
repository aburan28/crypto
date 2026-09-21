#!/usr/bin/env python3
"""Matched native S3 experiment. Preserve every timeout, failure and raw output."""
import argparse, hashlib, itertools, json, os, pathlib, resource, socket, subprocess, time, uuid
ROOT=pathlib.Path(__file__).resolve().parents[2]
HERE=pathlib.Path(__file__).resolve().parent

def limit():
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=pathlib.Path,required=True)
    p.add_argument('--reference',type=pathlib.Path,required=True)
    p.add_argument('--redis-server',type=pathlib.Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    contract=json.loads((HERE/'contract.json').read_text())
    (a.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    env={k:v for k,v in os.environ.items() if not k.startswith('IC_')}
    env['RAYON_NUM_THREADS']='1'
    stage=ROOT/'target/release/examples/polynomial_reuse_bench';ic=ROOT/'target/release/ic'
    hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [ROOT/'src/cryptanalysis/polynomial_reuse.rs',ROOT/'src/cryptanalysis/algebra_cache.rs',ROOT/'src/cryptanalysis/koblitz_groebner.rs',ROOT/'src/cryptanalysis/koblitz_index_calculus.rs',ROOT/'examples/polynomial_reuse_bench.rs',HERE/'run.py',HERE/'contract.json']}
    hashes['candidate_binary']=hashlib.sha256(ic.read_bytes()).hexdigest()
    hashes['reference_binary']=hashlib.sha256(a.reference.read_bytes()).hexdigest()
    (a.output/'hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    records=[]
    def run(name,cmd,extra=None):
        start=time.monotonic();status='finished';code=None
        try:
            r=subprocess.run(list(map(str,cmd)),env=env| (extra or {}),capture_output=True,timeout=30,preexec_fn=limit)
            stdout,stderr,code=r.stdout,r.stderr,r.returncode
        except subprocess.TimeoutExpired as e:
            status='timeout';stdout=e.stdout or b'';stderr=e.stderr or b''
        (a.output/(name+'.stdout')).write_bytes(stdout);(a.output/(name+'.stderr')).write_bytes(stderr)
        rec={'name':name,'command':list(map(str,cmd)),'status':status,'returncode':code,'wall_seconds':time.monotonic()-start}
        records.append(rec)
        (a.output/'processes.json').write_text(json.dumps(records,indent=2)+'\n')
        print(name,status,code,flush=True)
        return stdout
    for (n,ell,m),rep,variant in itertools.product(contract['stages']['cases'],range(3),contract['stages']['variants']):
        run(f'stage-n{n}-l{ell}-m{m}-r{rep}-{variant}',[stage,n,ell,m,variant])
    with socket.socket() as s:s.bind(('127.0.0.1',0));port=s.getsockname()[1]
    server=subprocess.Popen([str(a.redis_server),'--bind','127.0.0.1','--port',str(port),'--save','','--appendonly','no','--maxmemory','128mb','--maxmemory-policy','allkeys-lru'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    try:
        for _ in range(100):
            try:
                with socket.create_connection(('127.0.0.1',port),timeout=.1):break
            except OSError:time.sleep(.05)
        else:raise RuntimeError('test Redis failed to start')
        for n,solver,seed in itertools.product(contract['dlp']['degrees'],contract['dlp']['solvers'],contract['dlp']['seeds']):
            name=f'dlp-n{n}-{solver}-s{seed}'
            args=['run','--degree',n,'--solver',solver,'--seed',seed,'--known-log',5,'--batch',1,'--max-trials',500,'--json']
            run(name+'-reference',[a.reference,*args])
            ns='reuse-'+uuid.uuid4().hex
            modes=[('off',{}),('preprocess-local',{'IC_PREPROCESS_CACHE':'local'}),('both-local',{'IC_PREPROCESS_CACHE':'local','IC_REDUCTION_CACHE':'local'}),
                ('both-redis-cold',{'IC_PREPROCESS_CACHE':'redis','IC_REDUCTION_CACHE':'redis'}),('both-redis-warm',{'IC_PREPROCESS_CACHE':'redis','IC_REDUCTION_CACHE':'redis'})]
            for mode,settings in modes:
                run(name+'-'+mode,[ic,*args],settings|{'IC_REDIS_URL':f'redis://127.0.0.1:{port}/0','IC_REDIS_NAMESPACE':ns})
    finally:
        server.terminate();server.wait(timeout=5)
    print('Saved',len(records),'process records',flush=True)
if __name__=='__main__':main()
