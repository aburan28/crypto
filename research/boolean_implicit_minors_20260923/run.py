#!/usr/bin/env python3
"""Run the fixed, standalone Boolean implicit-minor protocol once into a new directory."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import time
import sys
sys.dont_write_bytecode=True

HERE = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def child(command: list[str], out: Path, err: Path, timeout: float) -> dict:
    started = time.perf_counter_ns()
    with out.open("wb") as stdout, err.open("wb") as stderr:
        proc = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        timed_out = False
        while True:
            pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if (time.perf_counter_ns() - started) / 1e9 >= timeout:
                timed_out = True
                proc.kill()
                _, status, usage = os.wait4(proc.pid, 0)
                break
            time.sleep(0.01)
        proc.returncode = os.waitstatus_to_exitcode(status)
    return {
        "command": command,
        "exit_code": proc.returncode,
        "timed_out": timed_out,
        "process_wall_ns": time.perf_counter_ns() - started,
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
        "peak_rss_bytes": usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024),
        "meter": "wait4 for this fresh worker, not cumulative RUSAGE_CHILDREN",
        "stdout_sha256": sha(out),
        "stderr_sha256": sha(err),
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    out=parser.parse_args().out.resolve();out.mkdir(parents=True,exist_ok=False)
    names=sorted(p.name for p in HERE.glob('*.rs'))+['REFERENCE_PROTOCOL.json','SOURCE_LINEAGE.json','protocol.json','run.py','analyze.py','reference.py']
    for name in names:shutil.copyfile(HERE/name,out/name)
    protocol=json.loads((out/'protocol.json').read_text())
    rustc=shutil.which('rustc');assert rustc
    metadata={'schema_version':1,'complete':False,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_hashes':{name:sha(out/name) for name in names},'rustc':subprocess.check_output([rustc,'--version','--verbose'],text=True),
        'git_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True).strip(),
        'git_status':subprocess.check_output(['git','status','--short'],cwd=HERE,text=True),
        'host':{'system':platform.system(),'machine':platform.machine(),'platform':platform.platform(),'cpus':os.cpu_count()},'scope':protocol['scope']}
    receipts=[];dump(out/'metadata.json',metadata)
    try:
        with tempfile.TemporaryDirectory(prefix='minor-scratch-',dir=HERE) as directory:
            tmp=Path(directory)
            for name,flags in [('worker',[]),('tests',['--test'])]:
                compiled=subprocess.run([rustc,'--edition','2021','-O',*flags,str(out/'measure.rs'),'-o',str(tmp/name)],capture_output=True,timeout=120)
                (out/(name+'-compile.txt')).write_bytes(compiled.stdout+compiled.stderr)
                assert compiled.returncode==0,'Compilation failed; retained'
                metadata[name+'_sha256']=sha(tmp/name)
            so,se=tmp/'test.stdout',tmp/'test.stderr';check=child([str(tmp/'tests')],so,se,120)
            check.update(stdout=so.read_text(),stderr=se.read_text());dump(out/'test_receipt.json',check)
            assert check['exit_code']==0 and not check['timed_out']
            dump(out/'metadata.json',metadata);started=time.monotonic()
            for n in protocol['variables']:
                for seed in protocol['seeds']:
                    for family in protocol['families']:
                        assert time.monotonic()-started<protocol['limits']['campaign_seconds'],'Campaign cap; partial evidence retained'
                        cell=f'n{n}-{seed}-{family}';order_seed=(protocol['order_seed']^(n<<48)^(seed<<8)^protocol['families'].index(family))&((1<<64)-1)
                        command=[str(tmp/'worker'),str(n),str(seed),family,'8',str(protocol['limits']['nodes']),str(order_seed)]
                        so,se=out/(cell+'.jsonl'),out/(cell+'.stderr');receipt=child(command,so,se,protocol['limits']['worker_seconds'])
                        receipt.update(cell=cell,n=n,seed=seed,family=family,order_seed=order_seed);receipts.append(receipt);dump(out/'receipts.json',receipts)
                        assert receipt['exit_code']==0 and not receipt['timed_out'],'Producer failure retained: '+cell
                        print('completed',cell,len(receipts),'of 12',flush=True)
            metadata.update(complete=True,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),campaign_seconds=time.monotonic()-started)
            dump(out/'metadata.json',metadata)
        subprocess.run([sys.executable,str(out/'analyze.py'),str(out)],check=True)
    except BaseException as error:
        dump(out/'EXECUTION_STATUS.json',{'status':'FAILED','exception_type':type(error).__name__,'message':str(error),
            'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'completed_workers':len(receipts)})
        raise
    finally:
        dump(out/'manifest.json',{'scope':protocol['scope'],'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'}})


if __name__=='__main__':main()

