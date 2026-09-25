#!/usr/bin/env python3
"""Run the fixed, standalone Boolean schedule protocol once into a new directory."""
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
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--protocol',type=Path,default=HERE/'protocol.json')
    args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=json.loads(args.protocol.read_text())
    names=['worker.rs','kernel.rs','protocol.json','run.py','analyze.py']
    for name in names:shutil.copyfile(args.protocol if name=='protocol.json' else HERE/name,out/name)
    rustc=shutil.which('rustc');assert rustc
    metadata={'schema_version':1,'complete':False,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_hashes':{name:sha(out/name) for name in names},'rustc':subprocess.check_output([rustc,'--version','--verbose'],text=True),
              'git_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True).strip(),
              'git_status':subprocess.check_output(['git','status','--short'],cwd=HERE,text=True),
              'host':{'system':platform.system(),'machine':platform.machine(),'platform':platform.platform(),'cpus':os.cpu_count()},'scope':protocol['scope']}
    dump(out/'metadata.json',metadata);receipts=[]
    with tempfile.TemporaryDirectory(prefix='boolean-solve-trace-') as tmp:
        build=Path(tmp);executable=build/'worker';tests=build/'tests'
        for target,flags in [(executable,[]),(tests,['--test'])]:
            result=subprocess.run([rustc,'--edition','2021','-O',*flags,str(out/'worker.rs'),'-o',str(target)],capture_output=True,timeout=60)
            (out/(target.name+'-compile.txt')).write_bytes(result.stdout+result.stderr)
            if result.returncode:raise SystemExit('compile failed; retained evidence')
        metadata['worker_sha256']=sha(executable);metadata['test_executable_sha256']=sha(tests)
        so,se=build/'test.stdout',build/'test.stderr';receipt=child([str(tests)],so,se,30)
        receipt.update(stdout=so.read_text(),stderr=se.read_text());dump(out/'test_receipt.json',receipt)
        if receipt['exit_code']!=0:raise SystemExit('tests failed; experiment not launched')
        dump(out/'metadata.json',metadata);started=time.monotonic()
        for n in protocol['variables']:
            with (out/f'raw-n{n}.jsonl').open('x') as raw:
                for split in ['discovery','holdout']:
                    for seed in protocol[split+'_seeds']:
                        for family in protocol['families']:
                            if time.monotonic()-started>protocol['limits']['campaign_seconds']:raise SystemExit('campaign cap reached; incomplete evidence retained')
                            cell=f'n{n}-{split}-{seed}-{family}'
                            command=[str(executable),str(n),str(seed),family,str(protocol['repetitions']),str(protocol['limits']['nodes'])]
                            so,se=build/'cell.stdout',build/'cell.stderr';receipt=child(command,so,se,protocol['limits']['worker_seconds'])
                            receipt.update(cell=cell,n=n,split=split,seed=seed,family=family);receipts.append(receipt);dump(out/'receipts.json',receipts)
                            if receipt['exit_code']!=0:
                                shutil.copyfile(so,out/'failed.stdout');shutil.copyfile(se,out/'failed.stderr');raise SystemExit(cell+' producer failure; evidence retained')
                            for line in so.read_text().splitlines(keepends=True):
                                row=json.loads(line);row.update(cell=cell,split=split,raw_line=line);raw.write(json.dumps(row,separators=(',',':'))+'\n')
                            raw.flush()
                print(f'completed n={n}: {len(receipts)} fixed cells',flush=True)
        metadata.update(complete=True,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),campaign_seconds=time.monotonic()-started)
        dump(out/'metadata.json',metadata)
    subprocess.run([shutil.which('python3'),str(out/'analyze.py'),str(out)],check=True)
    dump(out/'manifest.json',{'files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'},'scope':protocol['scope']})

if __name__=='__main__':main()
