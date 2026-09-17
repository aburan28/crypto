"""Finite client executions with immutable logs and process metadata."""
from types import SimpleNamespace
import json,os
from experiment import ROOT,OUT,finite,sha
def capture(command,log,timeout=300):
 command=[str(v) for v in command]
 row=finite(command,log,timeout,env=dict(os.environ,OMP_NUM_THREADS='8',CUDA_DISABLE_PTX_JIT='1'))
 row['binary_sha256']=sha(__import__('pathlib').Path(command[0]))
 target=log.with_name(log.stem+'-process.json');assert not target.exists();target.write_text(json.dumps(row,indent=2)+'\n')
 return SimpleNamespace(returncode=row['returncode'] if not row['timed_out'] else 124,stdout=log.read_text())
