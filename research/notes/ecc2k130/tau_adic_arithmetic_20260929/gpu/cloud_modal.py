"""Ephemeral one-GPU known-scalar diagnostic; modal run, not modal deploy."""
import json
from pathlib import Path
import modal

HERE=Path(__file__).resolve().parent
PREFIX='/study/research/notes/ecc2k130/tau_adic_arithmetic_20260929'
app=modal.App('tau-adic-known-scalar-diagnostic')
image=(modal.Image.from_registry('nvidia/cuda:12.8.1-devel-ubuntu22.04',add_python='3.11')
       .apt_install('libnuma1','util-linux','procps')
       .pip_install('cupy-cuda12x==13.6.0','numpy==2.2.6')
       .add_local_dir(HERE.parent,PREFIX,ignore=['**/__pycache__/**','**/*.log','gpu/results/**'])
       .add_local_file(HERE.parents[4]/'tools/curve_identity.py','/study/tools/curve_identity.py'))


@app.function(image=image,gpu='RTX-PRO-6000',cpu=2,memory=2048,
              timeout=600,retries=0,max_containers=1,scaledown_window=2)
def execute():
    import os
    import subprocess
    import sys
    output=Path('/tmp/tau-gpu-result.json')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    try:
        proc=subprocess.run([sys.executable,PREFIX+'/gpu/run_gpu.py','--output',str(output)],
                            capture_output=True,text=True,timeout=540,env=env)
        result=json.loads(output.read_text()) if output.exists() else {'status':'failed','gpu_executed':False}
        return {'benchmark':result,'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr}
    except subprocess.TimeoutExpired:
        partial=json.loads(output.read_text()) if output.exists() else {'gpu_executed':False}
        return {'status':'process_timeout','timeout_seconds':540,'partial':partial}


@app.local_entrypoint()
def main(output: str='modal-gpu-result.json'):
    destination=Path(output);destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('x') as f:f.write('{}\n')
    try:
        result=execute.remote()
    except Exception as e:
        destination.write_text(json.dumps({'status':'cloud_error','gpu_executed':False,
                                           'error_type':type(e).__name__,'error':str(e)},indent=2)+'\n')
        raise
    destination.write_text(json.dumps(result,indent=2)+'\n')
    if result.get('benchmark',{}).get('status')!='passed':
        raise RuntimeError('GPU diagnostic did not pass; the failure/partial receipt was retained')
    print('Verified GPU receipt:',destination)
