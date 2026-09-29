"""Ephemeral one-GPU known-scalar diagnostic; modal run, not modal deploy."""
from pathlib import Path
import modal
from receipt_io import reserve, save

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
def execute(mode: str='smoke'):
    import os
    import subprocess
    import sys
    import tempfile
    sys.path.insert(0,PREFIX+'/gpu')
    from receipt_io import read_partial
    if mode not in ('smoke','benchmark'):raise ValueError(mode)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with tempfile.TemporaryDirectory(prefix='tau-gpu-') as folder:
        output=Path(folder)/'result.json'
        try:
            proc=subprocess.run([sys.executable,PREFIX+'/gpu/run_gpu.py','--mode',mode,'--output',str(output)],
                                capture_output=True,text=True,timeout=540,env=env)
            return {'benchmark':read_partial(output),'returncode':proc.returncode,
                    'stdout':proc.stdout,'stderr':proc.stderr}
        except subprocess.TimeoutExpired:
            return {'status':'process_timeout','timeout_seconds':540,'partial':read_partial(output)}


@app.local_entrypoint()
def main(output: str='modal-gpu-result.json',mode: str='smoke'):
    if mode not in ('smoke','benchmark'):raise ValueError(mode)
    destination=Path(output);destination.parent.mkdir(parents=True,exist_ok=True)
    reserve(destination,{'status':'cloud_pending','mode':mode,'gpu_executed':None})
    try:
        result=execute.remote(mode)
    except Exception as e:
        save(destination,{'status':'cloud_error','mode':mode,'gpu_executed':None,
                          'error_type':type(e).__name__,'error':str(e)})
        raise
    save(destination,result)
    if result.get('benchmark',{}).get('status')!='passed':
        raise RuntimeError('GPU diagnostic did not pass; the failure/partial receipt was retained')
    print('Verified GPU receipt:',destination)
