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
       .add_local_file(HERE.parents[4]/'tools/isolated_bench.py','/study/tools/isolated_bench.py')
       .add_local_file(HERE.parents[4]/'tools/curve_identity.py','/study/tools/curve_identity.py'))


@app.function(image=image,gpu='RTX-PRO-6000',cpu=2,memory=2048,
              timeout=600,retries=0,max_containers=1,scaledown_window=2)
def execute(mode: str='smoke',study: str='original'):
    import sys
    import tempfile
    sys.path.insert(0,PREFIX+'/gpu')
    from receipt_io import read_partial
    from launch_gpu import launch
    with tempfile.TemporaryDirectory(prefix='tau-gpu-') as folder:
        output=Path(folder)/'result.json'
        launch(output,mode,study)
        return read_partial(output)


@app.local_entrypoint()
def main(output: str='modal-gpu-result.json',mode: str='smoke',study: str='original'):
    if mode not in ('smoke','benchmark'):raise ValueError(mode)
    if study not in ('original','square-unroll'):raise ValueError(study)
    destination=Path(output);destination.parent.mkdir(parents=True,exist_ok=True)
    reserve(destination,{'status':'cloud_pending','mode':mode,'study':study,'gpu_executed':None})
    try:
        result=execute.remote(mode,study)
    except Exception as e:
        save(destination,{'status':'cloud_error','mode':mode,'study':study,'gpu_executed':None,
                          'error_type':type(e).__name__,'error':str(e)})
        raise
    save(destination,result)
    if result.get('status')!='passed':
        raise RuntimeError('GPU diagnostic did not pass; the failure/partial receipt was retained')
    print('Verified GPU receipt:',destination)
