"""Execute the already compiled raw128 correctness fixture on one RTX GPU."""

if not __debug__:
    raise RuntimeError("This benchmark requires Python assertions; remove -O, -OO and PYTHONOPTIMIZE.")

from pathlib import Path
import csv,hashlib,json,os,signal,subprocess,time
import modal

HERE=Path(__file__).resolve().parent
RESULT=HERE/'build'
IMAGE='nvidia/cuda@sha256:03c372fd9c65fe7739279f8c65473b315dc61efaaffab03e1e65bc7be7aee61e'
app=modal.App('ecc2k-sparse-hybrid-validation')
image=(modal.Image.from_registry(IMAGE,add_python='3.12').entrypoint([])
       .apt_install('build-essential').env({'CUDA_DISABLE_PTX_JIT':'1'}))
volume=modal.Volume.from_name('ecc2k130',create_if_missing=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

@app.function(image=image,gpu='RTX-PRO-6000',cpu=2,memory=8192,timeout=300,
              block_network=True,retries=0,single_use_containers=True,volumes={'/data':volume})
def validate(compile_receipt,expected_binary_sha,expected_source_sha):
    result=dict(valid=False,kind='raw128 sparse/dense GPU correctness only',timed=False,commands=[])
    def run(label,command,timeout):
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
        timed=False
        try:output,_=p.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed=True
            try:os.killpg(p.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            output,_=p.communicate()
        row=dict(label=label,command=command,returncode=p.returncode,timedOut=timed,output=output)
        result['commands'].append(row)
        print(json.dumps(row),flush=True)
        if p.returncode or timed:raise RuntimeError(label+' failed')
        return output
    try:
        assert compile_receipt['valid'] and compile_receipt['compileAttempts']==1
        volume.reload()
        path=Path(compile_receipt['remoteArtifact'])
        assert str(path).startswith('/data/sparse-hybrid-compile/') and path.name=='result.json'
        assert sha(path)==compile_receipt['rawSha256']
        compiled=json.loads(path.read_text());result['compileReceipt']=compile_receipt
        assert compiled['valid'] and compiled['sourceSha256']==expected_source_sha
        assert compiled['binarySha256']==expected_binary_sha and compiled['gpuAllocated'] is False
        source=path.parent/'probe.cu';binary=path.parent/'probe'
        assert sha(source)==expected_source_sha and sha(binary)==expected_binary_sha
        target=Path('/tmp/sparse-hybrid-probe');target.write_bytes(binary.read_bytes());target.chmod(0o755)
        assert sha(target)==expected_binary_sha
        result['binarySha256']=expected_binary_sha;result['sourceSha256']=expected_source_sha
        result['jitEnvironment']={k:os.environ.get(k) for k in ('CUDA_DISABLE_PTX_JIT','CUDA_FORCE_PTX_JIT','CUDA_FORCE_JIT')}
        assert result['jitEnvironment']['CUDA_DISABLE_PTX_JIT']=='1'
        assert all(result['jitEnvironment'][k] in (None,'','0') for k in ('CUDA_FORCE_PTX_JIT','CUDA_FORCE_JIT'))
        command=['nvidia-smi','--query-gpu=name,uuid,driver_version','--format=csv']
        before=run('GPU before',command,15);rows=list(csv.reader(before.splitlines()))
        assert len(rows)==2 and rows[1][0].strip()=='NVIDIA RTX PRO 6000 Blackwell Server Edition'
        result['gpu']=rows[1]
        result['programOutput']=run('raw128 validation',[str(target)],120)
        summary=json.loads(result['programOutput'].splitlines()[-1])
        expected=dict(valid=True,basisPairs=16384,edgePairs=64,densePairs=64,cases=16512,
                      coefficientComparisons=4227072,rawProductComparisons=33024,timed=False)
        assert summary==expected and 'sm_120, 188 SMs' in result['programOutput']
        result['validation']=summary
        after=run('GPU after',command,15);assert list(csv.reader(after.splitlines()))[1]==rows[1]
        assert sha(target)==sha(binary)==expected_binary_sha and sha(source)==expected_source_sha
        result['valid']=True
    except Exception as exc:result['error']=str(exc)
    folder=Path('/data/sparse-hybrid-validation',str(time.time_ns()));folder.mkdir(parents=True)
    path=folder/'result.json';result['remoteArtifact']=str(path)
    path.write_text(json.dumps(result,indent=2)+'\n');volume.commit()
    return dict(valid=result['valid'],error=result.get('error'),remoteArtifact=str(path),
                rawSha256=sha(path),artifactBytes=path.stat().st_size)

@app.local_entrypoint()
def main():
    RESULT.mkdir(parents=True,exist_ok=True)
    receipt=json.loads((RESULT/'compile-return.json').read_text())
    compiled=json.loads((RESULT/'compile-result.json').read_text())
    assert sha(RESULT/'compile-result.json')==receipt['rawSha256']
    assert sha(HERE/'probe.cu')==compiled['sourceSha256']
    answer=validate.remote(receipt,compiled['binarySha256'],compiled['sourceSha256'])
    (RESULT/'gpu-return.json').write_text(json.dumps(answer,indent=2)+'\n')
    path=RESULT/'gpu-result.json'
    with path.open('wb') as f:
        for chunk in volume.read_file(answer['remoteArtifact'][len('/data'):]):f.write(chunk)
    assert sha(path)==answer['rawSha256'] and path.stat().st_size==answer['artifactBytes']
    print(json.dumps(answer),flush=True)
    if not answer['valid']:raise RuntimeError(answer.get('error','GPU validation failed'))
