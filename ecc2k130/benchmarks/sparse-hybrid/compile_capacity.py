"""One bounded CPU-only CUDA compile; no executable or GPU run."""

if not __debug__:
    raise RuntimeError("This benchmark requires Python assertions; remove -O, -OO and PYTHONOPTIMIZE.")

from pathlib import Path
import hashlib,json,os,signal,subprocess,time
import modal
import capacity_gate

HERE=Path(__file__).resolve().parent
RESULT=HERE/'build'
IMAGE='nvidia/cuda@sha256:03c372fd9c65fe7739279f8c65473b315dc61efaaffab03e1e65bc7be7aee61e'
app=modal.App('ecc2k-sparse-hybrid-capacity-compile')
image=(modal.Image.from_registry(IMAGE,add_python='3.12').entrypoint([])
       .apt_install('build-essential')
       .add_local_file(str(HERE/'probe.cu'),'/root/probe.cu',copy=True)
       .add_local_file(str(HERE/'bench.cu'),'/root/bench.cu',copy=True)
       .add_local_file(str(HERE/'capacity_gate.py'),'/root/capacity_gate.py',copy=True))
volume=modal.Volume.from_name('ecc2k130',create_if_missing=True)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

@app.function(image=image,cpu=2,memory=4096,timeout=360,block_network=True,
              retries=0,single_use_containers=True,volumes={'/data':volume})
def build(source_hashes):
    folder=Path('/data/sparse-hybrid-capacity-compile',str(time.time_ns()));folder.mkdir(parents=True)
    result=dict(valid=False,gpuAllocated=False,executableRun=False,compileAttempts=0,
                image=IMAGE,expectedSourceHashes=source_hashes,sourceHashes={n:sha('/root/'+n) for n in source_hashes},commands=[])
    def run(label,command,timeout):
        started=time.monotonic()
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
        timed=False
        try:output,_=p.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed=True
            try:os.killpg(p.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            output,_=p.communicate()
        row=dict(label=label,command=command,returncode=p.returncode,timedOut=timed,
                 elapsedSeconds=time.monotonic()-started,output=output)
        result['commands'].append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='output'}),flush=True)
        if p.returncode or timed:raise RuntimeError(label+' failed')
        return row
    try:
        assert result['sourceHashes']==source_hashes
        for tool in ('nvcc','ptxas','cuobjdump'):
            row=run(tool+' version',[tool,'--version'],15)
            assert 'V13.3.73' in row['output']
        result['compileAttempts']=1
        run('compile',['nvcc','-O3','-std=c++17','-arch=sm_120','-lineinfo','-Xptxas','-v',
                       '/root/bench.cu','-o','/tmp/sparse-hybrid-capacity-probe'],120)
        run('SASS',['cuobjdump','--dump-sass','--function','denseKernel,hybridKernel,dense4,hybrid4,dense8,hybrid8','/tmp/sparse-hybrid-capacity-probe'],30)
        run('resources',['cuobjdump','--dump-resource-usage','--function','denseKernel,hybridKernel,dense4,hybrid4,dense8,hybrid8','/tmp/sparse-hybrid-capacity-probe'],30)
        run('ELF dependencies',['readelf','-d','/tmp/sparse-hybrid-capacity-probe'],15)
        result['binarySha256']=sha('/tmp/sparse-hybrid-capacity-probe')
        (folder/'probe').write_bytes(Path('/tmp/sparse-hybrid-capacity-probe').read_bytes())
        result['binaryArtifact']=str(folder/'probe')
        assert {n:sha('/root/'+n) for n in source_hashes}==source_hashes
        result['valid']=True
    except Exception as exc:result['error']=str(exc)
    for n in source_hashes:(folder/n).write_bytes(Path('/root/'+n).read_bytes())
    path=folder/'result.json';result['remoteArtifact']=str(path)
    path.write_text(json.dumps(result,indent=2)+'\n');volume.commit()
    return dict(valid=result['valid'],error=result.get('error'),remoteArtifact=str(path),
                rawSha256=sha(path),artifactBytes=path.stat().st_size,
                compileAttempts=result['compileAttempts'],gpuAllocated=False)

@app.local_entrypoint()
def main():
    RESULT.mkdir(parents=True,exist_ok=True)
    source_hashes={n:sha(HERE/n) for n in ('probe.cu','bench.cu')}
    receipt=build.remote(source_hashes)
    (RESULT/'capacity-compile-return.json').write_text(json.dumps(receipt,indent=2)+'\n')
    path=RESULT/'capacity-compile-result.json'
    with path.open('wb') as f:
        for chunk in volume.read_file(receipt['remoteArtifact'][len('/data'):]):f.write(chunk)
    assert sha(path)==receipt['rawSha256'] and path.stat().st_size==receipt['artifactBytes']
    print(json.dumps(receipt),flush=True)
    if not receipt['valid']:raise RuntimeError(receipt.get('error','compile failed'))
    capacity_gate.bind_review(path,RESULT/'capacity-code-review.json')
