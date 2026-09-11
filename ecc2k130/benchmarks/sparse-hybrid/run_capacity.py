"""Bounded paired raw-MMA capacity measurement, not a field/walk benchmark."""
from pathlib import Path
import csv,hashlib,json,math,os,signal,statistics,subprocess,time
import modal

HERE=Path(__file__).resolve().parent
RESULT=HERE/'build'
IMAGE='nvidia/cuda@sha256:03c372fd9c65fe7739279f8c65473b315dc61efaaffab03e1e65bc7be7aee61e'
app=modal.App('ecc2k-sparse-hybrid-capacity')
image=(modal.Image.from_registry(IMAGE,add_python='3.12').entrypoint([])
       .apt_install('build-essential').env({'CUDA_DISABLE_PTX_JIT':'1'})
       .add_local_file(str(RESULT/'capacity-code-review.json'),'/root/code-review.json',copy=True))
volume=modal.Volume.from_name('ecc2k130',create_if_missing=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

@app.function(image=image,gpu='RTX-PRO-6000',cpu=2,memory=8192,timeout=300,
              block_network=True,retries=0,single_use_containers=True,volumes={'/data':volume})
def measure(compile_receipt,review_sha):
    result=dict(valid=False,kind='paired rawMMA capacity only',commands=[],walkBenchmark=False)
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
        print(json.dumps(dict(label=label,returncode=p.returncode,timedOut=timed,outputBytes=len(output))),flush=True)
        if p.returncode or timed:raise RuntimeError(label+' failed')
        return output
    try:
        assert sha('/root/code-review.json')==review_sha
        review=json.loads(Path('/root/code-review.json').read_text())
        assert review['valid'] and review['compileRawSha256']==compile_receipt['rawSha256']
        assert compile_receipt['valid'] and compile_receipt['compileAttempts']==1
        volume.reload();path=Path(compile_receipt['remoteArtifact'])
        assert str(path).startswith('/data/sparse-hybrid-capacity-compile/') and path.name=='result.json'
        assert sha(path)==compile_receipt['rawSha256']
        compiled=json.loads(path.read_text());result['compileReceipt']=compile_receipt
        assert compiled['valid'] and compiled['sourceHashes']==review['sourceHashes']
        assert compiled['binarySha256']==review['binarySha256']
        for name,value in review['sourceHashes'].items():assert sha(path.parent/name)==value
        binary=path.parent/'probe';assert sha(binary)==review['binarySha256']
        target=Path('/tmp/sparse-hybrid-capacity-probe');target.write_bytes(binary.read_bytes());target.chmod(0o755)
        assert sha(target)==review['binarySha256']
        result.update(binarySha256=review['binarySha256'],sourceHashes=review['sourceHashes'],codeReviewSha256=review_sha)
        result['jitEnvironment']={k:os.environ.get(k) for k in ('CUDA_DISABLE_PTX_JIT','CUDA_FORCE_PTX_JIT','CUDA_FORCE_JIT')}
        assert result['jitEnvironment']['CUDA_DISABLE_PTX_JIT']=='1'
        assert all(result['jitEnvironment'][k] in (None,'','0') for k in ('CUDA_FORCE_PTX_JIT','CUDA_FORCE_JIT'))
        query=['nvidia-smi','--query-gpu=name,uuid,driver_version','--format=csv']
        before=run('GPU before',query,15);inventory=list(csv.reader(before.splitlines()))
        assert len(inventory)==2 and inventory[1][0].strip()=='NVIDIA RTX PRO 6000 Blackwell Server Edition'
        result['gpu']=inventory[1]
        output=run('correctness and capacity',[str(target),'--bench'],180)
        result['programOutput']=output
        assert 'sm_120, 188 SMs' in output
        records=[json.loads(line) for line in output.splitlines() if line.startswith('{')]
        checks=[r for r in records if 'basisPairs' in r]
        assert checks==[dict(valid=True,basisPairs=16384,edgePairs=64,densePairs=64,cases=16512,
                             coefficientComparisons=4227072,rawProductComparisons=33024,timed=False)]
        result['correctness']=checks[0]
        resources=[r for r in records if r.get('kind')=='resources']
        modes=['dense4','hybrid4','dense8','hybrid8']
        assert [r['mode'] for r in resources]==modes
        expected_regs={'dense4':40,'hybrid4':38,'dense8':56,'hybrid8':52}
        for row in resources:
            assert row['registers']==expected_regs[row['mode']] and row['localBytes']==row['sharedBytes']==0
            assert row['occupancyApiBlocksPerSm']>0
        result['resources']=resources
        rows=[r for r in records if r.get('kind')=='sample']
        order=[(name,'warmup',i) for name in modes for i in range(2)]
        for pair in range(2):
            for repeat in range(5):
                order += [(modes[2*pair+(repeat&1)],'measurement',repeat),
                          (modes[2*pair+1-(repeat&1)],'measurement',repeat)]
        assert [(r['mode'],r['phase'],r['repeat']) for r in rows]==order
        for row in rows:
            chains=4 if row['mode'].endswith('4') else 8
            sparse=row['mode'].startswith('hybrid')
            cores=6016*1024*chains
            assert row['valid'] and row['workers']==192512 and row['warps']==6016
            assert row['rounds']==1024 and row['chains']==chains
            assert row['rawCoreEquivalents']==cores and row['outputsChecked']==192512*chains*4
            assert row['denseMatrixInstructions']==cores*(1 if sparse else 3)
            assert row['sparseMatrixInstructions']==(cores if sparse else 0)
            assert 0<=row['maximumOutput']<=1124139008
            assert math.isfinite(row['milliseconds']) and row['milliseconds']>=0.05
            seconds=row['milliseconds']/1000
            row.update(rawCoreEquivalentsPerSecond=cores/seconds,
                       denseMatrixInstructionsPerSecond=row['denseMatrixInstructions']/seconds,
                       sparseMatrixInstructionsPerSecond=row['sparseMatrixInstructions']/seconds)
        result['samples']=rows
        summaries={}
        for name in modes:
            values=[r['rawCoreEquivalentsPerSecond'] for r in rows if r['phase']=='measurement' and r['mode']==name]
            assert len(values)==5
            summaries[name]=dict(median=statistics.median(values),minimum=min(values),maximum=max(values))
        result['rawCoreRateSummaries']=summaries
        result['hybridToDenseRatios']={str(c):summaries['hybrid'+str(c)]['median']/summaries['dense'+str(c)]['median'] for c in (4,8)}
        final=[r for r in records if r.get('kind')=='capacitySummary']
        assert len(final)==1 and final[0]['valid'] and final[0]['warmups']==8 and final[0]['measurements']==20
        assert final[0]['workers']==192512 and final[0]['rounds']==1024 and len(records)==34
        result['summary']=final[0]
        after=run('GPU after',query,15);assert list(csv.reader(after.splitlines()))[1]==inventory[1]
        assert sha(target)==sha(binary)==review['binarySha256']
        for name,value in review['sourceHashes'].items():assert sha(path.parent/name)==value
        result['limits']=['Per-launch event rates with full CPU verification between launches; not uninterrupted sustained capacity.',
                          'Includes device fragment load, seed initialization and output store; excludes host packing/checks, parity reconstruction, field reduction, inversion and walk.',
                          'Raw-core equivalents and matrix instruction rates are not complete scalar walk rates or universal hardware ceilings.']
        result['valid']=True
    except Exception as exc:result['error']=str(exc)
    folder=Path('/data/sparse-hybrid-capacity',str(time.time_ns()));folder.mkdir(parents=True)
    path=folder/'result.json';result['remoteArtifact']=str(path)
    path.write_text(json.dumps(result,indent=2)+'\n');volume.commit()
    return dict(valid=result['valid'],error=result.get('error'),remoteArtifact=str(path),rawSha256=sha(path),artifactBytes=path.stat().st_size)

@app.local_entrypoint()
def main():
    RESULT.mkdir(parents=True,exist_ok=True)
    receipt=json.loads((RESULT/'capacity-compile-return.json').read_text())
    assert sha(RESULT/'capacity-compile-result.json')==receipt['rawSha256']
    answer=measure.remote(receipt,sha(RESULT/'capacity-code-review.json'))
    (RESULT/'capacity-gpu-return.json').write_text(json.dumps(answer,indent=2)+'\n')
    path=RESULT/'capacity-gpu-result.json'
    with path.open('wb') as f:
        for chunk in volume.read_file(answer['remoteArtifact'][len('/data'):]):f.write(chunk)
    assert sha(path)==answer['rawSha256'] and path.stat().st_size==answer['artifactBytes']
    result=json.loads(path.read_text())
    print(json.dumps({k:result.get(k) for k in ('valid','error','rawCoreRateSummaries','hybridToDenseRatios','limits')}),flush=True)
    if not answer['valid']:raise RuntimeError(answer.get('error','capacity measurement failed'))
