"""CPU deployment build of vector4-sigma clients, arithmetic tests and storage tests; no executable or GPU run."""
from pathlib import Path
import json,os,signal,shlex,subprocess,sys
sys.path.insert(0, '/root')
from toolchain_gates import MODE_FLAGS, MODES, WALK, INIT, source_identity, validate_sources, require, sha, validate_code
ROOT=Path('/root/ecc2k130')
SOURCE_ROOTS={name:ROOT for name in MODES}

def captured(command, timeout):
    process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, start_new_session=True)
    try:
        output, _ = process.communicate(timeout=timeout)
        return dict(command=command, returncode=process.returncode, output=output, timedOut=False)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        output, _ = process.communicate()
        return dict(command=command, returncode=process.returncode, output=output, timedOut=True)


expected=json.loads(Path('/root/toolchain-expected.json').read_text())
result=dict(valid=False,kind='fixed-source vector4-sigma full client/arithmetic/storage/shared-probe builds',sourceFiles={name:source_identity(SOURCE_ROOTS[name]) for name in MODES},image=expected['image'],commands={},compiler={},compileAttempts=0,gpuAllocated=False,executableRun=False)
try:
    require(expected.get('cpuEvidenceStatus') == 'complete', 'CPU native/resource evidence pending; no compiler may run')
    result['sourceHashes']=validate_sources(result['sourceFiles'],expected)
    for tool in ('nvcc','ptxas','cuobjdump'):
        row=captured([tool,'--version'],30);result['compiler'][tool]=row
        require(row['returncode']==0 and not row['timedOut'] and 'V13.3.73' in row['output'],'pinned toolchain differs')
    for name in MODES:
        ROOT=SOURCE_ROOTS[name]
        flags=MODE_FLAGS[name];args=[f'{k}={v}' for k,v in flags.items()]
        command=['make','-B','gpu','ARCH=-gencode arch=compute_120,code=sm_120',*args,'GLOBAL_CG=0']
        require(not any(x.startswith('GPUFLAGS=') for x in command),'unexpected Make override')
        item={};result['commands'][name]=item
        dry=captured(['make','-n',*command[1:]],60);item['dryRun']=dry
        require(dry['returncode']==0 and not dry['timedOut'],'Make rendering failed')
        tokens=shlex.split(dry['output'].replace('\\\n',' '))
        require(tokens.count(f'-DECC_PACKED_STATE_TILE={flags["PACKED_STATE_TILE"]}')==1 and '--def-load-cache=cg' not in tokens,'layout flag/cache recipe differs')
        require(tokens.count('-DECC_PACKED_GENERATED_PRODUCT=1')==1 and tokens.count('-DECC_PACKED_CLMAD=1')==1,'arithmetic build flag missing or duplicated')
        require(tokens.count(f'-DECC_PACKED_WEIGHTED_PREFIX={flags["PACKED_WEIGHTED_PREFIX"]}')==1,'weighted-prefix build flag missing or duplicated')
        require(tokens.count(f'-DECC_PACKED_COMPACT_STATE={flags["PACKED_COMPACT_STATE"]}')==1,'vector4-sigma build flag missing or duplicated')
        require(tokens.count(f'-DECC_PACKED_SHARED_SIGMA={flags["PACKED_SHARED_SIGMA"]}')==1,'vector4-sigma flag missing or duplicated')
        require(tokens.count(f'-DECC_PACKED_VECTOR_SIGMA={flags["PACKED_VECTOR_SIGMA"]}')==1,'vector mode flag missing or duplicated')
        result['compileAttempts']+=1
        built=captured(command,1200);item['client']=built
        require(built['returncode']==0 and not built['timedOut'] and dry['output'].strip() in built['output'],'full client build differs')
        target=Path('/root',name);target.write_bytes((ROOT/'ecc2k130').read_bytes());target.chmod(0o755)
        defs=[f'-DECC_{k}={v}' for k,v in flags.items()]
        arithmetic=['nvcc','-O3','-std=c++17','-gencode','arch=compute_120,code=sm_120',*defs,'-Xptxas','-v','src/testpackedcuda.cu','-o','/root/check-'+name]
        result['compileAttempts']+=1
        item['arithmetic']=captured(arithmetic,300)
        require(item['arithmetic']['returncode']==0 and not item['arithmetic']['timedOut'],'arithmetic build failed')
        storage=['nvcc','-O3','-std=c++17','-gencode','arch=compute_120,code=sm_120',*defs,'-Xptxas','-v','src/testpackedstatecuda.cu','-o','/root/storage-'+name]
        result['compileAttempts']+=1
        item['storage']=captured(storage,300)
        require(item['storage']['returncode']==0 and not item['storage']['timedOut'],'storage test build failed')
        shared_probe=['nvcc','-O3','-std=c++17','-gencode','arch=compute_120,code=sm_120',*defs,'-Xptxas','-v','src/testsharedsigmacuda.cu','-o','/root/shared-probe-'+name]
        result['compileAttempts']+=1
        item['sharedProbe']=captured(shared_probe,300)
        require(item['sharedProbe']['returncode']==0 and not item['sharedProbe']['timedOut'],'shared sigma probe build failed')

    result['binaries']={name:sha('/root/'+name) for name in MODES}
    result['arithmeticBinaries']={name:sha('/root/check-'+name) for name in MODES}
    result['storageBinaries']={name:sha('/root/storage-'+name) for name in MODES}
    result['sharedProbeBinaries']={name:sha('/root/shared-probe-'+name) for name in MODES}
    result['codeInspection']={};result['elfDependencies']={}
    for name in MODES:
        sass=captured(['cuobjdump','--dump-sass','--function',WALK+','+INIT,'/root/'+name],120)
        resources=captured(['cuobjdump','--dump-resource-usage','--function',WALK+','+INIT,'/root/'+name],120)
        row=dict(binarySha256=sha('/root/'+name),sass=sass,resources=resources)
        result['codeInspection'][name]=row;row['gate']=validate_code(name,sass,resources,expected)
        for filename in (name,'check-'+name,'storage-'+name,'shared-probe-'+name):
            row=captured(['readelf','-d','/root/'+filename],30);result['elfDependencies'][filename]=row
            require(row['returncode']==0 and not row['timedOut'],'ELF dependency inspection failed')
    result['sourceHashesAfter']=validate_sources({name:source_identity(SOURCE_ROOTS[name]) for name in MODES},expected)
    require(result['sourceHashesAfter']==result['sourceHashes'],'source changed during build')
    require(result['compileAttempts']==8,'exactly eight deployment builds')
    result['valid']=True
except Exception as exc:
    result['error']=str(exc)
Path('/root/toolchain-build.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result.get(k) for k in ('valid','error','binaries','arithmeticBinaries')}),flush=True)
if not result['valid']:raise RuntimeError(result.get('error','incomplete build'))
