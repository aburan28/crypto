"""Independent offline final vector4 audit; only owner-confirmed terminal input."""
from pathlib import Path
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
import re
import statistics

HERE=Path(__file__).parent
ROOT=Path('/private/tmp/ecc2k-vector4-sigma-walk-20260912')
SOURCE=ROOT/'source'
MODES=('control','vector1')
UPDATES=201863462912
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def canonical(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def exactly(text,pattern):
    values=re.findall(pattern,text,re.M);assert len(values)==1,(pattern,values)
    return values[0]
def passed(row):assert row['returncode']==0 and row.get('timedOut',False) is False

def markers(name,text,resources):
    vals={'denominator cache':1,'multiply by value':1,'Frobenius network':3,'polynomial chain':1,
          'polynomial state':1,'unrolled inversion':1,'paired products':1,'direct reduction':1,
          'generated product':1,'native carryless multiply':1,'weighted prefix':2,'compact state':1,
          'shared sigma':1,'vector sigma':int(name=='vector1'),'state tile':256}
    for key,value in vals.items():assert exactly(text,'^packed '+re.escape(key)+r': (.*)$')==str(value)
    expect=tuple(str(resources[k]) for k in ('registers','localBytes','sharedBytes'))
    assert exactly(text,r'^packed kernel: (\d+) registers/thread, (\d+) local bytes/thread, (\d+) shared bytes/block, single-product multiplier$')==expect
    assert exactly(text,r'^packed driver reserved shared bytes/block: (\d+), device (\d+)$')==('1024','0')
    assert 'packed composed native:' not in text

def client(name,row,resources,prior=None,allow_auto=False):
    passed(row)
    cmd=row['command'];text=row.get('output',row.get('raw'))
    assert cmd[0]=='/root/'+name and '--packed' in cmd and '--curve' in cmd and cmd[cmd.index('--curve')+1]=='131'
    assert 'MISMATCH' not in text and 'stopping:' not in text and 'warning: could not write checkpoint' not in text
    workers=int(cmd[cmd.index('--threads')+1]);steps=int(cmd[cmd.index('--steps')+1]);launches=int(cmd[cmd.index('--launches')+1])
    if workers==0:
        assert allow_auto;workers=188*256*2
    weight=0 if '--bench' in cmd else int(cmd[cmd.index('--dp-weight')+1])
    markers(name,text,resources)
    assert exactly(text,r'^backend cuda-packed131: (.*)$')==f'{workers} threads x 16 slots x 1 lanes = {workers*16} walks, dp weight {weight}, {steps} steps per launch'
    progress=[]
    for line in text.splitlines():
        if ' iterations ' not in line:continue
        m=re.fullmatch(r'\s*[\d.]+ s\s+[\d.]+ M it/s\s+(\d+) iterations\s+(\d+) dp\s+(\d+) stored\s+(\d+) dropped\s*',line)
        assert m,line;progress.append(tuple(map(int,m.groups())))
    expected=workers*16*steps*launches
    assert progress and progress[-1][0]==expected
    assert all(0<x[0]<=expected and x[0]%(workers*16*steps)==0 and x[3]==0 for x in progress)
    assert all(a[0]<b[0] and a[1]<=b[1] and a[2]<=b[2] for a,b in zip(progress,progress[1:]))
    reloads=re.findall(r'^reloaded (\d+) points from (\d+) file\(s\), (\d+) distinct orbits$',text,re.M)
    prior_stored=0
    if prior is None:assert not reloads
    else:
        assert weight==60 and '--dp-file' in cmd
        assert reloads==[(str(prior['reports']),'1',str(prior['stored']))]
        prior_stored=prior['stored']
    assert all(prior_stored<=x[2]<=prior_stored+x[1] for x in progress)
    final=exactly(text,r'^\s*finished:\s+(\S+) M it/s, (\d+) distinguished points \((\d+) verified against the reference, (\d+) dropped\)$')
    rate=float(final[0]);reports,verified,dropped=map(int,final[1:])
    assert math.isfinite(rate) and rate>0 and dropped==0 and reports==progress[-1][1]
    return {'rate':rate,'updates':expected,'reports':reports,'verified':verified,'stored':progress[-1][2],
            'workers':workers,'steps':steps,'launches':launches,'weight':weight,'progressRows':len(progress)}

def audit(raw_path,expected_sha,session,app):
    owner=read(HERE/'local-process-exit.json')
    assert owner['exitCode']==0 and owner['artifactSha256']==expected_sha and sha(raw_path)==expected_sha
    assert owner['artifactBytes']==Path(raw_path).stat().st_size
    ret=read(str(raw_path)+'.return.json')
    assert ret['valid'] and ret['rawSha256']==expected_sha and ret['artifactBytes']==Path(raw_path).stat().st_size
    assert session==15541 and app=='ap-8wtaL9FddF0MIdQZJffxYc'
    assert sha(HERE/'predispatch-review.json')=='3608fa937bb73648754042eb2a2ca294e4358a69c6fce17154127a21cfa849c1'
    expected=read(ROOT/'plan/expected.json');r=read(raw_path)
    assert r['valid'] and not r.get('error') and r['remoteArtifact']==ret['remoteArtifact']
    for name,digest in read(ROOT/'driver-freeze.json').items():assert sha(ROOT/name)==digest
    sources={str(p.relative_to(SOURCE)):sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file() and p.name not in ('ecc2k130','ecc2k130-cpu')
             and not {'build','__pycache__'}.intersection(p.relative_to(SOURCE).parts) and p.suffix not in ('.pyc','.o')}
    assert len(sources)==160 and canonical(sources)=='5b629d80a63f176807f05e8fc202aee0968291ec7691923a233752b9e9705c21'
    assert r['sourceFiles']==expected['sourceManifests']=={'control':sources,'vector1':sources}
    assert r['sourceSha256ByMode']==expected['sourceManifestCanonicalSha256']
    assert hashlib.sha256((json.dumps(r['build'],indent=2)+'\n').encode()).hexdigest()=='740f1b2300d1070d9f44d50ef33d874d4f8d1c47d1d9b91d06a28c301a99ce76'
    build=r['build'];assert build['valid'] and build['compileAttempts']==8 and build['gpuAllocated'] is build['executableRun'] is False
    assert r['flags']==expected['flags'] and r['modeFlags']==expected['modeFlags'] and r['image']==expected['image']
    assert r['toolchains']=={'control':'13.3.73','vector1':'13.3.73'} and 'V13.3.73' in r['compiler']['output'];passed(r['compiler'])
    assert r['jitEnvironment']['CUDA_DISABLE_PTX_JIT']=='1' and all(r['jitEnvironment'][k] in (None,'','0') for k in ('CUDA_FORCE_PTX_JIT','CUDA_FORCE_JIT'))
    assert r['actualFileHashes']==r['expectedFileHashes']==r['finalFileHashes']
    assert r['gateSha256']==expected['layoutGateSha256']==sha('/private/tmp/test_polytune_states.py')
    for kind in ('binaries','arithmeticBinaries','storageBinaries','sharedProbeBinaries'):assert r[kind]==build[kind]
    spec=importlib.util.spec_from_file_location('independent_prepare',HERE/'prepare_v2.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    for name in MODES:
        raw_code=build['codeInspection'][name]
        assert raw_code['binarySha256']==build['binaries'][name] and r['compiledCodeGates'][name]==raw_code['gate']
        for fn,code in module.native_sections(raw_code['sass']['output']).items():
            words=[[off,f'0x{lo:016x}',f'0x{hi:016x}'] for off,_,lo,hi in code]
            assert canonical(words)==expected['cpuCode'][name]['sections'][fn]['encodingSha256']
    passed(r['gpuBefore']);passed(r['gpuAfter'])
    before,after=r['gpuInventory'],r['gpuInventoryAfter']
    assert before['name']=='NVIDIA RTX PRO 6000 Blackwell Server Edition'
    assert all(before[k]==after[k] for k in ('name','uuid','driverVersion'))
    assert re.fullmatch(r'GPU-[0-9a-fA-F-]{36}',before['uuid'])
    assert r['runtimeResourceCalibrationComplete'] is True
    cal=r['runtimeCalibrationBinding'];assert cal==r['finalRuntimeCalibrationBinding']
    assert hashlib.sha256((json.dumps(cal['payload'],indent=2)+'\n').encode()).hexdigest()==cal['sha256']
    payload=cal['payload'];assert payload['valid'] and payload['records']==r['runtimeCalibrations']
    assert payload['expectedFileSha256']==sha(ROOT/'plan/expected.json') and payload['binarySha256']==r['binaries'] and payload['sourceHashes']==r['sourceSha256ByMode']
    resources={}
    for name,row in r['runtimeCalibrations'].items():
        assert name in MODES
        assert row['binarySha256Before']==row['binarySha256After']==r['binaries'][name]
        assert row['sourceSha256Before']==row['sourceSha256After']==r['sourceSha256ByMode'][name]
        attrs=tuple(map(int,exactly(row['output'],r'^packed kernel: (\d+) registers/thread, (\d+) local bytes/thread, (\d+) shared bytes/block, single-product multiplier$')))
        assert attrs[0:2]==(104,0) and attrs[2] in (1792,2816)
        resources[name]=dict(registers=attrs[0],localBytes=0,sharedBytes=attrs[2])
        parsed=client(name,row,resources[name]);assert parsed['updates']==128 and parsed['reports']==parsed['verified']==parsed['stored']==0
        assert row['calibration']['valid'] and row['calibration']['includedInRanking'] is False and row['calibration']['resources']==resources[name]
    assert set(resources)==set(MODES) and r['calibratedRuntimeResources']==cal['runtimeResources']==resources
    assert r['modeBindings']==r['finalModeBindings']
    for name,binding in r['modeBindings'].items():
        assert binding['sourceSha256']==r['sourceSha256ByMode'][name] and binding['sourceFileCount']==160
        assert binding['flags']==expected['modeFlags'][name] and binding['runtimeResources']==resources[name]
        assert binding['runtimeCalibrationSha256']==cal['sha256']
        for field,key in [('binarySha256','binaries'),('arithmeticBinarySha256','arithmeticBinaries'),('storageBinarySha256','storageBinaries'),('sharedProbeBinarySha256','sharedProbeBinaries')]:assert binding[field]==r[key][name]
    assert len(r['linkedLibraries'])==8
    for row in r['linkedLibraries'].values():passed(row);assert 'not found' not in row['output']
    for name in MODES:
        for key in ('arithmetic','storage','sharedProbe'):
            row=r[key][name];passed(row);assert row['modeBinding']==r['modeBindings'][name]
        arithmetic=r['arithmetic'][name]['output']
        assert exactly(arithmetic,r'^packed arithmetic vector sigma: (.*)$')==str(int(name=='vector1'))
        passes=[line for line in arithmetic.splitlines() if line.startswith('PASS:')]
        assert [int(exactly(line,r'^PASS: (\d+) GPU')) for line in passes]==[3120,2526,18194,18194,1157,6240]
        assert r['storage'][name]['completed']==dict(valid=True,mode=name,compactState=1,batch=16,cases=128,records=297344)
        assert r['sharedProbe'][name]['output'].splitlines()==[f'packed shared sigma probe: 1',f'packed vector sigma probe: {int(name=="vector1")}',
            'PASS: 21 GPU sigma scenarios, 21036 input pairs, global and selected helpers against independent routing',
            'PASS: 114 complete block mask snapshots, 51072 words, output guards and inactive blocks']
        passed(r['integration'][name]);assert r['integration'][name]['modeBinding']==r['modeBindings'][name]
        assert r['integration'][name]['output'].splitlines()==['PASS: DP replay across launch boundaries, restart and checkpoint resume',
            'PASS: byte-identical resume, scalar iteration count, incompatible checkpoint preserved','PASS: overdue walks restart without false distinguished-point reports']
    nested_count=0
    for label,total in [('small',64),('fullBlocks',16384)]:
        row=r['layoutChecks'][label];passed(row);details=row['details']
        assert details['valid'] and details['totalScalarSlots']==total and details['binarySha256']==r['binaries']
        assert row['modeBindings']==r['modeBindings']
        hashes=[];corpus=[]
        for entry in details['results']:
            name=entry['config']['name'];assert name in MODES
            assert entry['config']['threads']*16==total
            prior=None
            for nested in entry['runs']:
                stage=nested['stage']
                checked=client(name,nested,resources[name],prior=prior if stage=='reports2' else None)
                if stage=='reports1':prior=checked
                if stage.startswith('reports'):assert checked['reports']>0 and checked['verified']==checked['reports']
                else:assert checked['reports']==checked['verified']==0
                assert nested['vectorSigmaIdentityGate']['valid'] and nested['vectorSigmaIdentityGate']['runtimeVectorSigma']==int(name=='vector1')
                nested_count+=1
            hashes.append(entry['stateHashes']);corpus.append(entry['sortedCorpusSha256'])
        assert len(hashes)==2 and hashes[0]==hashes[1] and corpus[0]==corpus[1]
    checkpoints=r['crossLayoutCheckpoints'];passed(checkpoints);cp=checkpoints['details']
    assert cp['valid'] and cp['binarySha256']==cp['binarySha256After']==r['binaries']
    assert checkpoints['modeBindings']==r['modeBindings']
    assert len(cp['rows'])==28 and len(cp['comparisons'])==12
    expected_kinds={'uninterrupted':6,'prefix':6,'same-mode-resume':6,'same-batch-cross-mode-resume':6,'cross-batch-prefix':2,'worker-geometry-rejection':2}
    assert Counter(row['kind'] for row in cp['rows'])==Counter(expected_kinds)
    for row in cp['rows']:
        name=Path(row['command'][0]).name
        if row['kind']=='worker-geometry-rejection':
            assert row['returncode']==6 and row['timedOut'] is False and 'incompatible or incomplete' in row['output']
            assert not re.search(r'^resumed from |finished:',row['output'],re.M)
            rejection=row['rejection'];assert rejection['valid'] and rejection['inputSha256Before']==rejection['inputSha256After']==rejection['sourceSha256After']
            markers(name,row['output'],resources[name])
        else:
            checked=client(name,row,resources[name]);assert checked['reports']==checked['verified']==checked['stored']==0
            assert row['checkpoint']['iteration'] in (16,32) and re.fullmatch('[0-9a-f]{64}',row['checkpoint']['sha256'])
    for name in MODES:
        auto=client(name,r['occupancyProbes'][name],resources[name],allow_auto=True);assert auto['workers']==96256 and auto['updates']==1540096
        worker=client(name,r['workerProbes'][name],resources[name]);assert worker['workers']==385024 and worker['updates']==6160384
        assert r['occupancyProbes'][name]['occupancy']['residentBlocksPerSm']==2
    assert r['compatibilityGatesComplete'] is True and r['finalReseedSynchronized'] is True
    assert r['workersByMode']=={'control':385024,'vector1':385024} and r['batchesByMode']=={'control':16,'vector1':16}
    assert (r['scalarWalks'],r['steps'],r['launches'],r['runId'],r['expectedUpdates'])==(6160384,1024,32,1,UPDATES)
    screening=r['screen'];assert len(r['warmup'])==2 and len(screening)==3
    assert [row['name'] for row in screening]==['control','vector1','control']
    threshold=max(screening[0]['rate'],screening[2]['rate'])*1.005
    assert r['qualificationThresholdM']==threshold
    assert r['controlDriftPercent']==100*(screening[2]['rate']/screening[0]['rate']-1)
    qualified=screening[1]['rate']>threshold
    assert (r.get('screenWinner')=='vector1')==qualified
    planned=[('control','warmup',0),('vector1','warmup',0),('control','screen',0),('vector1','screen',1),('control','screen',2)]
    all_samples=list(r['warmup'])+list(screening)
    if qualified:
        for phase in ('confirmation','collection'):
            for repeat in range(3):
                order=(('vector1','control') if repeat%2==0 else MODES) if phase=='confirmation' else (MODES if repeat%2==0 else ('vector1','control'))
                planned += [(name,phase,repeat) for name in order]
            assert set(r[phase])==set(MODES) and all(len(r[phase][name])==3 for name in MODES)
            all_samples += [row for name in MODES for row in r[phase][name]]
    else:
        assert r['confirmation']==r['collection']=={} and 'summary' not in r and 'collectionSummary' not in r
    assert [(row['name'],row['phase'],row['repeat']) for row in r['timedAttempts']]==planned
    assert len(all_samples)==len(planned)==(17 if qualified else 5)
    keyed={(row['name'],row['phase'],row['repeat']):row for row in all_samples};assert len(keyed)==len(planned)
    counts=[];fresh_paths=[];corpus_hashes=set();corpus_sizes=set()
    for attempt in r['timedAttempts']:
        key=(attempt['name'],attempt['phase'],attempt['repeat']);sample=keyed[key]
        assert attempt['command']==sample['command'] and attempt['output']==sample['raw'] and sample['valid']
        checked=client(attempt['name'],attempt,resources[attempt['name']])
        assert checked['updates']==UPDATES and checked['verified']==0 and checked['stored']==checked['reports']
        assert sample['rate']==checked['rate'] and sample['modeBinding']==r['modeBindings'][attempt['name']]
        if attempt['phase']=='collection':
            assert checked['weight']==34 and checked['reports']>0
            assert sample['corpusRecords']==checked['reports'] and sample['corpusBytes']==32*checked['reports']
            assert re.fullmatch('[0-9a-f]{64}',sample['sortedCorpusSha256'])
            fresh_paths.append(attempt['command'][attempt['command'].index('--dp-file')+1])
            corpus_hashes.add(sample['sortedCorpusSha256']);corpus_sizes.add((sample['corpusBytes'],sample['corpusRecords']))
        else:assert checked['weight']==0 and checked['reports']==0
        counts.append(dict(name=attempt['name'],phase=attempt['phase'],repeat=attempt['repeat'],**checked))
    if qualified:assert len(fresh_paths)==len(set(fresh_paths))==6 and len(corpus_hashes)==len(corpus_sizes)==1
    metrics={}
    if qualified:
        for phase,summary_key in [('confirmation','summary'),('collection','collectionSummary')]:
            metrics[phase]={}
            for name in MODES:
                values=[row['rate'] for row in r[phase][name]];summary=r[summary_key][name]
                assert summary['valid'] and summary['samples']==r[phase][name]
                assert summary['rate']==statistics.median(values) and summary['minRate']==min(values) and summary['maxRate']==max(values)
                metrics[phase][name]={'ratesM':values,'medianB':statistics.median(values)/1000,'minimumB':min(values)/1000,'maximumB':max(values)/1000}
            metrics[phase]['gainPercent']=(metrics[phase]['vector1']['medianB']/metrics[phase]['control']['medianB']-1)*100
    labels=[row['label'] for row in r['attemptedCommands']]
    first_timing=next(i for i,label in enumerate(labels) if label.startswith('WARMUP '))
    assert all(labels.index(label)<first_timing for label in ('ARITHMETIC control','ARITHMETIC vector1','STORAGE control','STORAGE vector1','SHARED SIGMA PROBE control','SHARED SIGMA PROBE vector1','INTEGRATION control','INTEGRATION vector1','LAYOUT small','LAYOUT fullBlocks','BATCH GEOMETRY CHECKPOINTS','WORKER PROBE control','WORKER PROBE vector1'))
    report={'valid':True,'status':'terminal_fixed_vector4_comparison_independently_validated','ownerSession':session,'app':app,
            'rawSha256':expected_sha,'rawBytes':Path(raw_path).stat().st_size,'sourceCanonicalSha256':canonical(sources),
            'driverFreezeSha256':sha(ROOT/'driver-freeze.json'),'prebuildSha256':sha(ROOT/'plan/full-client-build.json'),
            'predispatchReviewSha256':sha(HERE/'predispatch-review.json'),'qualified':qualified,'screenRatesM':[row['rate'] for row in screening],
            'screenThresholdM':threshold,'controlDriftPercent':r['controlDriftPercent'],'timedRows':len(counts),'samples':counts,'metrics':metrics,
            'runtimeResources':resources,'calibrationSha256':cal['sha256'],'gpuBefore':before,'gpuAfter':after,
            'validation':{'arithmeticModes':2,'arithmeticSuitesEach':6,'storageCasesEach':128,'storageRecordsEach':297344,
                          'sharedProbeEach':{'scenarios':21,'inputPairs':21036,'snapshots':114,'words':51072},
                          'integrationModes':2,'normalizedSlots':[64,16384],'nestedLayoutClients':nested_count,
                          'checkpointChildren':28,'checkpointSuccesses':26,'expectedRejections':2,'checkpointComparisons':12,
                          'allBeforeTiming':True},
            'corpora':{'timedCollections':len(fresh_paths),'sizes':list(corpus_sizes),'sortedHashes':list(corpus_hashes),'allMatch':qualified},
            'goal':{'targetB':26,'confirmedVectorBenchmarkMedianAtLeastTarget':qualified and metrics['confirmation']['vector1']['medianB']>=26,
                    'confirmedVectorCollectionMedianAtLeastTarget':qualified and metrics['collection']['vector1']['medianB']>=26},
            'findings':[],'limits':['One fixed configuration on one GPU allocation; no population-level performance claim.',
                                 'Compiler selected-path reductions are not runtime speed claims; only the complete scalar samples are timed.',
                                 'Runtime shared attributes are calibrated separately from2816-byte compiler extents.',
                                 'GPU clock/temperature data are before/after snapshots; no counter/profile attribution is established.',
                                 'Returned corpus hashes/counts are validated against the frozen producer gates; temporary corpus/checkpoint payloads are not independently retained here.',
                                 'No knobs, source edits, GPU retries or alternate designs were run.'],
            'independentAuditorSha256':sha(__file__)}
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--artifact',type=Path,required=True);parser.add_argument('--sha256',required=True);parser.add_argument('--session',type=int,required=True);parser.add_argument('--app',required=True);args=parser.parse_args()
    report=audit(args.artifact,args.sha256,args.session,args.app)
    target=HERE/'final-review.json'
    with target.open('x') as handle:json.dump(report,handle,indent=2);handle.write('\n')
    print(json.dumps({'valid':True,'reviewSha256':sha(target),'qualified':report['qualified'],'metrics':report['metrics'],'goal':report['goal']},indent=2))
