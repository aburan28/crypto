"""Pure source, complete-code, binary and sample gates for the vector4 shared-sigma complete-walk comparison."""
from pathlib import Path
import csv
import hashlib
import json
import math
import os
import re

WALK = '_ZN12eccPacked1314walkE10WalkParamsIjEPj'
INIT = '_ZN12eccPacked1314initE10WalkParamsIjEb'
GPU = 'NVIDIA RTX PRO 6000 Blackwell Server Edition'
MODES = ('control', 'vector1')
SCALAR_WALKS = 6160384
STEPS, LAUNCHES, RUN_ID = 1024, 32, 1
UPDATES = SCALAR_WALKS * STEPS * LAUNCHES
BATCH_BY_MODE = {'control':16, 'vector1':16}
WORKERS_BY_MODE = {name:SCALAR_WALKS // batch for name,batch in BATCH_BY_MODE.items()}
RESIDENT_BLOCKS = {'control':2, 'vector1':2}
WEIGHTED_BY_MODE = {'control':2, 'vector1':2}
COMPACT_BY_MODE = {'control':1, 'vector1':1}
SHARED_BY_MODE = {'control':1, 'vector1':1}
VECTOR_BY_MODE = {'control':0, 'vector1':1}
FLAGS = dict(THREADS=256, STREAM_KARAT=0, SMEM_SPILL=0,
             PACKED_SINGLE_PRODUCT=1, PACKED_CACHE_DENOM=1, PACKED_BY_VALUE=1,
             PACKED_PERM_SIGMA=3, PACKED_POLY_CHAIN=1, PACKED_UNROLL_INV=1,
             PACKED_PAIR_PRODUCTS=1, PACKED_POLY_STATE=1, PACKED_DIRECT_REDUCE=1,
             PACKED_GENERATED_PRODUCT=1, PACKED_CLMAD=1, PACKED_STATE_TILE=256, PACKED_WEIGHTED_PREFIX=2)
MODE_FLAGS = {name:dict(FLAGS, BATCH=BATCH_BY_MODE[name], MINBLOCKS=RESIDENT_BLOCKS[name], PACKED_COMPACT_STATE=COMPACT_BY_MODE[name], PACKED_SHARED_SIGMA=SHARED_BY_MODE[name], PACKED_VECTOR_SIGMA=VECTOR_BY_MODE[name]) for name in MODES}
MARKERS = {'denominator cache':1, 'multiply by value':1, 'Frobenius network':3,
           'polynomial chain':1, 'polynomial state':1, 'unrolled inversion':1,
           'paired products':1, 'direct reduction':1, 'native carryless multiply':1}
# Offline compiler extents and runtime function requirements are separate.
REGISTERS_BY_MODE = {'control':104, 'vector1':104}
COMPILED_RESOURCES = {'control':dict(registers=104, localBytes=0, sharedBytes=2816),
                      'vector1':dict(registers=104, localBytes=0, sharedBytes=2816)}
RESOURCES = {'control':dict(registers=104, localBytes=0, sharedBytes=None),
             'vector1':dict(registers=104, localBytes=0, sharedBytes=None)}
CALIBRATION_PATH = Path('/tmp/vector4-sigma-runtime-calibration.json')
CALIBRATION_ENV = 'ECC_VECTOR4_SIGMA_CALIBRATION_SHA256'



def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def source_identity(root):
    root = Path(root)
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() and p.name not in ('ecc2k130', 'ecc2k130-cpu')
            and not set(p.relative_to(root).parts).intersection(('build', '__pycache__'))
            and p.suffix not in ('.pyc', '.o')}


def validate_sources(actual, expected):
    require(set(actual) == set(MODES), 'missing source binding')
    require(actual == expected['sourceManifests'], 'frozen source differs')
    require(all(value == actual['control'] for value in actual.values()), 'all configurations must use the same source')
    require({name: canonical_sha(value) for name, value in actual.items()} == expected['sourceManifestCanonicalSha256'], 'source digest differs')
    require(expected['modeFlags'] == MODE_FLAGS, 'physical layout flags differ')
    return {name: canonical_sha(value) for name, value in actual.items()}


def sass_sections(text):
    sections = {}
    for section in text.split('Function : ')[1:]:
        name = section.splitlines()[0].strip()
        if name not in (WALK, INIT):
            continue
        require(name not in sections, 'duplicate packed code section')
        words, instructions, pending = [], [], None
        for line in section.splitlines():
            m = re.match(r'\s*/\*([0-9a-f]+)\*/\s*(.*?);\s*/\* (0x[0-9a-f]+) \*/', line)
            if m:
                require(pending is None, 'missing high instruction word')
                pending = [int(m[1], 16), m[3]]
                instructions.append([int(m[1], 16), ' '.join(m[2].split())])
                continue
            m = re.match(r'\s*/\* (0x[0-9a-f]+) \*/', line)
            if m and pending is not None:
                words.append(pending + [m[1]])
                pending = None
        require(words and pending is None, 'empty or incomplete code section')
        require([r[0] for r in words] == list(range(0, words[-1][0] + 16, 16)), 'code/helper instruction gap')
        require(all(re.fullmatch(r'0x[0-9a-f]{16}', v) for row in words for v in row[1:]), 'instruction width differs')
        sections[name] = dict(instructionSlots=len(words), encodingSha256=canonical_sha(words),
                              instructionTextSha256=canonical_sha(instructions))
    require(set(sections) == {WALK, INIT}, 'missing packed walk/init or their complete helper sections')
    return sections


def resource_sections(text):
    result = {}
    for name in (WALK, INIT):
        rows = re.findall(r'Function ' + re.escape(name) + r':\s*REG:(\d+) STACK:(\d+) SHARED:(\d+) LOCAL:(\d+)', text)
        require(len(rows) == 1, 'missing or duplicate packed resource row')
        result[name] = dict(zip(('registers', 'stackBytes', 'sharedBytes', 'declaredLocalBytes'), map(int, rows[0])))
    return result


def validate_code(name, sass, resource, expected):
    require(expected.get('cpuEvidenceStatus') == 'complete', 'CPU native code/resource bindings are pending')
    require(name in MODES and sass['returncode'] == resource['returncode'] == 0 and sass.get('timedOut') is False and resource.get('timedOut') is False, 'compiled client inspection failed or timed out')
    for row, flag in ((sass, '--dump-sass'), (resource, '--dump-resource-usage')):
        require(row['command'] == ['cuobjdump', flag, '--function', WALK + ',' + INIT, '/root/' + name], 'wrong inspected client or section filter')
    actual = sass_sections(sass['output'])
    want = expected['cpuCode'][name]
    for function in (WALK, INIT):
        require(actual[function]['encodingSha256'] == want['sections'][function]['encodingSha256']
                and actual[function]['instructionSlots'] == want['sections'][function]['instructionSlots'],
                'actual full-client code differs from reviewed CPU instructions: ' + name + ':' + function)
    resources = resource_sections(resource['output'])
    require(resources == want['resources'], 'actual full-client resources differ from reviewed CPU resources')
    return dict(valid=True, sections=actual, resources=resources, cpuReceiptSha256=want['rawSha256'],
                comparison='Full walk/init instruction encodings, including scoped helper bodies; metadata and pretty-printer text need not be identical.')


def gpu_inventory(text):
    rows = [[v.strip() for v in row] for row in csv.reader(text.splitlines()) if any(v.strip() for v in row)]
    require(len(rows) == 2 and len(rows[0]) == len(rows[1]) == 8 and rows[0][:3] == ['name', 'uuid', 'driver_version'], 'expected exactly one GPU')
    require(rows[1][0] == GPU and re.fullmatch(r'GPU-[0-9a-fA-F]{8}(?:-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}', rows[1][1]), 'GPU identity differs')
    return dict(name=rows[1][0], uuid=rows[1][1], driverVersion=rows[1][2], fields=dict(zip(*rows)))


def _parse_client_result(name, command, returncode, output, workers, steps, launches, dp_weight, runtime_resources):
    require(name in MODES and command[0] == '/root/' + name and returncode == 0, 'client mode or exit status differs')
    require('MISMATCH' not in output and 'stopping:' not in output, 'client interrupted or mismatched')
    require(isinstance(runtime_resources['registers'], int) and isinstance(runtime_resources['sharedBytes'], int), 'selected CPU resource expectation is pending')
    find = lambda pattern: re.findall(pattern, output, re.MULTILINE)
    require(find(r'^packed kernel: (\d+) registers/thread, (\d+) local bytes/thread, (\d+) shared bytes/block, single-product multiplier$')
            == [tuple(str(runtime_resources[k]) for k in ('registers', 'localBytes', 'sharedBytes'))], 'runtime kernel resource/multiplier mismatch')
    for label, value in MARKERS.items():
        require(find(r'^packed ' + re.escape(label) + r': (\d+)$') == [str(value)], 'native flag identity mismatch: ' + label)
    require(find(r'^packed state tile: (\d+)$') == [str(MODE_FLAGS[name]['PACKED_STATE_TILE'])], 'physical state tile marker differs')
    require(find(r'^packed generated product: (.*)$') == ['1'], 'G1 arithmetic marker differs')
    require(find(r'^packed weighted prefix: (.*)$') == [str(WEIGHTED_BY_MODE[name])], 'weighted prefix marker differs')
    require(find(r'^packed compact state: (.*)$') == [str(COMPACT_BY_MODE[name])], 'compact-state marker differs')
    require(find(r'^packed shared sigma: (.*)$') == [str(SHARED_BY_MODE[name])], 'shared-sigma mode marker differs')
    require(find(r'^packed vector sigma: (.*)$') == [str(VECTOR_BY_MODE[name])], 'vector-sigma mode marker differs')
    reserved=find(r'^packed driver reserved shared bytes/block: (\d+), device (\d+)$')
    require(reserved == [('1024','0')], 'driver reserved shared diagnostic differs')
    require(not find(r'^packed composed native: .*?$'), 'private arithmetic variant marker')
    geometry = tuple(map(str, (workers, BATCH_BY_MODE[name], workers * BATCH_BY_MODE[name], dp_weight, steps)))
    require(find(r'^backend cuda-packed131: (\d+) threads x (\d+) slots x 1 lanes = (\d+) walks, dp weight (\d+), (\d+) steps per launch$') == [geometry], 'scalar geometry mismatch')
    counts = list(map(int, find(r'(\d+) iterations\s+\d+ dp')))
    expected = workers * BATCH_BY_MODE[name] * steps * launches
    require(counts and counts[-1] == expected and all(0<x<=expected and x%(workers*BATCH_BY_MODE[name]*steps)==0 for x in counts) and all(a < b for a, b in zip(counts, counts[1:])), 'incomplete or inflated scalar counter')
    progress=find(r'M it/s\s+(\d+) iterations\s+(\d+) dp\s+(\d+) stored\s+(\d+) dropped')
    require(len(progress)==len(counts) and all(int(q[3])==0 for q in progress), 'progress/drop accounting differs')
    reloads=find(r'^reloaded (\d+) points from (\d+) file\(s\), (\d+) distinct orbits$')
    require(len(reloads)<=1, 'duplicate corpus reload')
    if reloads:
        require(dp_weight==60 and '--dp-file' in command and len(reloads)==1 and 0<int(reloads[0][2])<=int(reloads[0][0]) and reloads[0][1]=='1', 'unexpected corpus history')
    prior=int(reloads[0][2]) if reloads else 0
    require(all(prior<=int(q[2])<=prior+int(q[1]) for q in progress), 'stored counters exceed validated process/history bounds')
    finals = find(r'^\s*finished:\s+(\S+) M it/s, (\d+) distinguished points \((\d+) verified against the reference, (\d+) dropped\)$')
    require(len(finals) == 1, 'missing or duplicate completion')
    rate, reports, verified, dropped = finals[0]
    rate = float(rate)
    require(math.isfinite(rate) and rate > 0 and int(dropped) == 0 and int(reports)==int(progress[-1][1]), 'invalid rate or dropped reports')
    return dict(valid=True, name=name, command=command, returncode=returncode, raw=output, rate=rate,
                expectedIterations=expected, reportedIterations=counts[-1], reports=int(reports),
                verified=int(verified), dropped=int(dropped), runtimeResources=dict(runtime_resources), runtimeWeightedPrefix=WEIGHTED_BY_MODE[name], runtimeCompactState=COMPACT_BY_MODE[name], runtimeSharedSigma=SHARED_BY_MODE[name], runtimeVectorSigma=VECTOR_BY_MODE[name])




def calibration_result(name, raw):
    command=['/root/'+name,'--packed','--curve','131','--run-id','1','--threads','8','--steps','1','--launches','1','--verify','0','--bench']
    require(name in MODES and raw['command']==command and raw['returncode']==0 and raw.get('timedOut') is False, 'untimed calibration executable/geometry failed')
    matches=re.findall(r'^packed kernel: (\d+) registers/thread, (\d+) local bytes/thread, (\d+) shared bytes/block, single-product multiplier$',raw['output'],re.MULTILINE)
    require(len(matches)==1,'calibration function attributes cardinality')
    regs,local,shared=map(int,matches[0])
    require(regs==REGISTERS_BY_MODE[name] and local==0,'calibration registers/local differ')
    reserved=re.findall(r'^packed driver reserved shared bytes/block: (\d+), device (\d+)$',raw['output'],re.MULTILINE)
    require(SHARED_BY_MODE=={'control':1,'vector1':1} and reserved==[('1024','0')], 'both calibration modes require shared1 and driver reservation/device0')
    require(shared==2816 or shared+1024==2816,'calibration static requirement does not match either explicit compiled-extent convention')
    resource=dict(registers=regs,localBytes=local,sharedBytes=shared)
    parsed=_parse_client_result(name,command,raw['returncode'],raw['output'],8,1,1,0,resource)
    require(parsed['reportedIterations']==128 and parsed['reports']==parsed['verified']==parsed['dropped']==0,'calibration work/count differs')
    return dict(valid=True,mode=name,observedFunctionSharedBytes=shared,observedDeviceReservedBytes=1024,compiledSharedExtentBytes=COMPILED_RESOURCES[name]['sharedBytes'],resources=resource,includedInRanking=False,completed=parsed)


def _validate_calibration_payload(payload):
    require(payload.get('valid') is True and set(payload['records'])==set(MODES),'complete calibration inventory')
    expected_path=Path('/root/toolchain-expected.json')
    require(payload['expectedFileSha256']==sha(expected_path),'calibration compiled/source expectation changed')
    expected=json.loads(expected_path.read_text())
    require(payload['sourceHashes']==expected['sourceManifestCanonicalSha256'],'calibration source identity')
    build=json.loads(Path('/root/toolchain-build.json').read_text())
    require(build['valid'] is True and payload['binarySha256']==build['binaries'],'calibration prebuilt binary binding')
    resources={}
    for name,record in payload['records'].items():
        require(record['binarySha256Before']==record['binarySha256After']==payload['binarySha256'][name]==sha('/root/'+name),'calibration executable changed')
        require(record['sourceSha256Before']==record['sourceSha256After']==payload['sourceHashes'][name],'calibration source changed')
        parsed=calibration_result(name,record)
        require(record['calibration']==parsed,'calibration raw/parsed binding')
        resources[name]=dict(parsed['resources'])
    return resources


def publish_runtime_calibration(records,binaries,source_hashes):
    require(not CALIBRATION_PATH.exists() and CALIBRATION_ENV not in os.environ,'calibration must be a fresh first observation')
    payload=dict(valid=True,records=records,binarySha256=binaries,sourceHashes=source_hashes,expectedFileSha256=sha('/root/toolchain-expected.json'))
    resources=_validate_calibration_payload(payload)
    with CALIBRATION_PATH.open('x') as handle:handle.write(json.dumps(payload,indent=2)+'\n')
    os.environ[CALIBRATION_ENV]=sha(CALIBRATION_PATH)
    RESOURCES.update(resources)
    return dict(path=str(CALIBRATION_PATH),sha256=sha(CALIBRATION_PATH),payload=payload,runtimeResources=resources)


def runtime_calibration():
    expected=os.environ.get(CALIBRATION_ENV)
    require(isinstance(expected,str) and re.fullmatch('[0-9a-f]{64}',expected) and CALIBRATION_PATH.is_file() and sha(CALIBRATION_PATH)==expected,'runtime resources are pending or calibration changed')
    payload=json.loads(CALIBRATION_PATH.read_text());resources=_validate_calibration_payload(payload)
    RESOURCES.update(resources)
    return dict(path=str(CALIBRATION_PATH),sha256=expected,payload=payload,runtimeResources=resources)


def client_result(name, command, returncode, output, workers, steps, launches, dp_weight):
    resources=runtime_calibration()['runtimeResources']
    return _parse_client_result(name,command,returncode,output,workers,steps,launches,dp_weight,resources[name])


def shared_probe_result(name,raw):
    require(name in MODES and raw['command']==['/root/shared-probe-'+name] and raw['returncode']==0 and raw.get('timedOut') is False,'shared sigma probe executable/exit differs')
    want=[f'packed shared sigma probe: {SHARED_BY_MODE[name]}',
          f'packed vector sigma probe: {VECTOR_BY_MODE[name]}',
          'PASS: 21 GPU sigma scenarios, 21036 input pairs, global and selected helpers against independent routing',
          'PASS: 114 complete block mask snapshots, 51072 words, output guards and inactive blocks']
    require(raw['output'].splitlines()==want,'shared sigma probe mode/coverage differs')
    return dict(valid=True,mode=name,sharedSigma=SHARED_BY_MODE[name],vectorSigma=VECTOR_BY_MODE[name],scenarios=21,inputPairs=21036,blockSnapshots=114,maskWords=51072)


def storage_result(name, raw):
    require(name in MODES and raw['command'] == ['/root/storage-' + name], 'storage executable identity differs')
    require(raw['returncode'] == 0 and raw.get('timedOut') is False, 'storage validation failed or timed out')
    expected = [f'packed storage compact state: {COMPACT_BY_MODE[name]}',
                'packed storage batch: 16',
                'PASS: 128 GPU storage cases, 297344 records, independent physical images and logical reads with canaries']
    require(raw['output'].splitlines() == expected, 'storage mode, batch or complete validation counts differ')
    return dict(valid=True, mode=name, compactState=COMPACT_BY_MODE[name], batch=16,
                cases=128, records=297344)


def timed_result(name, raw, corpus=None):
    require(raw.get('timedOut') is False, 'timed sample timed out or lacks explicit completion status')
    command = raw['command']
    prefix = ['/root/' + name, '--packed', '--curve', '131', '--run-id', '1', '--threads', str(WORKERS_BY_MODE[name]),
              '--steps', str(STEPS), '--launches', str(LAUNCHES), '--verify', '0']
    if corpus is None:
        require(command == prefix + ['--bench'], 'timed benchmark argv differs')
    else:
        require(command == prefix + ['--dp-weight', '34', '--dp-file', str(corpus)], 'timed collection argv differs')
    row = client_result(name, command, raw['returncode'], raw['output'], WORKERS_BY_MODE[name], STEPS, LAUNCHES, 0 if corpus is None else 34)
    require(row['expectedIterations'] == UPDATES and row['verified'] == 0, 'timed workload differs')
    if corpus is None:
        require(row['reports'] == 0, 'benchmark emitted reports')
    else:
        data = Path(corpus).read_bytes() if Path(corpus).exists() else b''
        require(row['reports'] > 0 and len(data) == row['reports'] * 32, 'collection file/count mismatch')
        row.update(corpusBytes=len(data), corpusRecords=len(data)//32,
                   sortedCorpusSha256=hashlib.sha256(b''.join(sorted(data[i:i+32] for i in range(0, len(data), 32)))).hexdigest())
    return row
