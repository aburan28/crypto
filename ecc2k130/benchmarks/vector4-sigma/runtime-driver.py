"""Vector sigma0/1 comparison at fixed shared1/compact1/WP2/B16. Eight CPU deployment builds; first runtime calibrations precede correctness gates and ranking. Goal26B/s."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import math
import os
import re
import signal
import statistics
import subprocess
import sys
import tempfile
import time
import modal

PLAN = Path('/private/tmp/ecc2k-vector4-sigma-walk-20260912/plan')
SOURCE = Path('/private/tmp/ecc2k-vector4-sigma-walk-20260912/source')
DRIVER = Path('/private/tmp/ecc2k-vector4-sigma-walk-20260912/run.py')
IMAGE = 'nvidia/cuda@sha256:03c372fd9c65fe7739279f8c65473b315dc61efaaffab03e1e65bc7be7aee61e'
app = modal.App('ecc2k130-vector4-sigma-walk')
image = (modal.Image.from_registry(IMAGE, add_python='3.12').entrypoint([])
         .apt_install('build-essential').env({'CUDA_DISABLE_PTX_JIT': '1'})
         .add_local_dir(str(SOURCE), remote_path='/root/ecc2k130', copy=True)
         .add_local_file('/private/tmp/test_polytune_states.py', remote_path='/root/test_polytune_states.py', copy=True)
         .add_local_file(str(PLAN/'expected.json'), remote_path='/root/toolchain-expected.json', copy=True)
         .add_local_file(str(PLAN/'gates.py'), remote_path='/root/toolchain_gates.py', copy=True)
         .add_local_file(str(PLAN/'build.py'), remote_path='/root/toolchain-build.py', copy=True)
         .add_local_file(str(PLAN/'checkpoint_gate.py'), remote_path='/root/checkpoint-gate.py', copy=True)
         .add_local_file(str(DRIVER), remote_path='/root/toolchain-driver.py', copy=True)
         .run_commands('python3 /root/toolchain-build.py'))
volume = modal.Volume.from_name('ecc2k130', create_if_missing=True)


@app.function(image=image, cpu=2, memory=4096, timeout=180, retries=0,
              block_network=True, single_use_containers=True, volumes={'/data': volume})
def build_receipt():
    """Export the image's validated CPU build without allocating a GPU."""
    build = json.loads(Path('/root/toolchain-build.json').read_text())
    if build.get('valid') is not True:
        raise RuntimeError('image build did not pass its native/source checks')
    path = Path('/data/vector4-sigma-prebuild', str(time.time_ns())+'.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(build, indent=2)+'\n')
    volume.commit()
    return dict(valid=True, gpuAllocated=False, executableRun=False,
                remoteArtifact=str(path), rawSha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                artifactBytes=path.stat().st_size)


def validate_local_freeze():
    frozen=json.loads((DRIVER.parent/'driver-freeze.json').read_text())
    if not all(hashlib.sha256((DRIVER.parent/name).read_bytes()).hexdigest()==want for name,want in frozen.items()):
        raise RuntimeError('driver freeze changed')
    expected=json.loads((PLAN/'expected.json').read_text())
    actual={str(p.relative_to(SOURCE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SOURCE.rglob('*')) if p.is_file() and p.name not in ('ecc2k130','ecc2k130-cpu') and not {'build','__pycache__'}.intersection(p.relative_to(SOURCE).parts) and p.suffix not in ('.pyc','.o')}
    if any(manifest!=actual for manifest in expected['sourceManifests'].values()):
        raise RuntimeError('runtime source changed')


@app.local_entrypoint()
def compile_only():
    validate_local_freeze()
    local_expected=json.loads((PLAN/'expected.json').read_text())
    if local_expected.get('cpuEvidenceStatus')!='complete':
        raise RuntimeError('CPU native/resource evidence pending; no remote prebuild submitted')
    answer = build_receipt.remote()
    path = PLAN/'full-client-build.json'
    with path.open('wb') as handle:
        for chunk in volume.read_file(answer['remoteArtifact'][len('/data'):]): handle.write(chunk)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=answer['rawSha256'] or path.stat().st_size!=answer['artifactBytes']:
        raise RuntimeError('downloaded prebuild differs')
    (PLAN/'full-client-build-return.json').write_text(json.dumps(answer, indent=2)+'\n')
    print(json.dumps(answer), flush=True)


def captured(command, timeout, cwd):
    process = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
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


@app.function(image=image, gpu='RTX-PRO-6000', cpu=4, memory=8192, timeout=5400, retries=0, single_use_containers=True, block_network=True, volumes={'/data': volume})
def run(expected_files):
    sys.path.insert(0, '/root')
    from toolchain_gates import (MODES, FLAGS, MODE_FLAGS, MARKERS, RESOURCES, COMPILED_RESOURCES, WORKERS_BY_MODE, BATCH_BY_MODE, RESIDENT_BLOCKS, SCALAR_WALKS, STEPS, LAUNCHES,
                                UPDATES, RUN_ID, WALK, INIT, require, sha, canonical_sha, source_identity, validate_sources,
                                gpu_inventory, client_result, timed_result, validate_code, storage_result, shared_probe_result, calibration_result, publish_runtime_calibration, runtime_calibration)
    root = Path('/root/ecc2k130')
    source_roots = {name:root for name in MODES}
    sys.path.insert(0, str(root/'codegen'))
    from benchreport import summarizeSamples
    result = dict(valid=False, kind='complete scalar walk vector sigma0/1 comparison with shared1/compact1/weighted2 fixed; fixed B16/compiler/TILE256/min2 and canonical logical population',
                  image=IMAGE, toolchains={name:'13.3.73' for name in MODES},
                  flags=FLAGS, modeFlags=MODE_FLAGS, expectedRuntimeResources={name:dict(value) for name,value in RESOURCES.items()}, compiledResources=COMPILED_RESOURCES, initialRuntimeResourceStatus='pending first calibration', targetBillionScalarUpdatesPerSecond=26, workersByMode=WORKERS_BY_MODE, batchesByMode=BATCH_BY_MODE, scalarWalks=SCALAR_WALKS,
                  steps=STEPS, launches=LAUNCHES, runId=RUN_ID, expectedUpdates=UPDATES,
                  finalReseedSynchronized=True, compatibilityGatesComplete=False,
                  attemptedCommands=[], runtimeCalibrations={}, arithmetic={}, storage={}, sharedProbe={}, integration={}, layoutChecks={}, occupancyProbes={},
                  workerProbes={}, warmup=[], screen=[], confirmation={}, collection={}, timedAttempts=[])
    def execute(label, command, timeout):
        print(label, flush=True)
        row = captured(command, timeout, root)
        result['attemptedCommands'].append(dict(label=label, **row))
        print(json.dumps(dict(label=label,returncode=row['returncode'],timedOut=row['timedOut'],outputTail=row['output'][-1800:])),flush=True)
        return row
    def gate(label, command, timeout):
        row = execute(label, command, timeout)
        require(row['returncode'] == 0 and not row['timedOut'], label + ' failed: ' + row['output'])
        return row
    try:
        result['expectedFileHashes'] = expected_files
        result['actualFileHashes'] = {name: sha('/root/' + name) for name in expected_files}
        require(result['actualFileHashes'] == expected_files, 'driver/helper/provenance hashes differ')
        expected = json.loads(Path('/root/toolchain-expected.json').read_text())
        require(expected['image'] == IMAGE and expected['toolchains'] == result['toolchains'] and expected['flags'] == FLAGS and expected['modeFlags'] == MODE_FLAGS, 'reviewed configuration differs')
        require(expected.get('cpuEvidenceStatus') == 'complete', 'CPU native/resource evidence is pending')
        result['sourceFiles'] = {name: source_identity(path) for name, path in source_roots.items()}
        result['sourceSha256ByMode'] = validate_sources(result['sourceFiles'], expected)
        result['sourceComparisonScope'] = expected['comparisonScope']
        result['gateSha256'] = sha('/root/test_polytune_states.py')
        require(result['gateSha256'] == expected['layoutGateSha256'], 'layout helper differs')
        result['jitEnvironment'] = {key: os.environ.get(key) for key in ('CUDA_DISABLE_PTX_JIT', 'CUDA_FORCE_PTX_JIT', 'CUDA_FORCE_JIT')}
        require(result['jitEnvironment']['CUDA_DISABLE_PTX_JIT'] == '1' and all(result['jitEnvironment'][k] in (None, '', '0') for k in ('CUDA_FORCE_PTX_JIT', 'CUDA_FORCE_JIT')), 'native-only JIT environment differs')
        result['build'] = json.loads(Path('/root/toolchain-build.json').read_text())
        build = result['build']
        require(build['valid'] is True and build['sourceFiles'] == expected['sourceManifests'], 'incomplete or mismatched image build')
        result['binaries'] = {name: sha('/root/' + name) for name in MODES}
        result['arithmeticBinaries'] = {name: sha('/root/check-' + name) for name in MODES}
        result['storageBinaries'] = {name: sha('/root/storage-' + name) for name in MODES}
        result['sharedProbeBinaries'] = {name: sha('/root/shared-probe-' + name) for name in MODES}
        require(result['sharedProbeBinaries']==build['sharedProbeBinaries'], 'shared-probe binary custody differs')
        require(build['compileAttempts']==8 and build['gpuAllocated'] is False and build['executableRun'] is False, 'exact CPU-only deployment build inventory differs')
        require(result['binaries'] == build['binaries'] and result['arithmeticBinaries'] == build['arithmeticBinaries'] and result['storageBinaries'] == build['storageBinaries'], 'runtime binary custody differs')
        result['compiledCodeGates'] = {}
        for name in MODES:
            row = build['codeInspection'][name]
            require(row['binarySha256'] == result['binaries'][name], 'inspected binary hash differs')
            result['compiledCodeGates'][name] = validate_code(name, row['sass'], row['resources'], expected)
            require(result['compiledCodeGates'][name] == row['gate'], 'build code gate differs on recheck')
        timer_source = (root/'src/main.cu').read_text()
        require('eng.synchronize();\n    const double el = nowSeconds() - t0;' in timer_source, 'final synchronized timer missing')
        result['compiler'] = gate('TOOLCHAIN', ['nvcc', '--version'], 30)
        require('V13.3.73' in result['compiler']['output'], 'runtime image toolchain differs')
        inventory_command = ['nvidia-smi', '--query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,temperature.gpu', '--format=csv']
        result['gpuBefore'] = gate('GPU INVENTORY BEFORE', inventory_command, 30)
        result['gpuInventory'] = gpu_inventory(result['gpuBefore']['output'])
        # First executable client observations: fixed eight-worker, one-step bench.
        # Their elapsed rates are diagnostic output and are never ranked.
        for name in MODES:
            before_sources=validate_sources({mode:source_identity(path) for mode,path in source_roots.items()},expected)
            before_binary=sha('/root/'+name)
            command=['/root/'+name,'--packed','--curve','131','--run-id','1','--threads','8','--steps','1','--launches','1','--verify','0','--bench']
            row=gate('RUNTIME SHARED CALIBRATION '+name,command,120)
            after_sources=validate_sources({mode:source_identity(path) for mode,path in source_roots.items()},expected)
            row.update(binarySha256Before=before_binary,binarySha256After=sha('/root/'+name),sourceSha256Before=before_sources[name],sourceSha256After=after_sources[name])
            row['calibration']=calibration_result(name,row)
            require(before_sources==after_sources and row['binarySha256Before']==row['binarySha256After']==result['binaries'][name],'calibration source/binary changed')
            result['runtimeCalibrations'][name]=row
        result['runtimeCalibrationBinding']=publish_runtime_calibration(result['runtimeCalibrations'],result['binaries'],result['sourceSha256ByMode'])
        require(result['runtimeCalibrationBinding']==runtime_calibration(),'published runtime calibration differs')
        result['runtimeResourceCalibrationComplete']=True
        result['calibratedRuntimeResources']={name:dict(value) for name,value in RESOURCES.items()}

        result['linkedLibraries'] = {}
        for filename in [item for mode in MODES for item in (mode,'check-'+mode,'storage-'+mode,'shared-probe-'+mode)]:
            result['linkedLibraries'][filename] = gate('LINKED LIBRARIES ' + filename, ['ldd', '/root/' + filename], 30)
            require('not found' not in result['linkedLibraries'][filename]['output'], 'linked dependency unavailable')
        def binding(name):
            require(runtime_calibration()==result['runtimeCalibrationBinding'],'runtime calibration changed')
            current_sources = {mode: source_identity(path) for mode, path in source_roots.items()}
            source_hashes = validate_sources(current_sources, expected)
            require(sha('/root/' + name) == result['binaries'][name] and sha('/root/check-' + name) == result['arithmeticBinaries'][name] and sha('/root/storage-' + name) == result['storageBinaries'][name] and sha('/root/shared-probe-' + name)==result['sharedProbeBinaries'][name], 'selected binaries changed')
            return dict(mode=name, toolchain=result['toolchains'][name], binarySha256=result['binaries'][name],
                        arithmeticBinarySha256=result['arithmeticBinaries'][name], storageBinarySha256=result['storageBinaries'][name], sharedProbeBinarySha256=result['sharedProbeBinaries'][name], runtimeCalibrationSha256=result['runtimeCalibrationBinding']['sha256'], sourceSha256=source_hashes[name], sourceFileCount=len(current_sources[name]),
                        flags=MODE_FLAGS[name], runtimeResources=RESOURCES[name], entrypoint='/root/' + name)
        result['modeBindings'] = {name: binding(name) for name in MODES}
        arithmetic_lines = ['packed arithmetic direct reduction: 1', 'packed arithmetic generated product: 1',
            'packed arithmetic native carryless multiply: 1',
            'PASS: 3120 GPU Frobenius vectors, every field basis vector for all selected powers plus dense cases',
            'PASS: 2526 GPU polynomial reductions against long division, including ignored upper-word bits and canonical outputs',
            'PASS: 18194 GPU polynomial products, including all 17161 basis pairs',
            'PASS: 18194 GPU paired polynomial products against independent multiplication',
            'PASS: 1157 GPU polynomial squares against independent multiplication and long division']
        for name in MODES:
            before = binding(name)
            row = gate('ARITHMETIC ' + name, ['/root/check-' + name], 180)
            row['modeBinding'] = binding(name)
            result['arithmetic'][name] = row
            expected_arithmetic = arithmetic_lines[:3] + [f'packed arithmetic weighted prefix: {MODE_FLAGS[name]["PACKED_WEIGHTED_PREFIX"]}',f'packed arithmetic vector sigma: {MODE_FLAGS[name]["PACKED_VECTOR_SIGMA"]}'] + arithmetic_lines[3:] + ['PASS: 6240 GPU paired Frobenius vectors, both inputs against independent routing']
            require(row['modeBinding'] == before and row['output'].splitlines() == expected_arithmetic, 'arithmetic did not complete the required vectors or weighted-prefix marker')
        for name in MODES:
            before = binding(name)
            row = gate('STORAGE ' + name, ['/root/storage-' + name], 180)
            row['modeBinding'] = binding(name)
            result['storage'][name] = row
            row['completed'] = storage_result(name, row)
            require(row['modeBinding'] == before, 'storage validation changed selected code/source')
        for name in MODES:
            before=binding(name)
            row=gate('SHARED SIGMA PROBE '+name,['/root/shared-probe-'+name],180)
            row['modeBinding']=binding(name);row['completed']=shared_probe_result(name,row)
            result['sharedProbe'][name]=row
            require(row['modeBinding']==before,'shared-probe source/binary changed')
        integration_lines = ['PASS: DP replay across launch boundaries, restart and checkpoint resume',
            'PASS: byte-identical resume, scalar iteration count, incompatible checkpoint preserved',
            'PASS: overdue walks restart without false distinguished-point reports']
        print('FULL INTEGRATION: TWO CONFIGURATIONS, AT MOST TWO PROCESSES', flush=True)
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = {name: pool.submit(captured, ['python3', 'codegen/testpackedclient.py', '/root/' + name], 900, root) for name in MODES}
            for name, future in futures.items():
                row = future.result(); row['modeBinding'] = binding(name)
                result['integration'][name] = row
                result['attemptedCommands'].append(dict(label='INTEGRATION ' + name, **row))
                print(json.dumps(dict(label='INTEGRATION '+name,returncode=row['returncode'],timedOut=row['timedOut'],output=row['output'])),flush=True)
                require(row['returncode'] == 0 and not row['timedOut'] and row['output'].splitlines() == integration_lines, 'full integration failed: ' + name)
        arguments = [name+':'+str(BATCH_BY_MODE[name])+':/root/'+name for name in MODES]
        for label, extra, total in [('small', ['--total-slots','64'], 64), ('fullBlocks', ['--total-slots', '16384', '--state-only'], 16384)]:
            output = Path('/tmp/layout-' + label + '.json')
            row = gate('LAYOUT ' + label, ['python3', '/root/test_polytune_states.py', *arguments, *extra, '--output', str(output)], 1200)
            row['details'] = json.loads(output.read_text()); result['layoutChecks'][label] = row
            row['modeBindings'] = {name: binding(name) for name in MODES}
            require(row['details']['valid'] is True and row['details']['totalScalarSlots'] == total and row['details']['binarySha256'] == result['binaries'], 'layout gate identity differs')
            for entry in row['details']['results']:
                name = entry['config']['name']
                for nested in entry['runs']:
                    command = nested['command']
                    nested_steps = int(command[command.index('--steps') + 1])
                    nested_launches = int(command[command.index('--launches') + 1])
                    nested_dp = 0 if '--bench' in command else int(command[command.index('--dp-weight') + 1])
                    nested['vectorSigmaIdentityGate'] = client_result(name, command, nested['returncode'], nested['output'], entry['config']['threads'], nested_steps, nested_launches, nested_dp)
        checkpoint_gate = gate('BATCH GEOMETRY CHECKPOINTS', ['python3', '/root/checkpoint-gate.py', '--output', '/tmp/cross-layout-checkpoints.json'], 900)
        checkpoint_gate['details'] = json.loads(Path('/tmp/cross-layout-checkpoints.json').read_text())
        require(checkpoint_gate['details']['valid'] is True and checkpoint_gate['details']['binarySha256'] == result['binaries'], 'cross-layout checkpoint gate failed')
        checkpoint_gate['modeBindings'] = {name: binding(name) for name in MODES}
        result['crossLayoutCheckpoints'] = checkpoint_gate
        for name in MODES:
            row = gate('OCCUPANCY ' + name, ['/root/' + name, '--packed', '--curve', '131', '--run-id', '1', '--threads', '0', '--steps', '1', '--launches', '1', '--verify', '0', '--bench'], 60)
            result['occupancyProbes'][name] = row; row['modeBinding'] = binding(name)
            device = re.findall(r'^device: (.+), (\d+) SMs, (\d+) block\(s\) of (\d+) packed threads resident per SM$', row['output'], re.MULTILINE)
            blocks=RESIDENT_BLOCKS[name]
            auto_workers=188*256*blocks
            require(device == [(result['gpuInventory']['name'], '188', str(blocks), '256')], 'resident capacity differs')
            row['completed'] = client_result(name, row['command'], row['returncode'], row['output'], auto_workers, 1, 1, 0)
            row['occupancy'] = dict(smCount=188, residentBlocksPerSm=blocks, blockThreads=256, autoWorkers=auto_workers, includedInRanking=False)
            row = gate('WORKER PROBE ' + name, ['/root/' + name, '--packed', '--curve', '131', '--run-id', '1', '--threads', str(WORKERS_BY_MODE[name]), '--steps', '1', '--launches', '1', '--verify', '0', '--bench'], 60)
            result['workerProbes'][name] = row; row['modeBinding'] = binding(name)
            row['completed'] = client_result(name, row['command'], row['returncode'], row['output'], WORKERS_BY_MODE[name], 1, 1, 0)
        result['compatibilityGatesComplete'] = True
        def sample(name, phase, repeat, corpus=None):
            before = binding(name)
            command = ['/root/' + name, '--packed', '--curve', '131', '--run-id', '1', '--threads', str(WORKERS_BY_MODE[name]), '--steps', str(STEPS), '--launches', str(LAUNCHES), '--verify', '0']
            command += ['--bench'] if corpus is None else ['--dp-weight', '34', '--dp-file', str(corpus)]
            raw = execute(phase.upper() + ' ' + name + ' ' + str(repeat), command, 240)
            result['timedAttempts'].append(dict(name=name, phase=phase, repeat=repeat, **raw))
            after = binding(name)
            require(before == after, 'sample changed code or source')
            row = timed_result(name, raw, corpus)
            row.update(phase=phase, repeat=repeat, modeBinding=after)
            return row
        for name in MODES:
            result['warmup'].append(sample(name, 'warmup', 0))
        for index, name in enumerate(('control', *MODES[1:], 'control')):
            result['screen'].append(sample(name, 'screen', index))
        controls = [r['rate'] for r in result['screen'] if r['name'] == 'control']
        result['controlDriftPercent'] = 100 * (controls[-1] / controls[0] - 1)
        result['qualificationThresholdM'] = max(controls) * 1.005
        winner=max((row for row in result['screen'] if row['name']!='control'),key=lambda row:row['rate'])
        if winner['rate'] > result['qualificationThresholdM']:
            result['screenWinner'] = winner['name']
            matched=('control',winner['name'])
            result['confirmation'] = {name: [] for name in matched}
            result['collection'] = {name: [] for name in matched}
            for repeat in range(3):
                for name in (matched if repeat % 2 else tuple(reversed(matched))):
                    result['confirmation'][name].append(sample(name, 'confirmation', repeat))
            result['summary'] = {name: summarizeSamples(rows) for name, rows in result['confirmation'].items()}
            for repeat in range(3):
                for name in (tuple(reversed(matched)) if repeat % 2 else matched):
                    with tempfile.TemporaryDirectory(prefix='shared-sigma-collection-') as directory:
                        result['collection'][name].append(sample(name, 'collection', repeat, Path(directory)/'points.bin'))
            hashes = {row['sortedCorpusSha256'] for rows in result['collection'].values() for row in rows}
            sizes = {(row['corpusBytes'], row['corpusRecords']) for rows in result['collection'].values() for row in rows}
            require(len(hashes) == len(sizes) == 1, 'identical scalar workloads produced different corpora')
            result['collectionSummary'] = {name: summarizeSamples(rows) for name, rows in result['collection'].items()}
        result['gpuAfter'] = gate('GPU INVENTORY AFTER', inventory_command, 30)
        after_gpu = gpu_inventory(result['gpuAfter']['output'])
        require((after_gpu['name'], after_gpu['uuid'], after_gpu['driverVersion']) == (result['gpuInventory']['name'], result['gpuInventory']['uuid'], result['gpuInventory']['driverVersion']), 'GPU changed during comparison')
        result['gpuInventoryAfter'] = after_gpu
        result['finalModeBindings'] = {name: binding(name) for name in MODES}
        require(result['finalModeBindings'] == result['modeBindings'], 'panel source/binary bindings changed')
        result['finalFileHashes'] = {name: sha('/root/' + name) for name in expected_files}
        require(result['finalFileHashes'] == expected_files, 'panel driver/helper/provenance files changed')
        result['finalRuntimeCalibrationBinding']=runtime_calibration()
        require(result['finalRuntimeCalibrationBinding']==result['runtimeCalibrationBinding'],'calibration payload changed after panel')
        result['valid'] = True
    except Exception as exc:
        result['error'] = str(exc)
    path = Path('/data/vector4-sigma-walk', str(time.time_ns()) + '.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    result['remoteArtifact'] = str(path)
    path.write_text(json.dumps(result, indent=2) + '\n'); volume.commit()
    return dict(valid=result['valid'],error=result.get('error'),remoteArtifact=str(path),rawSha256=hashlib.sha256(path.read_bytes()).hexdigest(),artifactBytes=path.stat().st_size)


@app.local_entrypoint()
def main(output: str = '/private/tmp/ecc2k-vector4-sigma-walk-20260912/result.json'):
    validate_local_freeze()
    local_expected=json.loads((PLAN/'expected.json').read_text())
    if local_expected.get('cpuEvidenceStatus')!='complete':
        raise RuntimeError('CPU native/resource evidence pending; no GPU request submitted')
    files = {'toolchain-expected.json': PLAN/'expected.json', 'toolchain_gates.py': PLAN/'gates.py',
             'toolchain-build.py': PLAN/'build.py', 'toolchain-driver.py': DRIVER,
             'checkpoint-gate.py': PLAN/'checkpoint_gate.py'}
    expected_files = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in files.items()}
    answer = run.remote(expected_files)
    destination = Path(output); destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('wb') as handle:
        for chunk in volume.read_file(answer['remoteArtifact'][len('/data'):]): handle.write(chunk)
    if hashlib.sha256(destination.read_bytes()).hexdigest()!=answer['rawSha256'] or destination.stat().st_size!=answer['artifactBytes']: raise RuntimeError('downloaded result differs')
    result=json.loads(destination.read_text())
    Path(output+'.return.json').write_text(json.dumps(answer,indent=2)+'\n')
    print(json.dumps({k: result.get(k) for k in ('valid', 'error', 'screenWinner', 'controlDriftPercent', 'qualificationThresholdM', 'summary', 'collectionSummary', 'remoteArtifact')}), flush=True)
    if not result['valid']:
        raise RuntimeError(result.get('error', 'incomplete comparison'))
