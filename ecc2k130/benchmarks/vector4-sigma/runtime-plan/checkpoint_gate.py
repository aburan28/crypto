"""Bounded checkpoint gate for vector-sigma modes at fixed compact B16 geometry.

Run only after the separately bound executable correctness gates. Checkpoint v2
stores logical SoA coordinates and requires the original batch/worker geometry.
This gate never treats cross-batch checkpoint byte equality as a requirement.
"""
from pathlib import Path
import argparse
import json
import os
import re
import signal
import subprocess
import sys
import tempfile

sys.path.insert(0, '/root')
from toolchain_gates import MODES, BATCH_BY_MODE, MODE_FLAGS, require, sha, client_result
import test_polytune_states as states


STEPS = 16
RUN_ID = 1
THREAD_CASES = (8, 256, 257)
SMALL_WALKS = 64


def captured(command, timeout):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
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


def state_from_file(path, name, threads, iteration):
    batch = BATCH_BY_MODE[name]
    states.TOTAL = threads * batch
    data = path.read_bytes()
    parsed = states.normalized_checkpoint(data, dict(threads=threads, batch=batch), RUN_ID, iteration)
    for index, seed, dead, start, _, _ in parsed[2]:
        require(seed == states.initial_seed(RUN_ID, index) and dead == 0 and start == 0,
                'no-report checkpoint changed logical seed/dead/start state')
    return dict(fileBytes=len(data), sha256=sha(path), normalizedSha256=states.digest_state(parsed),
                logicalThreads=threads, batch=batch, scalarSlots=threads*batch, iteration=iteration)


def checkpoint_gate(result):
    require(tuple(MODES) == ('control', 'vector1'),
            'checkpoint panel modes differ')
    require(BATCH_BY_MODE == dict(control=16, vector1=16),
            'checkpoint panel batch geometry differs')
    require(all(MODE_FLAGS[name]['BATCH'] == BATCH_BY_MODE[name]
                and MODE_FLAGS[name]['THREADS'] == 256
                and MODE_FLAGS[name]['PACKED_STATE_TILE'] == 256
                and MODE_FLAGS[name]['PACKED_COMPACT_STATE'] == 1
                and MODE_FLAGS[name]['PACKED_WEIGHTED_PREFIX'] == 2
                and MODE_FLAGS[name]['PACKED_SHARED_SIGMA'] == 1
                and MODE_FLAGS[name]['PACKED_VECTOR_SIGMA'] == (0 if name == 'control' else 1)
                for name in MODES),
            'checkpoint mode flags differ')
    result.update(binarySha256={name:sha('/root/'+name) for name in MODES},
                  batchesByMode=dict(BATCH_BY_MODE), rows=[], comparisons=[])
    with tempfile.TemporaryDirectory(prefix='batch-checkpoint-') as directory:
        root = Path(directory)

        def invoke(name, threads, launches, path, kind):
            before = sha('/root/'+name)
            command = ['/root/'+name, '--packed', '--curve', '131', '--run-id', str(RUN_ID),
                       '--threads', str(threads), '--steps', str(STEPS), '--launches', str(launches),
                       '--verify', '0', '--bench', '--checkpoint', str(path)]
            raw = captured(command, 120)
            row = dict(mode=name, threads=threads, batch=BATCH_BY_MODE[name], kind=kind, **raw)
            result['rows'].append(row)
            require(not raw['timedOut'], 'checkpoint child timed out')
            require(before == sha('/root/'+name) == result['binarySha256'][name],
                    'checkpoint client bytes changed')
            return row

        def success(name, threads, launches, path, kind, resume_at=None):
            row = invoke(name, threads, launches, path, kind)
            require(row['returncode'] == 0, 'checkpoint client did not complete')
            # The gate owns validation/loading of the parent's calibration,
            # bound by ECC_VECTOR4_SIGMA_CALIBRATION_SHA256 in this subprocess.
            # Never derive an expected shared-memory size from this child.
            complete = client_result(name, row['command'], row['returncode'], row['output'],
                                     threads, STEPS, launches, 0)
            require(complete['reports'] == complete['verified'] == complete['dropped'] == 0,
                    'no-report checkpoint control emitted or dropped reports')
            resumes = re.findall(r'^resumed from .* at iteration (\d+)$', row['output'], re.MULTILINE)
            require(resumes == ([] if resume_at is None else [str(resume_at)]),
                    'checkpoint resumed at an unexpected iteration')
            require('warning: could not write checkpoint' not in row['output'], 'checkpoint write failed')
            row['completed'] = complete
            row['checkpoint'] = state_from_file(path, name, threads, (resume_at or 0)+STEPS*launches)
            row['expectedPhysicalThreads'] = ((threads+255)//256)*256
            return row

        def rejection(first, second, source, source_threads, target_threads, kind):
            source_before = sha(source)
            rejected = root / (kind+'-'+first+'-to-'+second+'.ck')
            rejected.write_bytes(source.read_bytes())
            before = sha(rejected)
            row = invoke(second, target_threads, 1, rejected, kind)
            require(row['returncode'] == 6, 'incompatible checkpoint was not rejected with exit 6')
            require(('checkpoint '+str(rejected)+' is incompatible or incomplete;') in row['output'],
                    'missing incompatible-checkpoint diagnostic')
            require(not re.search(r'^resumed from |^\s*finished:', row['output'], re.MULTILINE),
                    'incompatible checkpoint unexpectedly resumed or completed')
            require(before == sha(rejected) == source_before == sha(source),
                    'rejected checkpoint or its source bytes changed')
            row['rejection'] = dict(valid=True, sourceMode=first, sourceThreads=source_threads,
                sourceBatch=BATCH_BY_MODE[first], targetThreads=target_threads,
                targetBatch=BATCH_BY_MODE[second], inputSha256Before=before,
                inputSha256After=sha(rejected), sourceSha256After=sha(source),
                sourceScalarSlots=source_threads*BATCH_BY_MODE[first],
                targetScalarSlots=target_threads*BATCH_BY_MODE[second])

        prefixes = {}
        for threads in THREAD_CASES:
            whole, prefix = {}, {}
            for name in MODES:
                whole[name] = root / f'{threads}-{name}-whole.ck'
                prefix[name] = root / f'{threads}-{name}-prefix.ck'
                success(name, threads, 2, whole[name], 'uninterrupted')
                success(name, threads, 1, prefix[name], 'prefix')
                split = root / f'{threads}-{name}-split.ck'
                split.write_bytes(prefix[name].read_bytes())
                success(name, threads, 1, split, 'same-mode-resume', resume_at=STEPS)
                require(split.read_bytes() == whole[name].read_bytes(),
                        'same-geometry resumed checkpoint differs from uninterrupted reference')
                result['comparisons'].append(dict(kind='same-mode split/resume', mode=name,
                    threads=threads, batch=BATCH_BY_MODE[name], finalSha256=sha(split)))
                prefixes[threads, name] = prefix[name]
            for first in MODES:
                for second in MODES:
                    if first == second or BATCH_BY_MODE[first] != BATCH_BY_MODE[second]:
                        continue
                    require(whole[first].read_bytes() == whole[second].read_bytes()
                            and prefix[first].read_bytes() == prefix[second].read_bytes(),
                            'same-batch mode checkpoints differ before cross-mode resume')
                    cross = root / f'{threads}-{first}-to-{second}.ck'
                    cross.write_bytes(prefix[first].read_bytes())
                    success(second, threads, 1, cross, 'same-batch-cross-mode-resume', resume_at=STEPS)
                    require(cross.read_bytes() == whole[second].read_bytes(),
                            'same-batch cross-mode resume differs from uninterrupted reference')
                    result['comparisons'].append(dict(kind='same-batch cross-mode resume',
                        sourceMode=first, targetMode=second, threads=threads,
                        batch=BATCH_BY_MODE[first], finalSha256=sha(cross)))

        # Match total payload size across batches while retaining distinct v2 headers.
        small = {}
        for name in MODES:
            threads = SMALL_WALKS // BATCH_BY_MODE[name]
            require(threads * BATCH_BY_MODE[name] == SMALL_WALKS, 'small population is not divisible')
            if (threads, name) in prefixes:
                small[name] = prefixes[threads, name]
            else:
                small[name] = root / f'small-{name}.ck'
                success(name, threads, 1, small[name], 'cross-batch-prefix')
        for first in MODES:
            for second in MODES:
                if BATCH_BY_MODE[first] != BATCH_BY_MODE[second]:
                    rejection(first, second, small[first], SMALL_WALKS//BATCH_BY_MODE[first],
                              SMALL_WALKS//BATCH_BY_MODE[second], 'cross-batch-rejection')
        for name in MODES:
            rejection(name, name, prefixes[8, name], 8, 9, 'worker-geometry-rejection')

    result['binarySha256After'] = {name:sha('/root/'+name) for name in MODES}
    require(result['binarySha256After'] == result['binarySha256'], 'checkpoint binaries changed')
    expected = {'uninterrupted':6, 'prefix':6, 'same-mode-resume':6,
                'same-batch-cross-mode-resume':6, 'cross-batch-prefix':2,
                'cross-batch-rejection':0, 'worker-geometry-rejection':2}
    require({kind:sum(row['kind'] == kind for row in result['rows']) for kind in expected} == expected
            and len(result['rows']) == sum(expected.values()), 'checkpoint gate coverage differs')
    expected_comparisons = {'same-mode split/resume':6, 'same-batch cross-mode resume':6}
    require(len(result['comparisons']) == 12
            and {kind:sum(row['kind'] == kind for row in result['comparisons'])
                 for kind in expected_comparisons} == expected_comparisons,
            'checkpoint comparison coverage differs')
    result['coverage'] = expected
    result['limits'] = ['No cross-batch checkpoint import or byte equality is permitted.',
        'Cross-batch walk-state equivalence is checked separately on normalized logical states.',
        'This gate exercises no-report resume and geometry rejection; separate client/state gates cover reports, guards and replay.']
    result['valid'] = True
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = dict(valid=False)
    try:
        checkpoint_gate(result)
    except Exception as exc:
        result['error'] = str(exc)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    if not result['valid']:
        raise RuntimeError(result.get('error', 'incomplete checkpoint gate'))


if __name__ == '__main__':
    main()
