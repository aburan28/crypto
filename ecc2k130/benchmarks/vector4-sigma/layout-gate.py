"""Compare the same deterministic scalar walks across polynomial build layouts.

Usage: python3 test_polytune_states.py base:32:/root/base b16:16:/root/b16
By default, four bounded runs of 64 slots per binary check checkpoint state,
DP/restart replay across resume, corpus equality, and overdue guard behavior.
Use --total-slots 8192 --state-only for a larger no-report state comparison.
These small-worker runs are correctness gates, not throughput measurements.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import tempfile

TOTAL = 64
HEADER = struct.Struct('<8s6IQ')
FINISHED = re.compile(r'finished: (\S+) M it/s, (\d+) distinguished points '
                      r'\((\d+) verified against the reference, (\d+) dropped\)')


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def config(value):
    parts = value.split(':', 2)
    require(len(parts) == 3, 'Expected name:batch:/path/to/binary')
    name, batch, binary = parts
    require(re.fullmatch(r'[A-Za-z0-9_-]+', name), 'Use a simple unique configuration name')
    batch = int(batch)
    require(batch > 0 and TOTAL % batch == 0, f'Batch must be a positive divisor of total slots {TOTAL}')
    return dict(name=name, batch=batch, threads=TOTAL // batch, binary=str(Path(binary).resolve()))


def initial_seed(run_id, index):
    # Independently reflects the public seed layout: 16 run bits, 32 walk
    # index bits, then 16 restart bits. Packed init uses id=slot*threads+tid.
    return (run_id << 48) | (index << 16)


def normalized_checkpoint(data, cfg, run_id, iteration):
    require(len(data) >= HEADER.size, 'Truncated checkpoint header')
    expected = (b'ECC2K130', 2, 131, cfg['threads'], cfg['batch'], 1, run_id, iteration)
    require(HEADER.unpack_from(data) == expected, 'Checkpoint header differs from requested settings')
    require(len(data) == HEADER.size + TOTAL * 60, 'Checkpoint payload length is not exact')
    points = []
    for index in range(TOTAL):
        slot, tid = divmod(index, cfg['threads'])
        coordinates = []
        for coordinate in range(2):
            words = tuple(struct.unpack_from('<I', data, HEADER.size + coordinate * TOTAL * 20 +
                                            4 * ((slot * 5 + word) * cfg['threads'] + tid))[0]
                          for word in range(5))
            require(words[4] & ~7 == 0, 'Checkpoint coordinate contains non-field bits')
            coordinates.append(words)
        dead = struct.unpack_from('<I', data, HEADER.size + TOTAL * 40 + 4 * index)[0]
        seed = struct.unpack_from('<Q', data, HEADER.size + TOTAL * 44 + 8 * index)[0]
        start = struct.unpack_from('<Q', data, HEADER.size + TOTAL * 52 + 8 * index)[0]
        require(dead in (0, 1), 'Packed dead value is not boolean')
        require(seed & ~0xffff == initial_seed(run_id, index), 'Checkpoint walk index/seed mapping changed')
        require(start <= iteration, 'Checkpoint start iteration exceeds current iteration')
        points.append((index, seed, dead, start, coordinates[0], coordinates[1]))
    # Index order is intentional: sorting could hide a wrong seed-to-slot map.
    return (run_id, iteration, tuple(points))


def digest_state(state):
    return hashlib.sha256(json.dumps(state, separators=(',', ':')).encode()).hexdigest()


def validate_output(output, cfg, steps, launches):
    require('MISMATCH' not in output and 'stopping:' not in output, output)
    require('warning: could not write checkpoint' not in output, output)
    require(re.findall(r'^packed polynomial state: ([01])$', output, re.MULTILINE) == ['1'],
            'Polynomial storage identity missing or wrong\n' + output)
    backend = re.findall(r'backend cuda-packed131: (\d+) threads x (\d+) slots x 1 lanes = (\d+) walks,', output)
    require(backend == [(str(cfg['threads']), str(cfg['batch']), str(TOTAL))], 'Scalar worker identity changed\n' + output)
    counts = re.findall(r'M it/s\s+(\d+) iterations', output)
    expected = TOTAL * steps * launches
    require(counts and int(counts[-1]) == expected, 'Incorrect newly executed scalar count\n' + output)
    final = FINISHED.findall(output)
    require(len(final) == 1, 'Missing or duplicate final summary\n' + output)
    rate, reports, verified, dropped = final[0]
    require(math.isfinite(float(rate)) and float(rate) > 0, 'Invalid completed-run rate')
    require(int(dropped) == 0, 'Reports were dropped')
    return dict(expectedIterations=expected, reportedIterations=int(counts[-1]),
                reports=int(reports), verified=int(verified), dropped=int(dropped))


def execute(cfg, arguments, steps, launches, run_id, timeout):
    command = [cfg['binary'], '--packed', '--curve', '131', '--run-id', str(run_id),
               '--threads', str(cfg['threads']), '--steps', str(steps), '--launches', str(launches)] + arguments
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    output = result.stdout + result.stderr
    require(result.returncode == 0, f'Failed command {command}\n{output}')
    row = dict(command=command, returncode=result.returncode, output=output)
    row.update(validate_output(output, cfg, steps, launches))
    return row


def corpus_records(path, run_id):
    data = path.read_bytes()
    require(len(data) % 32 == 0, 'Partial corpus record')
    result = Counter(data[offset:offset + 32] for offset in range(0, len(data), 32))
    for record in result:
        seed = struct.unpack_from('<Q', record)[0]
        index = (seed >> 16) & 0xffffffff
        require(index < TOTAL and seed >> 48 == run_id, 'Corpus contains an unexpected seed')
    return result


def check_config(cfg, directory, run_id, timeout, state_only=False):
    rows, states = [], {}
    cp = directory / 'whole.ck'
    row = execute(cfg, ['--bench', '--verify', '0', '--checkpoint', str(cp)], 16, 4, run_id, timeout)
    require((row['reports'], row['verified']) == (0, 0), 'No-report control produced a report')
    state = normalized_checkpoint(cp.read_bytes(), cfg, run_id, 64)
    for index, seed, dead, start, _, _ in state[2]:
        require(seed == initial_seed(run_id, index) and dead == 0 and start == 0,
                'No-report seed/dead/start state changed')
    states['whole'] = state
    row.update(stage='whole', normalizedSha256=digest_state(state))
    rows.append(row)
    if state_only:
        return (dict(config=cfg, runs=rows, stateOnly=True,
                     stateHashes={'whole': digest_state(state)}, reports=0,
                     sortedCorpusSha256=None), states, Counter())

    cp, corpus = directory / 'reports.ck', directory / 'reports.bin'
    total_reports = 0
    # A high cutoff yields short reference replays and repeated reports from
    # the same seed families. The second invocation restores the first.
    for part in (1, 2):
        row = execute(cfg, ['--dp-weight', '60', '--verify', str(2 * TOTAL), '--dp-cap', str(TOTAL),
                            '--checkpoint', str(cp), '--dp-file', str(corpus)], 32, 2, run_id, timeout)
        require(row['reports'] > 0 and row['verified'] == row['reports'],
                'Every produced DP must be replayed against the reference')
        if part == 2:
            require(re.search(r'resumed from .* at iteration 64\n', row['output']),
                    'Second report run did not restore the checkpoint')
        total_reports += row['reports']
        state = normalized_checkpoint(cp.read_bytes(), cfg, run_id, 64 * part)
        records = corpus_records(corpus, run_id)
        require(sum(records.values()) == total_reports, 'Corpus length differs from report counts')
        per_index = [[] for _ in range(TOTAL)]
        for record, multiplicity in records.items():
            seed = struct.unpack_from('<Q', record)[0]
            per_index[(seed >> 16) & 0xffffffff].extend([seed & 0xffff] * multiplicity)
        for index, seed, dead, start, _, _ in state[2]:
            restarts = sorted(per_index[index])
            require(restarts == list(range(len(restarts))), 'Missing or duplicate report/restart sequence')
            require(seed == initial_seed(run_id, index) + len(restarts) and dead == 0,
                    'Reports disagree with checkpoint reseed state')
            require(start % 32 == 0 and (bool(start) == bool(restarts)), 'Invalid report restart iteration')
        states[f'reports{part}'] = state
        row.update(stage=f'reports{part}', normalizedSha256=digest_state(state))
        rows.append(row)
    require(any(struct.unpack_from('<Q', record)[0] & 0xffff for record in records),
            'No restarted seed was independently replayed')

    guard = directory / 'guard.ck'
    row = execute(cfg, ['--dp-weight', '0', '--max-iters', '1', '--verify', '16',
                        '--checkpoint', str(guard)], 4096, 2, run_id, timeout)
    require((row['reports'], row['verified']) == (0, 0), 'Guard emitted false distinguished points')
    state = normalized_checkpoint(guard.read_bytes(), cfg, run_id, 8192)
    for index, seed, dead, start, _, _ in state[2]:
        require(seed == initial_seed(run_id, index) + 1 and dead == 0 and start == 8192,
                'Guard did not restart every overdue walk exactly once')
    states['guard'] = state
    row.update(stage='guard', normalizedSha256=digest_state(state))
    rows.append(row)
    return dict(config=cfg, runs=rows, stateOnly=False, stateHashes={k: digest_state(v) for k, v in states.items()},
                reports=total_reports, sortedCorpusSha256=hashlib.sha256(
                    b''.join(sorted(record for record, count in records.items() for _ in range(count)))).hexdigest()), states, records


def self_test(total_slots=64):
    global TOTAL
    TOTAL = total_slots
    # Serialize synthetic per-walk fields independently through explicit
    # slot/word/thread loops, then recover them with the production parser.
    run_id, iteration = 7, 64
    expected = tuple((i, initial_seed(run_id, i) + i % 4, i % 2, (i % 4) * 16,
                      tuple((i * 37 + word * 19 + 1) & (7 if word == 4 else 0xffffffff) for word in range(5)),
                      tuple((i * 53 + word * 23 + 2) & (7 if word == 4 else 0xffffffff) for word in range(5)))
                     for i in range(TOTAL))
    for batch in (16, 32, 64):
        cfg = dict(batch=batch, threads=TOTAL // batch)
        data = bytearray(HEADER.pack(b'ECC2K130', 2, 131, cfg['threads'], batch, 1, run_id, iteration))
        for coordinate in (4, 5):
            for slot in range(batch):
                for word in range(5):
                    for tid in range(cfg['threads']):
                        data.extend(struct.pack('<I', expected[slot * cfg['threads'] + tid][coordinate][word]))
        for field, fmt in ((2, '<I'), (1, '<Q'), (3, '<Q')):
            for point in expected:
                data.extend(struct.pack(fmt, point[field]))
        require(normalized_checkpoint(data, cfg, run_id, iteration) == (run_id, iteration, expected),
                'Synthetic SoA normalization failed')
        bad_seed = bytearray(data)
        struct.pack_into('<Q', bad_seed, HEADER.size + TOTAL * 44, initial_seed(run_id, 1))
        bad_top = bytearray(data)
        struct.pack_into('<I', bad_top, HEADER.size + 4 * 4 * cfg['threads'], 8)
        for bad in (data[:-1], data + b'\0', bad_seed, bad_top):
            try:
                normalized_checkpoint(bad, cfg, run_id, iteration)
            except AssertionError:
                pass
            else:
                raise AssertionError('Malformed synthetic checkpoint was accepted')
        iterations = TOTAL * 16 * 4
        output = (f"packed polynomial state: 1\nbackend cuda-packed131: {cfg['threads']} threads x {batch} slots x 1 lanes = {TOTAL} walks,\n"
                  f"0.1 s 1.000 M it/s {iterations} iterations 0 dp 0 stored 0 dropped\n"
                  "finished: 1.000 M it/s, 0 distinguished points (0 verified against the reference, 0 dropped)\n")
        require(validate_output(output, cfg, 16, 4)['reportedIterations'] == iterations, 'Synthetic scalar counter failed')
        for bad in (output.replace(f'{iterations} iterations', f'{iterations * 32} iterations'),
                    output.replace('polynomial state: 1', 'polynomial state: 0'),
                    output + output, output.replace('finished: 1.000', 'finished: nan')):
            try:
                validate_output(bad, cfg, 16, 4)
            except AssertionError:
                pass
            else:
                raise AssertionError('Malformed synthetic result was accepted')
    print(f'PASS: host {TOTAL}-slot B16/32/64 normalization and counters; malformed payload, seed, field and output controls rejected', flush=True)


def main():
    global TOTAL
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configs', nargs='*', help='name:batch:/path/to/binary')
    parser.add_argument('--run-id', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=240, help='Per-client timeout in seconds')
    parser.add_argument('--total-slots', type=int, default=64, help='Same scalar slots in every layout; must be divisible by every batch')
    parser.add_argument('--state-only', action='store_true', help='Only compare no-report whole state; skip DP/replay/resume/guard phases')
    parser.add_argument('--output', type=Path, help='Optional complete JSON correctness report')
    parser.add_argument('--self-test', action='store_true', help='Run host-only parser controls and exit')
    options = parser.parse_args()
    if options.self_test:
        for slots in (64, 8192):
            self_test(slots)
        return
    require(0 <= options.run_id <= 0xffff, 'Run id must fit the seed layout\'s 16 bits')
    require(options.timeout > 0, 'Timeout must be positive')
    require(0 < options.total_slots < (1 << 31), 'Total slots must be positive and fit the client worker counter')
    if not options.state_only:
        require(2 * options.total_slots < (1 << 31), 'Replay count must fit the client verification counter')
    TOTAL = options.total_slots
    configs = [config(value) for value in options.configs]
    require(len(configs) >= 2, 'At least two configurations are required for a comparison')
    require(len({c['name'] for c in configs}) == len(configs), 'Configuration names must be unique')
    require(all(Path(c['binary']).is_file() for c in configs), 'Every compiled binary must exist')
    report = dict(valid=False, kind='same-seed polynomial layout correctness', runId=options.run_id,
                  totalScalarSlots=TOTAL, stateOnly=options.state_only, results=[], binarySha256={c['name']: hashlib.sha256(
                      Path(c['binary']).read_bytes()).hexdigest() for c in configs})
    try:
        with tempfile.TemporaryDirectory(prefix='polytune-states-') as temporary:
            reference_states, reference_records = None, None
            for cfg in configs:
                print('CHECK ' + cfg['name'], flush=True)
                directory = Path(temporary) / cfg['name']
                directory.mkdir()
                try:
                    row, states, records = check_config(cfg, directory, options.run_id, options.timeout, options.state_only)
                except Exception as exc:
                    raise RuntimeError(cfg['name'] + ': ' + str(exc)) from exc
                report['results'].append(row)
                if reference_states is None:
                    reference_states, reference_records = states, records
                else:
                    require(states == reference_states, cfg['name'] + ': normalized checkpoint states differ')
                    require(records == reference_records, cfg['name'] + ': full report multisets differ')
                checked = ' exact scalar counts and normalized whole state' if options.state_only else ' exact scalar counts, normalized states and report/restart/guard checks'
                print('PASS: ' + cfg['name'] + checked, flush=True)
        report['valid'] = True
    except Exception as exc:
        report['error'] = str(exc)
    if options.output:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        options.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)
    require(report['valid'], report.get('error', 'Layout comparison failed'))


if __name__ == '__main__':
    main()
