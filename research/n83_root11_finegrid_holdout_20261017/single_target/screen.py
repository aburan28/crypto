#!/usr/bin/env python3
"""Bounded diagnostic screen of intermediate N83 root11 bases on exposed p00."""

import hashlib
import importlib.util
import json
import subprocess
import time
from pathlib import Path

ROOT = Path('/Volumes/SSD990/crypto/worktrees/n83-final-online-holdout-20261017')
OUT = Path(__file__).resolve().parent
STUDY = ROOT / 'research/n83_smallk_challenge_20261017'
BIN = Path('/Volumes/SSD990/crypto/artifacts/n83_k_intermediate_screen_20261017/target/release/examples/koblitz_orbit_dlp_fast')
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
K_VALUES = (2752, 2880, 2944)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


dynamic = module('n83_screen_dynamic',
                 ROOT / 'research/n83_root11_dynamic_20261016/run.py')
group = module('n83_screen_group',
               ROOT / 'research/n83_cold_factor_base_20261016/verify.py').group


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'existing receipt differs: {path}')
    path.write_text(blob)


def main():
    if sha256(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256:
        raise ValueError('current source differs from the five-arm panel')
    if not BIN.is_file():
        raise ValueError('offline release binary has not finished building')
    target = OUT / 'p00_target.jsonl'
    target_line = (STUDY / 'targets.jsonl').read_text().splitlines()[0] + '\n'
    if target.exists() and target.read_text() != target_line:
        raise ValueError('exposed target differs')
    target.write_text(target_line)
    answer = json.loads((STUDY / 'answers.jsonl').read_text().splitlines()[0])
    frozen = {'source_sha256': SOURCE_SHA256,
              'binary_sha256': sha256(BIN),
              'target_sha256': sha256(target),
              'answer_sidecar_sha256': sha256(STUDY / 'answers.jsonl'),
              'k_values': K_VALUES,
              'workload_role': 'exposed_point_algorithmic_screen',
              'hardware_tier': 'ordinary_local_exploratory',
              'cap_seconds_per_arm': 900}
    write_once(OUT / 'inputs.json', frozen)
    for k in K_VALUES:
        folder = OUT / f'k{k}'
        folder.mkdir(exist_ok=True)
        command = [str(BIN), f'construct:83:1:{k}', str(target), '29',
                   str(folder / 'targets_result.jsonl')]
        env = dynamic.environment(folder, 'root11_dynamic')
        write_once(folder / 'command.json',
                   {'argv': command, 'environment': env, 'cap_seconds': 900})
        receipt = folder / 'screen_result.json'
        if receipt.exists():
            print(f'K{k}: existing receipt; not rerunning', flush=True)
            continue
        print(f'K{k}: starting', flush=True)
        start = time.monotonic()
        try:
            process = subprocess.run(command, cwd=ROOT,
                                     env=env,
                                     capture_output=True, text=True,
                                     timeout=900, check=False)
            status = 'PROCESS_EXIT_0' if process.returncode == 0 else 'PROCESS_FAILED'
            stdout, stderr = process.stdout, process.stderr
            exit_code = process.returncode
        except subprocess.TimeoutExpired as error:
            status = 'TIMEOUT'
            stdout = (error.stdout or b'').decode(errors='replace')
            stderr = (error.stderr or b'').decode(errors='replace')
            exit_code = None
        elapsed = time.monotonic() - start
        (folder / 'stdout.jsonl').write_text(stdout)
        (folder / 'stderr.txt').write_text(stderr)
        result = {'k': k, 'producer_status': status, 'exit_code': exit_code,
                  'outer_wall_seconds': elapsed, 'correctness': 'UNVERIFIED'}
        if status == 'PROCESS_EXIT_0':
            summary = json.loads(stdout)
            answers = [json.loads(line) for line in
                       (folder / 'targets_result.jsonl').read_text().splitlines()]
            if len(answers) != 1 or summary['targets_solved'] != 1:
                raise ValueError(f'K{k}: scalar count differs')
            point = group.point_from_json(answers[0]['target'])
            recovered = int(answers[0]['recovered_scalar'])
            expected = int(answer['scalar'])
            base = group.replay_base(
                json.loads((folder / 'base_dump.jsonl').read_text()), summary)
            indices = answers[0]['point_indices']
            relation = None
            for index in indices:
                relation = group.add(relation, base[index])
            if (recovered != expected or
                    group.multiply(group.G, recovered) != point or
                    relation != point or len(base) != 166 * k):
                raise ValueError(f'K{k}: independent arithmetic replay differs')
            result.update({'correctness': 'PASS_INDEPENDENT_SCALAR_BASE_RELATION',
                           'online_ms': answers[0]['target_ms'],
                           'target_probes': answers[0]['probes'],
                           'cold_process_ms': summary['timing_ms']['process_total'],
                           'root_state_allocated_bytes':
                               summary['root_state_allocated_bytes'],
                           'base_hash': summary['base_hash']})
        write_once(receipt, result)
        print(f'K{k}: {result}', flush=True)


if __name__ == '__main__':
    main()
