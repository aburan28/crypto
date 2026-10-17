#!/usr/bin/env python3
"""Bind two five-arm provider archives and select validated rows."""

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


prior = module('n83_k3072_prior_aggregate',
               ROOT / 'research/n83_e2e_final_holdout_20261016/aggregate.py')
runner = module('n83_k3072_aggregate_runner', HERE / 'run.py')
selector = module('n83_k3072_aggregate_select', HERE / 'select.py')


def read(path):
    return json.loads(path.read_text())


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'sealed output differs: {path}')
    path.write_text(blob)


def bind_half(half, results, archive, provider_sha, run_id, attempt,
              head_sha, finalist_validated, finalist_selection,
              reuse_local_replay):
    cloud_analysis, archive_sha = prior.finegrid.archive_analysis(
        archive, provider_sha)
    if (results / 'analysis.json').read_bytes() != cloud_analysis:
        raise ValueError(f'half {half}: extracted analysis differs from ZIP')
    member_count = prior.finegrid.bind_extraction(archive, results)
    log = results.parent
    stdout_path = log / 'local_replay.stdout.txt'
    stderr_path = log / 'local_replay.stderr.txt'
    status_path = log / 'local_replay.exit_status.txt'
    if reuse_local_replay:
        status = int(status_path.read_text().strip())
        stdout, stderr = stdout_path.read_text(), stderr_path.read_text()
    else:
        process = subprocess.run(
            [sys.executable, str(HERE / 'verify.py'), str(results),
             '--half', str(half), '--validated', str(finalist_validated),
             '--selection', str(finalist_selection)],
            cwd=ROOT, capture_output=True, text=True, check=False)
        status, stdout, stderr = (process.returncode, process.stdout,
                                  process.stderr)
        stdout_path.write_text(stdout)
        stderr_path.write_text(stderr)
        status_path.write_text(f'{status}\n')
    analysis = json.loads(cloud_analysis)
    if (status not in (0, 2) or stderr or json.loads(stdout) != analysis or
            (results / 'analysis.json').read_bytes() != cloud_analysis):
        raise ValueError(f'half {half}: independent local replay differs')
    if (status == 0) != (analysis['status'] == 'HALF_ADMITTED'):
        raise ValueError(f'half {half}: replay exit and admission differ')
    host, source = read(results / 'host.json'), read(results / 'source.json')
    if (analysis['schema'] != 'n83.smallk-challenge-half-analysis/v1' or
            analysis['half'] != half or
            analysis['point_indices'] != list(range(half * 5, half * 5 + 5)) or
            source['source_sha256'] != runner.SOURCE_SHA256 or
            source['targets_sha256'] != runner.TARGETS_SHA256 or
            source['answers_sha256'] != runner.ANSWERS_SHA256 or
            source['runner_sha256'] != sha256(HERE / 'run.py') or
            source['choice']['validated_sha256'] != sha256(finalist_validated) or
            source['choice']['selection_sha256'] != sha256(finalist_selection) or
            host['git_sha'] != head_sha or host['provider'] != 'github_actions' or
            host['parent_instance_id'] != f'gh-{run_id}-{attempt}-smallk-challenge-half' or
            host['instance_id'] !=
            f'gh-{run_id}-{attempt}-smallk-challenge-half-smallk-challenge-half{half}' or
            host['stratum_half'] != half or not host['runner_name']):
        raise ValueError(f'half {half}: host, source or point coverage differs')
    return {
        'half': half, 'analysis': analysis, 'source': source, 'host': host,
        'raw_artifact_sha256': archive_sha,
        'raw_artifact_zip': str(archive),
        'bound_result_members': member_count,
        'analysis_sha256': hashlib.sha256(cloud_analysis).hexdigest(),
        'local_replay_exit': status,
        'local_replay_stdout_sha256': sha256(stdout_path),
    }


def aggregate(halves, run_id, attempt, head_sha):
    first, second = halves
    keys = ('hashes', 'source_sha256', 'binary_sha256', 'rho_binary_sha256',
            'runner_sha256', 'targets_sha256', 'answers_sha256',
            'fixture_manifest_sha256', 'arms', 'choice', 'rho_seed_rule',
            'calibration_seed')
    same_source = all(first['source'][key] == second['source'][key] for key in keys)
    different_instances = (first['host']['runner_name'] != second['host']['runner_name'] and
                           first['host']['parent_instance_id'] !=
                           second['host']['parent_instance_id'])
    points = sorted(first['analysis']['points'] + second['analysis']['points'],
                    key=lambda item: item['point_index'])
    ids = [item[3] for item in first['source']['arms']]
    same_base = all(len({arm['raw_replay'].get('base_hash') for point in points
                         for arm in point['arms']
                         if arm['candidate_id'] == candidate and
                         arm['status'] == 'PASS_INDEPENDENT_REPLAY'}) == 1
                    for candidate in ids)
    host_strata = [
        {'half': item['half'],
         'point_indices': item['analysis']['point_indices'],
         'tier_status': item['analysis']['tier_status'],
         'aa_status': item['analysis']['aa_status'],
         'instance_id': item['host']['instance_id'],
         'runner_name': item['host']['runner_name'],
         'raw_artifact_sha256': item['raw_artifact_sha256'],
         'analysis_sha256': item['analysis_sha256']}
        for item in halves]
    evidence = ('PASS_INDEPENDENT_REPLAY' if same_source and same_base and
                different_instances and
                len(points) == 10 and
                {point['point_index'] for point in points} == set(range(10)) and
                all(item['analysis']['evidence_status'] ==
                    'PASS_INDEPENDENT_REPLAY' for item in halves) else
                'VERIFICATION_FAILED')
    validated = {
        'schema': 'n83.smallk-challenge-validated/v1',
        'source_sha256': runner.SOURCE_SHA256,
        'targets_sha256': runner.TARGETS_SHA256,
        'answers_sha256': runner.ANSWERS_SHA256,
        'candidate_ids': ids,
        'choice': first['source']['choice'],
        'evidence_status': evidence,
        'run_id': run_id, 'run_attempt': attempt,
        'head_sha': head_sha,
        'same_source_and_binary': same_source,
        'different_cloud_instances': different_instances,
        'same_candidate_bases_across_strata': same_base,
        'host_strata': host_strata,
        'points': points,
    }
    return validated, selector.select(validated)


def main():
    parser = argparse.ArgumentParser()
    for half in (0, 1):
        parser.add_argument(f'--half{half}-results', type=Path, required=True)
        parser.add_argument(f'--half{half}-zip', type=Path, required=True)
        parser.add_argument(f'--half{half}-provider-sha256', required=True)
    parser.add_argument('--run-id', type=int, required=True)
    parser.add_argument('--run-attempt', type=int, default=1)
    parser.add_argument('--head-sha', required=True)
    parser.add_argument('--finalist-validated', type=Path, required=True)
    parser.add_argument('--finalist-selection', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reuse-local-replay', action='store_true')
    args = parser.parse_args()
    halves = []
    for half in (0, 1):
        halves.append(bind_half(
            half, getattr(args, f'half{half}_results').resolve(),
            getattr(args, f'half{half}_zip').resolve(),
            getattr(args, f'half{half}_provider_sha256'),
            args.run_id, args.run_attempt, args.head_sha,
            args.finalist_validated.resolve(),
            args.finalist_selection.resolve(), args.reuse_local_replay))
    validated, result = aggregate(halves, args.run_id, args.run_attempt,
                                  args.head_sha)
    args.output.mkdir(parents=True, exist_ok=True)
    write_once(args.output / 'validated_rows.json', validated)
    write_once(args.output / 'selection.json', result)
    print(json.dumps({'evidence_status': validated['evidence_status'],
                      'online_selected': result['online']['selected_candidate_id'],
                      'cold_selected': result['cold']['selected_candidate_id']},
                     sort_keys=True))


if __name__ == '__main__':
    main()
