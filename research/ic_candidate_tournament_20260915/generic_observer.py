"""Frozen native whole-mode controls; legacy observations never acquire admission."""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import sys
import time

from generic_admission import admit_rho
from generic_bases import verify_base
from generic_build import verify_binding, verify_build_record
from generic_driver import native_timing
from generic_phases import verify_native
from generic_queries import verify_queries
from generic_query_law import verify_query_law
from generic_stages import effective_config, verify_stages
from identity import natural, sha256
from oracle import require, verify
from tournament import digest, executed_job, execute, frozen_inputs, objhash, read, write

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / 'goal_20260924/generic-reference-qualification/PROTOCOL.md'


def entries(fixtures, arms):
    rng = random.Random(2026092664)
    cases = copy.deepcopy(fixtures['development'])
    rng.shuffle(cases)
    result = []
    for case in cases:
        initial = {a['id']: rng.randrange(2) for a in arms}
        for repetition in range(3):
            order = list(arms)
            rng.shuffle(order)
            for arm in order:
                modes = ['enabled', 'legacy']
                if (initial[arm['id']] + repetition) % 2:
                    modes.reverse()
                result.append(dict(case=case['id'], cell=case['cell'], arm=arm['id'],
                    repetition=repetition, order=modes,
                    name=f"{case['id']}-{arm['id']}-r{repetition}"))
    return result


def number_base(parent):
    # Each parent reserves a 2**16 block. Its observer uses the upper half.
    maximum = 6 * len(parent['run_aliases']) * parent['repetitions'] * parent.get('execution_number_stride', 1)
    require(maximum < 1 << 15, 'parent execution range overlaps observer reservation')
    return parent['run_number_base'] + (1 << 15)


def prepare(campaign, out, control):
    c, fixtures, ic = frozen_inputs(campaign)
    require(c['purpose'] == 'reference-qualification', 'observer needs a sealed development panel')
    arms = [a for a in ic + c['reference_arms'] if a.get('adapter') == 'generic-v1']
    expected = (2026092662, 3, 5, 45) if control else (2026092663, 15, 8, 360)
    require((c['seed'], len(fixtures['development']), len(arms)) == expected[:3],
            'observer panel differs from registration')
    plan = entries(fixtures, arms)
    require(len(plan) == expected[3] and len({p['name'] for p in plan}) == len(plan), 'invalid pair schedule')
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(campaign / 'evaluator', out / 'evaluator', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(Path(__file__), out / 'evaluator/generic_observer.py')
    shutil.copy2(PROTOCOL, out / 'PROTOCOL.md')
    contract = dict(schema_version=1, campaign_relative=os.path.relpath(campaign, out),
        campaign_contract_sha256=digest(campaign / 'contract.json'), control=control,
        host=platform.uname()._asdict(), resources=c['resources'], limits=c['limits'],
        arms=arms, entries=plan, seed=2026092664, bootstrap_seed=2026092665,
        run_number_base=number_base(c),
        run_number_policy='upper-half-of-parent-campaign-block-one-registered-observer-v1',
        maximum_pairs=50 if control else 400, expected_pairs=expected[3],
        measurement_scope='native whole-mode effects; legacy lacks scientific matrix/dispatch admission',
        promotion_eligible=False,
        evaluator_sha256={p.relative_to(out / 'evaluator').as_posix(): digest(p)
                          for p in (out / 'evaluator').rglob('*') if p.is_file()},
        protocol_sha256=digest(out / 'PROTOCOL.md'))
    write(out / 'contract.json', contract, exclusive=True)
    # Operational limits contain the parent's JSON wall-time float. Candidate
    # identities still use identity.sha256's stricter no-floats encoding.
    write(out / 'seal.json', dict(contract_sha256=objhash(contract)), exclusive=True)


def inputs(root):
    c = read(root / 'contract.json')
    require(objhash(c) == read(root / 'seal.json')['contract_sha256'], 'changed observer contract')
    require(digest(root / 'PROTOCOL.md') == c['protocol_sha256'], 'changed observer protocol')
    for name, expected in c['evaluator_sha256'].items():
        require(digest(HERE / name) == expected, 'changed observer evaluator: ' + name)
    campaign = (root / c['campaign_relative']).resolve()
    require(digest(campaign / 'contract.json') == c['campaign_contract_sha256'], 'changed parent campaign')
    parent, fixtures, _ = frozen_inputs(campaign)
    require(c['resources'] == parent['resources'] and c['limits'] == parent['limits'], 'observer resources differ')
    require(c['run_number_base'] == number_base(parent) and 2*c['maximum_pairs'] < 1 << 15,
            'observer execution range differs from reservation')
    require(entries(fixtures, c['arms']) == c['entries'] and len(c['entries']) == c['expected_pairs']
            and c['expected_pairs'] <= c['maximum_pairs'], 'changed observer schedule')
    for arm in c['arms']:
        manifest = campaign / arm['source_manifest_relative']
        source, build = read(manifest), read(manifest.parent / 'build-record.json')
        verify_build_record(build, source)
        require(digest(campaign / arm['binary_relative']) == build['worker_sha256'],
                'changed controlled observer executable')
    return c, campaign, {x['id']: x for x in fixtures['development']}, {a['id']: a for a in c['arms']}


def checked_observation(directory, *, mode, job, fixture, build, source, worker):
    process, report = read(directory / 'process.json'), read(directory / 'stdout.json')
    natural(process['process_wall_ns'], 'observer process nanoseconds', positive=True)
    require(type(process['exit_code']) is int, 'invalid process exit code')
    require(process['process_status'] == 'EXITED' and report['status'] in ('complete', 'incomplete')
            and process['exit_code'] == (0 if report['status'] == 'complete' else 2),
            'observer process failed, timed out or returned invalid output')
    require(report['fixture'] == fixture and report['mode'] == job['mode']
            and job['public_targets'] == fixture['targets'] and len(fixture['targets']) == 1,
            'changed observer workload')
    effective_config(job, report)
    require(report['target_input'] == 'supplied_public_point' and report['reusable_setup_excluded'] is True,
            'changed online boundary')
    signature = dict(status=report['status'], certificate=None, queries=None, base=None)
    if report['status'] == 'complete':
        signature['certificate'] = verify(report, fixture, expected_mode=job['mode'], summands=job['config']['summands'])
        require(report['scalar_replay_included'] is True, 'missing worker scalar replay')
    if job['mode'] == 'ic':
        signature['queries'] = verify_queries(report, fixture, job['config']['summands'])
        verify_query_law(report, fixture, job)
        signature['base'] = verify_base(report, fixture, job)
    else:
        signature['effective_walks'] = [s['effective_walks'] for s in report['solutions']]
    timing = audit = phases = None
    if mode == 'enabled':
        verify_binding(report, build, source, executable=worker)
        phases = verify_native(report, job, process_wall_ns=process['process_wall_ns'])
        audit = (verify_stages(report, fixture, job) if job['mode'] == 'ic' else
                 admit_rho(report, fixture, job, build, source, executable=worker,
                           process_wall_ns=process['process_wall_ns']))
        if report['status'] == 'complete':
            timing = native_timing(report, job, process['process_wall_ns'])
        outer = report['outer_online_wall_ns']
    else:
        require(job['exclusive_phases'] is False and report.get('generic_admission_schema') is None
                and report.get('generic_phase_timing') is None, 'legacy silently acquired exclusive admission')
        outer = report['online_wall_ns']
    if report['status'] == 'complete':
        natural(outer, 'observer outer nanoseconds', positive=True)
        require(outer <= process['process_wall_ns'], 'online interval exceeds process')
    else:
        require(outer is None and report['scalar_replay_included'] is False,
                'unstarted online work acquired a clock')
    return dict(signature=signature, outer_online_wall_ns=outer,
        process_wall_ns=process['process_wall_ns'], enabled_native_timing=timing,
        stage_audit=audit, phase_audit=phases, legacy_scientific_admission=None, instruction_cost=None,
        legacy_phase_costs=None, promotion_eligible=False)


def pair_record(root, c, campaign, case, arm, entry):
    directory = root / 'runs' / entry['name']
    manifest = campaign / arm['source_manifest_relative']
    source, build = read(manifest), read(manifest.parent / 'build-record.json')
    worker = campaign / arm['binary_relative']
    verify_build_record(build, source)
    require(digest(worker) == build['worker_sha256'], 'changed observer executable')
    record = dict(**entry, status='FAILED', reason=None, observations={}, promotion_eligible=False)
    try:
        for mode in entry['order']:
            job = executed_job(case, arm)
            job['exclusive_phases'] = mode == 'enabled'
            require(read(directory / mode / 'job.json') == job, 'changed observer job')
            process = read(directory / mode / 'process.json')
            require(process['memory_cap_bytes'] == c['limits']['memory_bytes']
                    and process['cpu'] == c['limits']['cpu'], 'changed observer process resources')
            require(type(process['process_wall_ns']) is int
                    and process['process_wall_seconds'] == process['process_wall_ns']/1e9,
                    'changed observer process clock units')
            record['observations'][mode] = checked_observation(directory / mode, mode=mode, job=job,
                fixture=case['fixture'], build=build, source=source, worker=worker)
        enabled, legacy = [record['observations'][m] for m in ('enabled', 'legacy')]
        require(sha256(enabled['signature']) == sha256(legacy['signature']), 'whole-mode semantics differ')
        record['status'] = 'MATCH_COMPLETE' if enabled['signature']['status'] == 'complete' else 'MATCH_INCOMPLETE'
    except (ValueError, KeyError, TypeError, OSError) as error:
        record['reason'] = f'{type(error).__name__}: {error}'
    admitted = read(campaign / 'admissions/development' / case['id'] / arm['id'] / 'admission.json')
    number = c['run_number_base'] + c['entries'].index(entry) * 2
    identity = admitted.get('candidate', admitted.get('reference'))
    key = 'candidate_id' if 'candidate' in admitted else 'reference_id'
    known = identity[key]
    workload = admitted['workload']['workload_id']
    legacy = 'OBS1h' + sha256(dict(build=build['build_sha256'], config=arm['config'], mode='legacy'))[:12]
    record.update(workload_id=workload, source_manifest_sha256=sha256(source), worker_sha256=build['worker_sha256'],
        run_ids=dict(enabled=f'{known}W{workload}R{number}', legacy=f'{legacy}W{workload}R{number+1}'),
        enabled_identity={key: known}, legacy_candidate_id=None,
        legacy_identity_scope='unadmitted observer diagnostic; OBS1 is not an IC candidate result')
    record['artifacts'] = {p.relative_to(directory).as_posix(): digest(p)
                           for p in directory.glob('*/*') if p.is_file()}
    return record


def effects(rows):
    if not rows or any(r['status'] != 'MATCH_COMPLETE' for r in rows):
        return None
    grouped = {}
    for row in rows:
        grouped.setdefault((row['cell'], row['case']), []).append(row)
    cells = {}
    for (cell, _), selected in sorted(grouped.items()):
        require(sorted(r['repetition'] for r in selected) == [0, 1, 2], 'incomplete observer repetitions')
        ratios = []
        for metric in ('outer_online_wall_ns', 'process_wall_ns'):
            a, b = [statistics.median(r['observations'][mode][metric] for r in selected)
                    for mode in ('enabled', 'legacy')]
            require(a > 0 and b > 0, 'nonpositive observer clock')
            ratios.append(math.log(a / b))
        cells.setdefault(cell, []).append(ratios)
    values = list(cells.values())
    estimate = [statistics.mean(statistics.mean(p[k] for p in cell) for cell in values) for k in (0, 1)]
    rng, draws = random.Random(2026092665), [[], []]
    for _ in range(10000):
        sampled = [[cell[rng.randrange(len(cell))] for _ in cell] for cell in values]
        for k in (0, 1):
            draws[k].append(statistics.mean(statistics.mean(p[k] for p in cell) for cell in sampled))
    result = {}
    for k, name in enumerate(('outer_online', 'process_wall')):
        draws[k].sort()
        result[name] = dict(enabled_over_legacy=math.exp(estimate[k]),
            descriptive_percentile_95=[math.exp(draws[k][i]) for i in (249, 9749)])
    return dict(metrics=result, target_count=len(grouped), cell_count=len(cells),
        scope='common outer interval whole-mode diagnostic; enabled scientific interval retained separately; no overhead subtraction, low-overhead test or performance promotion')


def summary(c, rows):
    require(len(rows) == c['expected_pairs'], 'partial observer schedule')
    keys = [key for r in rows for key in r['run_ids'].values()]
    require(len(keys) == len(set(keys)), 'observer identifiers collide')
    return dict(schema_version=1, pairs=len(rows), complete_pairs=sum(r['status'] == 'MATCH_COMPLETE' for r in rows),
        incomplete_pairs=sum(r['status'] == 'MATCH_INCOMPLETE' for r in rows),
        failed_pairs=sum(r['status'] == 'FAILED' for r in rows), promotion_eligible=False,
        all_observer_costs_charged=True, legacy_scientific_admission=None,
        effects={a['id']: effects([r for r in rows if r['arm'] == a['id']]) for a in c['arms']},
        interpretation='Semantic agreement supports fully charged instrumented-worker comparisons only; legacy remains diagnostic.')


def measured(root):
    c, campaign, cases, arms = inputs(root)
    require(platform.system() == 'Linux' and platform.machine() == 'x86_64'
            and platform.node() == c['host']['node'], 'observer must run on its declared Linux host')
    require(not (root / 'runs').exists(), 'partial observer runs are retained, never retried')
    rows = []
    for entry in c['entries']:
        case, arm = cases[entry['case']], arms[entry['arm']]
        directory = root / 'runs' / entry['name']
        directory.mkdir(parents=True, exist_ok=False)
        for mode in entry['order']:
            target = directory / mode
            target.mkdir()
            job = executed_job(case, arm)
            job['exclusive_phases'] = mode == 'enabled'
            write(target / 'job.json', job, exclusive=True)
            process = execute([str(campaign / arm['binary_relative'])], job, target,
                              c['limits']['timeout_seconds'], c['limits']['memory_bytes'], c['limits']['cpu'])
            write(target / 'process.json', process, exclusive=True)
        started = time.monotonic_ns()
        record = pair_record(root, c, campaign, case, arm, entry)
        write(directory / 'receipt.json', record, exclusive=True)
        write(directory / 'external-audit.json', dict(wall_ns=time.monotonic_ns()-started,
            boundary='external Python audit excluded from both worker intervals'), exclusive=True)
        rows.append(record)
        print(json.dumps(dict(pair=entry['name'], status=record['status'])), flush=True)
    result = summary(c, rows)
    write(root / 'summary.json', result, exclusive=True)
    require(result['failed_pairs'] == 0, 'observer mismatches/failures retained; qualification blocked')


def audit(root):
    c, campaign, cases, arms = inputs(root)
    rows = []
    for entry in c['entries']:
        expected = pair_record(root, c, campaign, cases[entry['case']], arms[entry['arm']], entry)
        require(read(root / 'runs' / entry['name'] / 'receipt.json') == expected, 'changed observer receipt')
        rows.append(expected)
    require(read(root / 'summary.json') == summary(c, rows), 'changed observer statistics')
    print(json.dumps(dict(status='VERIFIED', observer_pairs=len(rows), workers_executed=0)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'run', 'verify'))
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--round', type=Path)
    parser.add_argument('--control', action='store_true')
    args = parser.parse_args()
    if args.command == 'prepare':
        require(args.round is not None, 'supply the frozen parent campaign')
        prepare(args.round.resolve(), args.out.resolve(), args.control)
    else:
        (measured if args.command == 'run' else audit)(args.out.resolve())


if __name__ == '__main__':
    main()
