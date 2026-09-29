"""Independently audit natural-query outcomes in retained generic IC runs.

The qualification evaluator checks complete solves. This audit also checks
bounded *incomplete* reports before using their attempted-query histories.
Timeouts without a full report remain censored, never zero-yield observations.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import random

from generic_build import verify_binding
from generic_phases import verify_native
from generic_stages import verify_stages
from measurement import report_sha256
from oracle import require
from tournament import executed_job, frozen_inputs, read, stage_arms, trial_path, write

SEED = 2026092902
DRAWS = 10000


def observed_run(root, case, arm, stage, repetition):
    directory = trial_path(root, stage, case, arm['id'], repetition)
    receipt = read(directory/'receipt.json')
    raw = directory/'profile/stdout.json'
    base = dict(stage=stage, cell=case['cell'], case=case['id'], arm=arm['id'],
                repetition=repetition, execution_status=receipt['status'],
                profile_status=None, audited=False, query_count=None,
                witness_count=None, outcome_mix=None, accepted_rows=None,
                rank=None, base_points=None, folded_columns=None,
                report_sha256=None, censor_reason=None,
                profiled_phase_wall_ns=None, complete_instruction_phases=None)
    base['target_descent_queries'] = None
    base['target_outcome_mix'] = None
    if not raw.exists():
        base['censor_reason'] = 'no complete profiler report'
        return base
    try:
        report = read(raw)
    except (json.JSONDecodeError, UnicodeDecodeError):
        base['censor_reason'] = 'invalid or interrupted profiler JSON'
        return base
    base['profile_status'] = report.get('status')
    if report.get('status') not in {'complete', 'incomplete'}:
        base['censor_reason'] = 'profiler did not return a bounded IC report'
        return base
    if report['status'] == 'complete' and receipt['status'] != 'VERIFIED':
        base['censor_reason'] = 'complete JSON from an unverified process'
        return base
    require(report['status'] == 'complete' or receipt['status'] != 'VERIFIED',
            'incomplete report contradicts verified tournament receipt')
    manifest = read(root/arm['source_manifest_relative'])
    build = read(root/arm['build_record_relative'])
    verify_binding(report, build, manifest, executable=root/arm['binary_relative'])
    job = executed_job(case, arm)
    require(read(directory/'job.json') == job, 'changed frozen job')
    audit = verify_stages(report, case['fixture'], job)
    phases = verify_native(report, job)
    attempts = [item for batch in report['collection_reports'] for item in batch['attempts']]
    outcomes = Counter(item['pdp']['outcome'] for item in attempts)
    descents = [item for solution in report.get('solutions', [])
                for item in solution['attempts']]
    descent_outcomes = Counter(item['pdp']['outcome'] for item in descents)
    require(len(attempts) == audit['query_law']['collection_queries']
            and len(attempts) == len(audit['matrix']['rank_trajectory'])
            and len(descents) == audit['query_law']['descent_queries'],
            'attempt history differs from independent query and matrix checks')
    base.update(audited=True, query_count=len(attempts),
                witness_count=outcomes['witness'], outcome_mix=dict(sorted(outcomes.items())),
                accepted_rows=audit['matrix']['accepted_rows'], rank=audit['matrix']['rank'],
                base_points=audit['base']['inventory']['usable_point_count'],
                folded_columns=audit['base']['inventory']['effective_columns'],
                report_sha256=report_sha256(report),
                target_descent_queries=len(descents),
                target_outcome_mix=dict(sorted(descent_outcomes.items())),
                profiled_phase_wall_ns=phases['observed_phases_ns'],
                complete_instruction_phases=receipt['phase_costs'] if
                    receipt['status'] == 'VERIFIED' else None)
    return base


def case_rate_interval(cases):
    """Resample distinct public points, not deterministic process repetitions."""
    zero_query_points = sum(row['query_count'] == 0 for row in cases)
    cases = [row for row in cases if row['query_count'] > 0]
    if not cases:
        return dict(rate=None, ci95=None, distinct_points=0,
                    zero_query_points=zero_query_points,
                    point_mean=None, point_mean_hoeffding95=None,
                    uncertainty='no ordinary query exposure; rate is unknown')
    total = sum(row['query_count'] for row in cases)
    estimate = sum(row['witness_count'] for row in cases)/total
    point_rates = [row['witness_count']/row['query_count'] for row in cases]
    point_mean = sum(point_rates)/len(point_rates)
    radius = math.sqrt(math.log(40)/(2*len(cases)))
    conservative = [max(0, point_mean-radius), min(1, point_mean+radius)]
    if len(cases) == 1:
        return dict(rate=estimate, ci95=None, distinct_points=1,
                    zero_query_points=zero_query_points,
                    point_mean=point_mean, point_mean_hoeffding95=conservative,
                    uncertainty='one distinct point; bootstrap interval uninformative')
    rng = random.Random(SEED)
    draws = []
    for _ in range(DRAWS):
        sample = rng.choices(cases, k=len(cases))
        draws.append(sum(row['witness_count'] for row in sample)
                     / sum(row['query_count'] for row in sample))
    draws.sort()
    return dict(rate=estimate, ci95=[draws[int(.025*DRAWS)], draws[int(.975*DRAWS)]],
                distinct_points=len(cases), zero_query_points=zero_query_points,
                point_mean=point_mean,
                point_mean_hoeffding95=conservative,
                uncertainty='query-weighted point-cluster bootstrap is descriptive; '
                            'Hoeffding bound is for the distinct-point mean under the frozen point law')


def summarize(observations, repetitions=3):
    require(repetitions in (1, 3), 'unknown process repetition schedule')
    groups = defaultdict(list)
    for row in observations:
        groups[(row['stage'], row['arm'], row['cell'])].append(row)
    result = []
    for (stage, arm, cell), runs in sorted(groups.items()):
        by_case = defaultdict(list)
        for row in runs:
            by_case[row['case']].append(row)
        stable_cases = []
        for case, case_runs in sorted(by_case.items()):
            require({row['repetition'] for row in case_runs} == set(range(repetitions))
                    and len(case_runs) == repetitions,
                    'lost or repeated process repetition: '+case)
            if all(row['audited'] for row in case_runs):
                signatures = {(row['query_count'], row['witness_count'],
                               tuple(sorted(row['outcome_mix'].items())),
                               row['accepted_rows'], row['rank']) for row in case_runs}
                require(len(signatures) == 1,
                        'deterministic query/rank history differs across process repetitions')
                stable_cases.append(case_runs[0])
        outcomes = Counter()
        descent_outcomes = Counter()
        for row in runs:
            if row['audited']:
                outcomes.update(row['outcome_mix'])
                descent_outcomes.update(row['target_outcome_mix'])
        inventories = sorted({(row['base_points'], row['folded_columns'])
                              for row in runs if row['audited']})
        require(len(inventories) <= 1, 'factor-base inventory changed within an arm/cell')
        complete_costs = [sum(row['complete_instruction_phases'][phase] for phase in
                              ('queries', 'pdp', 'relation_check')) for row in runs
                          if row['complete_instruction_phases'] is not None]
        complete_rows = sum(row['accepted_rows'] for row in runs
                            if row['complete_instruction_phases'] is not None)
        complete_rank = sum(row['rank'] for row in runs
                            if row['complete_instruction_phases'] is not None)
        charged = sum(complete_costs) if complete_costs else None
        result.append(dict(stage=stage, arm=arm, cell=cell, scheduled_runs=len(runs),
            verified_complete_runs=sum(row['execution_status'] == 'VERIFIED' for row in runs),
            audited_bounded_runs=sum(row['audited'] for row in runs),
            censored_runs=sum(not row['audited'] for row in runs),
            censor_reasons=dict(Counter(row['censor_reason'] for row in runs if not row['audited'])),
            observed_ordinary_queries=sum(row['query_count'] or 0 for row in runs),
            independently_checked_outcomes=dict(sorted(outcomes.items())),
            verified_witness_queries=outcomes['witness'],
            observed_target_descent_queries=sum(row['target_descent_queries'] or 0 for row in runs),
            target_descent_outcomes=dict(sorted(descent_outcomes.items())),
            accepted_rows=sum(row['accepted_rows'] or 0 for row in runs),
            final_ranks=sorted({row['rank'] for row in runs if row['audited']}),
            actual_base_and_columns=inventories[0] if inventories else None,
            complete_query_pdp_relation_ir=charged,
            complete_cost_per_accepted_row_ir=charged/complete_rows if complete_rows else None,
            complete_cost_per_rank_increment_ir=charged/complete_rank if complete_rank else None,
            natural_witness_rate=case_rate_interval(stable_cases),
            qualification='all scheduled runs complete' if all(
                row['execution_status'] == 'VERIFIED' for row in runs) else 'incomplete'))
    return result


def audit_campaign(root):
    root = Path(root)
    contract, fixtures, arms = frozen_inputs(root)
    require(contract['qualification_reference_schema'] == 1,
            'not the frozen generic backend qualification')
    observations = []
    for stage in ('smoke', 'development'):
        active = [arm for arm in stage_arms(root, stage, arms)
                  if arm.get('adapter') == 'generic-v1' and arm['id'].startswith('generic_')]
        for case in fixtures[stage]:
            for arm in active:
                for repetition in range(contract['repetitions']):
                    observations.append(observed_run(root, case, arm, stage, repetition))
    return dict(schema_version=1, status='AUDITED', observations=observations,
                rows=summarize(observations, contract['repetitions']), bootstrap_seed=SEED, bootstrap_draws=DRAWS,
                rate_scope='natural sampled ordinary queries; one independent point per case; '
                           'censored runs and audited zero-query points excluded from rate '
                           'denominator and counted separately',
                planted_decompositions_used_for_yield=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--round', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = audit_campaign(args.round)
    write(args.out, result, exclusive=True)
    print(json.dumps(dict(status=result['status'], rows=len(result['rows']),
                          observations=len(result['observations']))))


if __name__ == '__main__':
    main()
