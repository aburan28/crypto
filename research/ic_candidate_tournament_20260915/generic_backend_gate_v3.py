"""F4/F5-only smoke gate for the third registration.

SAT is out of scope. Development is not scheduled. Input must already have
passed the retained tournament verifier and natural-query auditor.
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

from oracle import require
from tournament import read, write

CELLS = ('n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0')
PANEL_SHA256 = 'df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4'
LOST_V1_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'
LOST_V2_SHA256 = '0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54'
F4_F5 = ('generic_f4_subspace_dense', 'generic_f5_subspace_dense')
GENERIC = ('generic_pair_subspace_dense', *F4_F5)
STAGES = {'smoke': 1}


def evaluate_v3(summary, natural):
    require(summary.get('registration_panel_sha256') == PANEL_SHA256
            and summary.get('prior_censored_exposures_sha256') == LOST_V1_SHA256
            and summary.get('prior_v2_censored_exposures_sha256') == LOST_V2_SHA256
            and summary.get('qualification_schedule') == 'smoke'
            and summary['trial_slots'] == 45
            and summary['repetitions'] == 1
            and summary['exposed_points'] == 10
            and summary['verified_native_profile_pairs'] + summary['retained_failures'] == 45
            and summary['status'] == ('AUDITED_COMPLETE' if summary['retained_failures'] == 0
                                      else 'AUDITED_WITH_FAILURES')
            and summary['promotion_eligible'] is False
            and summary.get('sat_in_scope') is False,
            'not the registered v3 smoke schedule')
    require(natural['schema_version'] == 1 and natural['status'] == 'AUDITED'
            and natural.get('qualification_schedule') == 'smoke'
            and natural['planted_decompositions_used_for_yield'] is False,
            'natural-query auditor missing or changed')

    observations = defaultdict(list)
    seen = set()
    for item in natural['observations']:
        key = (item['stage'], item['arm'], item['cell'])
        require(item['stage'] in STAGES and item['arm'] in GENERIC
                and item['cell'] in CELLS and item['repetition'] == 0,
                'unexpected natural-query observation')
        identity = item['stage'], item['case'], item['arm'], item['repetition']
        require(identity not in seen, 'duplicate natural-query observation')
        seen.add(identity)
        observations[key].append(item)
    require(len(seen) == len(STAGES) * len(GENERIC) * len(CELLS),
            'missing natural-query process observations')
    yield_rows = {(row['stage'], row['arm'], row['cell']): row for row in natural['rows']}
    require(len(natural['rows']) == len(STAGES) * len(GENERIC) * len(CELLS)
            and len(yield_rows) == len(natural['rows']),
            'missing, repeated or extra natural-yield cells')

    gates = {}
    for arm in GENERIC:
        stage_rows = [yield_rows['smoke', arm, cell] for cell in CELLS]
        smoke_verified = sum(row['verified_complete_runs'] for row in stage_rows)
        require(smoke_verified + sum(row['scheduled_runs'] - row['verified_complete_runs']
                                     for row in stage_rows) == 5,
                'natural-query smoke schedule incomplete for '+arm)
        complete = smoke_verified == 5 and all(
            row['audited_bounded_runs'] == row['scheduled_runs']
            and row['censored_runs'] == 0
            and row['actual_base_and_columns'] is not None
            for row in stage_rows)
        outcomes = Counter()
        inventories = {}
        for row in stage_rows:
            outcomes.update(row['independently_checked_outcomes'])
            if row['actual_base_and_columns'] is not None:
                inventories[row['cell']] = list(row['actual_base_and_columns'])
        reasons = []
        if smoke_verified < 5:
            reasons.append('smoke failure')
        if not complete and not reasons:
            reasons.append('censored or incomplete natural-query audit')
        gates[arm] = dict(complete=complete, smoke_verified=smoke_verified,
                          smoke_scheduled=5, reasons=reasons,
                          ordinary_queries_observed_across_processes=sum(outcomes.values()),
                          ordinary_non_witness_queries_observed_across_processes=
                              sum(outcomes.values()) - outcomes['witness'],
                          ordinary_outcome_mix=dict(sorted(outcomes.items())),
                          censored_processes=sum(row['censored_runs'] for row in stage_rows),
                          actual_base_and_columns_by_cell=inventories)
    qualified = [arm for arm in F4_F5 if gates[arm]['complete']]
    return dict(schema_version=1,
                status='F4_F5_SMOKE_QUALIFIED' if qualified else 'NEGATIVE_FAMILY_QUALIFICATION',
                qualified_f4_f5_arms=qualified, arms=gates,
                registration_panel_sha256=PANEL_SHA256,
                prior_censored_exposures_sha256=LOST_V1_SHA256,
                prior_v2_censored_exposures_sha256=LOST_V2_SHA256,
                scope='F4/F5-only smoke completeness gate; SAT out of scope; no development '
                      'stage, held-out confirmation, or familywise promotion.',
                sat_in_scope=False,
                repetitions=1,
                process_repetitions_are_not_independent_points=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    result = evaluate_v3(read(bundle/'summary.json'), read(bundle/'natural-yield.json'))
    write(args.out, result, exclusive=True)
    print(result['status'])


if __name__ == '__main__':
    main()
