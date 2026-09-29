"""Apply the registered F4/F5 and SAT complete-solve gate to audited output.

`qualification.json` selects any eligible IC/rho reference. Its overall status
does not imply that either new algorithm family qualified. Run this read-only
gate after the pinned tournament verifier and natural-yield auditor; it keeps
failures and censored queries visible and never turns them into zero cost.
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

from oracle import require
from tournament import read, write

CELLS = ('n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0')
FAMILIES = {
    'f4_f5': ('generic_f4_dense', 'generic_f4_sparse', 'generic_f5_dense',
              'generic_inherited_f4_dense'),
    'sat': ('generic_sat_xor_dense', 'generic_sat_cnf_dense'),
}
GENERIC = ('generic_pair_dense', 'generic_pair_sparse', *FAMILIES['f4_f5'],
           *FAMILIES['sat'])
ALIASES = ('incumbent', *GENERIC, 'ic_online', 'rho', 'rho_online')
STAGES = {'smoke': 1, 'development': 3}  # distinct public points per cell
REPETITIONS = 3


def evaluate(summary, qualification, natural):
    require(summary['trial_slots'] == 750
            and summary['verified_native_profile_pairs'] + summary['retained_failures'] == 750
            and summary['distinct_record_ids'] == 1230
            and summary['exposed_points'] == 25
            and summary['status'] == ('AUDITED_COMPLETE' if summary['retained_failures'] == 0
                                      else 'AUDITED_WITH_FAILURES')
            and summary['promotion_eligible'] is False,
            'not the complete registered schedule')
    require(qualification['schema_version'] == 1
            and qualification['cells'] == list(CELLS)
            and qualification['cases'] == 15
            and qualification['repetitions'] == REPETITIONS
            and qualification['measured_stages'] == ['aa', 'smoke', 'development']
            and qualification['heldout_data_used'] is False
            and qualification['promotion_eligible'] is False,
            'not the registered development qualification')
    table = {row['alias']: row for row in qualification['table']}
    require(len(qualification['table']) == len(ALIASES) and set(table) == set(ALIASES),
            'missing, repeated or extra scheduled arms')
    require(natural['schema_version'] == 1 and natural['status'] == 'AUDITED'
            and natural['planted_decompositions_used_for_yield'] is False,
            'natural-query auditor missing or changed')

    observations = defaultdict(list)
    seen = set()
    for item in natural['observations']:
        key = (item['stage'], item['arm'], item['cell'])
        require(item['stage'] in STAGES and item['arm'] in GENERIC
                and item['cell'] in CELLS and item['repetition'] in range(REPETITIONS),
                'unexpected natural-query observation')
        identity = item['stage'], item['case'], item['arm'], item['repetition']
        require(identity not in seen, 'duplicate natural-query observation')
        seen.add(identity)
        observations[key].append(item)
    require(len(seen) == 480, 'missing natural-query process observations')
    yield_rows = {(row['stage'], row['arm'], row['cell']): row for row in natural['rows']}
    require(len(natural['rows']) == len(STAGES)*len(GENERIC)*len(CELLS)
            and len(yield_rows) == len(natural['rows']),
            'missing, repeated or extra natural-yield cells')

    for stage, cases_per_cell in STAGES.items():
        for arm in GENERIC:
            for cell in CELLS:
                key = stage, arm, cell
                require(key in yield_rows, 'missing natural-yield cell')
                row, items = yield_rows[key], observations[key]
                expected = cases_per_cell * REPETITIONS
                cases = defaultdict(set)
                for item in items:
                    cases[item['case']].add(item['repetition'])
                require(len(items) == expected and len(cases) == cases_per_cell
                        and all(reps == set(range(REPETITIONS)) for reps in cases.values()),
                        'natural-yield observations lost a point or repetition')
                audited = [item for item in items if item['audited']]
                outcomes = Counter()
                for item in audited:
                    outcomes.update(item['outcome_mix'])
                require(row['scheduled_runs'] == expected
                        and row['audited_bounded_runs'] == len(audited)
                        and row['censored_runs'] == expected-len(audited)
                        and row['verified_complete_runs'] == sum(
                            item['execution_status'] == 'VERIFIED' for item in items)
                        and row['observed_ordinary_queries'] == sum(
                            item['query_count'] or 0 for item in items)
                        and row['verified_witness_queries'] == outcomes['witness']
                        and row['independently_checked_outcomes'] == dict(sorted(outcomes.items()))
                        and sum(outcomes.values()) == row['observed_ordinary_queries'],
                        'natural-yield cell disagrees with process observations')

    gates = {}
    for arm in GENERIC:
        q = table[arm]
        require(q['mode'] == 'ic' and q['scheduled_runs'] == 45
                and q['verified_runs'] + len(q['failures']) == 45
                and type(q['qualified']) is bool
                and type(q['source_manifest_sha256']) is str
                and len(q['source_manifest_sha256']) == 64
                and all(ch in '0123456789abcdef' for ch in q['source_manifest_sha256']),
                'changed development qualification accounting')
        stage_rows = [yield_rows[stage, arm, cell] for stage in STAGES for cell in CELLS]
        smoke_verified = sum(row['verified_complete_runs'] for row in stage_rows
                             if row['stage'] == 'smoke')
        development_verified = sum(row['verified_complete_runs'] for row in stage_rows
                                   if row['stage'] == 'development')
        require(smoke_verified + len(q['smoke_failures']) == 15
                and development_verified == q['verified_runs'],
                'qualification and natural-query receipts disagree')
        comparison_eligible = q['comparison_to_archived_incumbent'].get('eligible', False)
        complete = bool(q['qualified'] and comparison_eligible
                        and smoke_verified == 15 and development_verified == 45)
        if complete:
            require(all(row['audited_bounded_runs'] == row['scheduled_runs']
                        and row['censored_runs'] == 0
                        and row['actual_base_and_columns'] is not None
                        for row in stage_rows),
                    'complete IC arm has censored or missing natural-query audit')
        outcomes = Counter()
        inventories = defaultdict(set)
        for row in stage_rows:
            outcomes.update(row['independently_checked_outcomes'])
            if row['actual_base_and_columns'] is not None:
                inventories[row['cell']].add(tuple(row['actual_base_and_columns']))
        require(all(len(shapes) <= 1 for shapes in inventories.values()),
                'factor-base inventory changed between smoke and development')
        observed = sum(outcomes.values())
        reasons = []
        if q['smoke_failures']:
            reasons.append('smoke failure')
        if q['failures']:
            reasons.append('development failure')
        if not comparison_eligible:
            reasons.append('ineligible paired complete-cost comparison')
        if not complete and not reasons:
            reasons.append('qualification report did not admit complete arm')
        gates[arm] = dict(complete=complete, smoke_verified=smoke_verified,
                          smoke_scheduled=15, development_verified=development_verified,
                          development_scheduled=45, reasons=reasons,
                          source_manifest_sha256=q['source_manifest_sha256'],
                          ordinary_queries_observed_across_processes=observed,
                          ordinary_non_witness_queries_observed_across_processes=
                              observed-outcomes['witness'],
                          ordinary_outcome_mix=dict(sorted(outcomes.items())),
                          censored_processes=sum(row['censored_runs'] for row in stage_rows),
                          actual_base_and_columns_by_cell={cell:list(next(iter(shapes)))
                              for cell, shapes in sorted(inventories.items()) if shapes})
    qualified = {family:[arm for arm in aliases if gates[arm]['complete']]
                 for family, aliases in FAMILIES.items()}
    return dict(schema_version=1,
                status='F4_F5_AND_SAT_QUALIFIED' if all(qualified.values())
                    else 'NEGATIVE_FAMILY_QUALIFICATION',
                qualified_families=qualified, arms=gates,
                reference_selection_status=qualification['status'],
                scope='Two-family completeness gate over the registered one-target development '
                      'report; source and raw-query replay are required separately. Paired '
                      'incumbent/rho performance and any promotion also require separate audit.',
                process_repetitions_are_not_independent_points=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True,
                        help='Finished registered runner output containing summary and natural yield.')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    result = evaluate(read(bundle/'summary.json'),
                      read(bundle/'tournament/qualification.json'),
                      read(bundle/'natural-yield.json'))
    write(args.out, result, exclusive=True)
    print(result['status'])


if __name__ == '__main__':
    main()
