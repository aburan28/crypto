"""Reference selection from development evidence in the existing tournament.

No confirmation data, promotion decision, or new execution path lives here.
"""
import math
import statistics

from oracle import require


def geometric_mean(values):
    require(bool(values) and all(math.isfinite(v) and v > 0 for v in values),
            'reference statistic needs positive finite measurements')
    return math.exp(statistics.mean(math.log(v) for v in values))


def reference_report(root, contract, fixtures, ic_arms):
    from tournament import comparison, is_rho, load_stage, read, trial_path

    require(contract.get('purpose') == 'reference-qualification', 'not a qualification contract')
    require(read(root/'summaries/aa.json')['passed'], 'reference A/A control failed')
    cases = fixtures['development']
    arms = ic_arms + contract['reference_arms']
    rows = load_stage(root, 'development', cases, arms, contract['repetitions'])
    smoke = load_stage(root, 'smoke', fixtures['smoke'], arms, contract['repetitions'])
    comparisons = {arm['id']: comparison(rows, arm['id'], draws=contract['bootstrap_draws'],
                    match_support=contract['comparison_kind'] != 'factor-base-policy') for arm in arms}
    table = []
    for arm in arms:
        alias = arm['id']
        selected = [row for row in rows if row['arm'] == alias]
        smoke_failures = [row for row in smoke if row['arm']==alias and row['status']!='VERIFIED']
        record = dict(qualified=comparisons[alias].get('eligible',False) and not smoke_failures,
            smoke_failures=[{key:row.get(key) for key in ('case','repetition','status','reason')} for row in smoke_failures],
            alias=alias, mode='rho' if is_rho(arm) else 'ic',
            source_manifest_sha256=arm['source_manifest_sha256'], config=arm['config'],
            comparison_to_archived_incumbent=comparisons[alias],
            verified_runs=sum(row['status']=='VERIFIED' for row in selected),
            scheduled_runs=len(cases)*contract['repetitions'],
            failures=[{key:row.get(key) for key in ('case','repetition','status','reason')}
                      for row in selected if row['status']!='VERIFIED'], per_cell={}, effective_rho_widths={})
        if comparisons[alias].get('eligible'):
            per_cell = {}
            widths = {}
            for case in cases:
                case_rows = [row for row in selected if row['case'] == case['id']]
                require(len(case_rows) == contract['repetitions'], 'incomplete reference repetitions')
                metrics = dict(
                    instructions=statistics.median(row['total_operations'] for row in case_rows),
                    cold_ns=statistics.median(row['measurement']['native_timing']['cold']['wall_ns'] for row in case_rows),
                    online_ns=statistics.median(row['measurement']['native_timing']['online']['wall_ns'] for row in case_rows))
                metrics['normalized_S'] = metrics['instructions']/math.sqrt(int(case['fixture']['subgroup_order']))
                per_cell.setdefault(case['cell'], []).append(metrics)
                if is_rho(arm):
                    for row in case_rows:
                        directory = trial_path(root,'development',case,alias,row['repetition'])
                        native = read(directory/'native/stdout.json')
                        profiled = read(directory/'profile/stdout.json')
                        width = native['solutions'][0]['effective_walks']
                        require(type(width) is int and 1 <= width <= arm['config']['rho_parallel_walks']
                                and width == profiled['solutions'][0]['effective_walks'],
                                'invalid or changed executed rho width')
                        widths.setdefault(case['cell'], set()).add(width)
            record['per_cell'] = {cell:{key:geometric_mean([value[key] for value in values])
                for key in values[0]} for cell,values in per_cell.items()}
            record['effective_rho_widths'] = {cell:sorted(values) for cell,values in widths.items()}
        table.append(record)

    eligible_ic = [row for row in table if row['mode']=='ic' and row['qualified']]
    eligible_rho = [row for row in table if row['mode']=='rho' and row['qualified']]

    def select(candidates, online):
        if not candidates:
            return None
        def key(row):
            result = row['comparison_to_archived_incumbent']
            return ((result['online']['candidate_over_baseline'] if online else result['candidate_over_baseline']),
                    result['native_wall_candidate_over_baseline'], row['alias'])
        return min(candidates, key=key)['alias']

    return dict(schema_version=1, status='DEVELOPMENT_REFERENCES_SELECTED' if eligible_ic and eligible_rho else 'INCOMPLETE',
        classification='reference qualification on a frozen development panel; no improvement claim',
        promotion_eligible=False, improvement_rounds_used=0,
        selected_ic_cold=select(eligible_ic,False), selected_ic_online=select(eligible_ic,True),
        selected_rho_cold=select(eligible_rho,False), selected_rho_online=select(eligible_rho,True),
        selection_rule='Minimum paired equal-cell geometric mean ratio; native cold ratio then alias break exact ties. Keep cold and online leaders if different.',
        uncertainty='Development paired bootstrap intervals are descriptive; no confirmation or familywise promotion test was performed.',
        cases=len(cases), cells=contract['cells'], repetitions=contract['repetitions'],
        measured_stages=contract['stages'], heldout_data_used=False,
        table=table,
        next_gate='Freeze qualified sources, reference settings, cross-campaign target exclusions and familywise confirmation before an improvement round.')
