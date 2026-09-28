"""Second/third-round protocol with separately bound cold and online references.

The first-round rules and archived evaluators stay unchanged. This version uses
the same three-attempt familywise budget, fixed cells and complete-cost gates.
"""
import copy
import math
import random
import statistics

import campaign_rules as legacy
from oracle import require

PURPOSE = 'bounded-improvement-20260924-v2'
RULE = copy.deepcopy(legacy.RULE)  # Same inference family across all three attempts.
CELLS, HOLDOUT_CELLS = legacy.CELLS, legacy.HOLDOUT_CELLS
CONFIG, HISTORY_SHA256 = legacy.CONFIG, legacy.HISTORY_SHA256
METRICS = legacy.METRICS
IC_SOURCE = legacy.IC_SOURCE
ONLINE_IC_SOURCE = legacy.COLD_RHO_SOURCE
GENERIC_SOURCE = '72da739cbc72e351e85c0184664ea60b6d2f23b54aab0e9ff45fdd1a60714d7d'
GENERIC_BUILD = '54316bcec91a68a5a08bcc8e92548a46ad75e2a27d8c056fe6520d79f3410962'
REFERENCE_ARCHIVE_SHA256 = '680369a9b822dcc3321380bff751bd2b94f0892c915f201a67916fc7bd49fc5a'
QUALIFICATION_SHA256 = '8479d9b9b777297d9c05dd860d5189ecc0791db9c70b6c641bc19ec986fae256'
OBSERVER_SHA256 = '96b687f7887442f9059d733225274dabd07b25900263eabb17f721636089dee6'
METRIC_REFERENCES = dict(instructions='incumbent', cold_ns='incumbent', online_ns='ic_online')
REFERENCE_ROLES = ('ic_online', 'rho', 'rho_online')
MAX_CANDIDATES = 11  # Incumbent plus ten challengers; three additional references.
validate_environment = legacy.validate_environment
selection_key = legacy.selection_key


def expected_binding():
    declarations = (
        ('incumbent', 'incumbent', IC_SOURCE, CONFIG, None, None),
        ('ic_online', 'prepared_both', ONLINE_IC_SOURCE, CONFIG, None, None),
        ('rho', 'rho_incumbent_8', IC_SOURCE, dict(CONFIG, rho_parallel_walks=8), None, None),
        ('rho_online', 'rho_generic_dense_16', GENERIC_SOURCE,
         dict(CONFIG, linear_algebra='dense', rho_parallel_walks=16), 'generic-v1', GENERIC_BUILD))
    return dict(schema_version=2, archive_sha256=REFERENCE_ARCHIVE_SHA256,
        qualification_sha256=QUALIFICATION_SHA256, observer_sha256=OBSERVER_SHA256,
        metric_references=copy.deepcopy(METRIC_REFERENCES),
        scope='fully charged instrumented pipelines; legacy observer remains diagnostic',
        bindings={name: dict(selected_alias=selected, source_manifest_sha256=source,
                            configuration=copy.deepcopy(config), adapter=adapter, build_sha256=build)
                  for name, selected, source, config, adapter, build in declarations})


def qualified_binding(report, observer, incumbent, references):
    require(legacy.object_hash(report) == QUALIFICATION_SHA256, 'unaccepted v2 qualification report')
    require(legacy.object_hash(observer) == OBSERVER_SHA256, 'unaccepted observer evidence')
    require(observer['complete_pairs'] == 360 and observer['failed_pairs'] == 0
            and observer['incomplete_pairs'] == 0 and observer['all_observer_costs_charged']
            and observer['legacy_scientific_admission'] is None,
            'observer evidence does not support fully charged comparisons')
    require([arm['id'] for arm in references] == list(REFERENCE_ROLES),
            'separate IC online and both rho references are required in canonical order')
    expected = expected_binding()
    table = {row['alias']: row for row in report['table']}
    selected = dict(incumbent=report['selected_ic_cold'], ic_online=report['selected_ic_online'],
                    rho=report['selected_rho_cold'], rho_online=report['selected_rho_online'])
    for arm in [incumbent] + references:
        name = arm['id']
        binding = expected['bindings'][name]
        row = table[selected[name]]
        require(row['qualified'] and row['verified_runs'] == row['scheduled_runs'] == 45
                and not row['failures'] and not row['smoke_failures'], 'incomplete bound reference')
        require(selected[name] == binding['selected_alias']
                and arm['source_manifest_sha256'] == row['source_manifest_sha256'] == binding['source_manifest_sha256']
                and arm['config'] == row['config'] == binding['configuration']
                and arm.get('adapter') == binding['adapter']
                and arm.get('build_sha256') == binding['build_sha256'],
                'changed qualified v2 reference ' + name)
        if name != 'incumbent':
            require(arm.get('kind') == ('ic-reference' if name == 'ic_online' else 'rho-reference'),
                    'changed reference role ' + name)
        else:
            require(arm.get('kind') is None, 'incumbent must remain an IC candidate')
    return expected


def scheduled_slot_bound(candidate_count):
    require(type(candidate_count) is int and 2 <= candidate_count <= MAX_CANDIDATES,
            'v2 permits at most ten challengers plus the incumbent')
    # 5 A/A points, 5 smoke, 15 development, 15 selection, 72 confirmation
    # and the same 72 points on replay. Every point gets three processes.
    return (30 + 60 * (candidate_count + 3) +
            45 * (1 + min(6, candidate_count - 1) + 3) + 432 * 5)


def validate_reference_declarations(references):
    require(isinstance(references, list) and [a['id'] for a in references] == list(REFERENCE_ROLES),
            'changed v2 reference registry')
    bindings = expected_binding()['bindings']
    for arm in references:
        name = arm['id']
        require(arm.get('kind') == ('ic-reference' if name == 'ic_online' else 'rho-reference')
                and arm.get('config') == bindings[name]['configuration']
                and arm.get('adapter') == bindings[name]['adapter'],
                'changed declared reference ' + name)


def validate_preparation(args):
    legacy.validate_preparation(args)
    require(args.attempt_number in (2, 3), 'v2 cannot reopen the sealed first round')
    require(getattr(args, 'reference_registry', None) and getattr(args, 'qualified_observer', None),
            'v2 needs the complete reference registry and accepted observer evidence')
    require(args.qualified_report and args.target_history,
            'v2 needs accepted qualification evidence and initial target history')
    require(not args.rho_source_root and not args.rho_config, 'v2 uses its explicit reference registry')


def validate_contract(contract):
    require(contract.get('purpose') == PURPOSE, 'not the versioned v2 protocol')
    require(contract.get('scientific_admission') is True and contract.get('schema_version') == 2
            and contract.get('unit') == 'valgrind-3.22-amd64-Ir'
            and contract.get('stages') == ['aa', 'smoke', 'development', 'selection', 'confirmation', 'replay'],
            'v2 requires complete scientific admission, phases and calibrated units')
    require(contract.get('familywise_rule') == RULE, 'changed original three-attempt inference budget')
    attempt = contract.get('attempt_number')
    require(type(attempt) is int and attempt in (2, 3)
            and contract.get('seed') == 2026092550 + attempt, 'invalid v2 round or seed')
    require(type(contract.get('target_exposure_schema')) is int and contract['target_exposure_schema'] == 1,
            'v2 requires sealed supplemental target exclusions')
    validate_environment(system=contract['host']['system'], machine=contract['host']['machine'],
                         compiler=contract['compiler'], profiler=contract['profiler_version'])
    require(contract['confirmation_cases'] == 72 and contract['repetitions'] == 3
            and contract['cells'] == CELLS and contract['holdout_cells'] == HOLDOUT_CELLS
            and contract['confirmation_cases_per_cell'] == {cell: 12 for cell in CELLS + HOLDOUT_CELLS},
            'changed fixed confirmation allocation')
    require(contract['confirmation_ratio'] == .8 and contract['max_cell_ratio'] == 1.1
            and contract['require_native_progress'] and contract['objective'] == 'incumbent',
            'changed complete-cost objective')
    require(contract['selection_width'] == 6 and contract['exploration_slots'] == 1,
            'changed diverse portfolio budget')
    limits = contract['limits']
    require(limits['max_profiled_jobs'] == 3500 and limits['memory_bytes'] == 8 * 1024**3
            and limits['timeout_seconds'] == 180 and limits['worker_threads'] == 1
            and contract['target_count'] == 1, 'changed resource or target limits')
    require(contract.get('reference_qualification') == expected_binding(),
            'missing or changed exact v2 reference binding')
    require(contract.get('metric_references') == METRIC_REFERENCES
            and contract.get('execution_number_stride') == 2, 'changed reference or execution policy')
    require(contract.get('scheduled_slot_bound') == scheduled_slot_bound(contract['candidate_count'])
            and contract['scheduled_slot_bound'] <= 3500, 'invalid full-schedule budget')


def comparison(rows, candidate, contract, *, draws=None):
    """Pair online observations with their own qualified IC reference."""
    from tournament import comparison as paired
    draws = contract['bootstrap_draws'] if draws is None else draws
    match_support = contract['comparison_kind'] != 'factor-base-policy'
    cold = paired(rows, candidate, baseline='incumbent', draws=draws, match_support=match_support)
    if not cold.get('eligible'):
        return cold
    online = paired(rows, candidate, baseline='ic_online', draws=draws, match_support=match_support)
    if not online.get('eligible'):
        return dict(candidate=candidate, eligible=False,
                    reason='online reference comparison: ' + online['reason'])
    require(cold['paired_cases'] == online['paired_cases']
            and set(cold['per_cell']) == set(online['per_cell']), 'unmatched reference workloads')
    cold['online'] = online['online']
    cold['metric_references'] = copy.deepcopy(METRIC_REFERENCES)
    return cold


def final_comparison(rows, candidate, contract):
    """Keep the original paired-target bootstrap, with a metric-specific baseline."""
    validate_contract(contract)
    result = comparison(rows, candidate, contract, draws=20)
    if not result.get('eligible'):
        return result
    grouped = {}
    for row in rows:
        if row['arm'] in (candidate, 'incumbent', 'ic_online'):
            grouped.setdefault((row['cell'], row['case'], row['arm']), []).append(row)
    logs = {}
    for cell, case in sorted({(row['cell'], row['case']) for row in rows}):
        samples = {alias: grouped[(cell, case, alias)] for alias in (candidate, 'incumbent', 'ic_online')}
        require(all(len(values) == 3 for values in samples.values()),
                'final inference requires three process repetitions')
        def costs(alias):
            values = samples[alias]
            return [statistics.median(item['total_operations'] for item in values),
                    statistics.median(item['measurement']['native_timing']['cold']['wall_ns'] for item in values),
                    statistics.median(item['measurement']['native_timing']['online']['wall_ns'] for item in values)]
        incumbent, online, challenger = costs('incumbent'), costs('ic_online'), costs(candidate)
        baseline = incumbent[:2] + [online[2]]
        require(all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in baseline + challenger),
                'missing positive complete timing/cost')
        logs.setdefault(cell, []).append([math.log(y / x) for x, y in zip(baseline, challenger)])
    cells = sorted(logs)
    require(cells == sorted(contract['confirmation_cases_per_cell']), 'changed final curve panel')
    require(all(len(logs[cell]) == contract['confirmation_cases_per_cell'][cell] for cell in cells),
            'changed final target allocation')
    estimates = [statistics.mean(statistics.mean(row[j] for row in logs[cell]) for cell in cells)
                 for j in range(3)]
    rng = random.Random(RULE['bootstrap_seed'] + contract['attempt_number'])
    samples = [[], [], []]
    draws = RULE['draws']
    for _ in range(draws):
        means = [[], [], []]
        for cell in cells:
            values = logs[cell]
            chosen = rng.choices(values, k=len(values))
            for j in range(3):
                means[j].append(math.fsum(row[j] for row in chosen) / len(chosen))
        for j in range(3):
            samples[j].append(math.fsum(means[j]) / len(cells))
    tail = max(0, math.floor(draws / 360) - 1)
    uncertainty = {}
    for j, metric in enumerate(METRICS):
        samples[j].sort()
        uncertainty[metric] = dict(ratio=math.exp(estimates[j]),
            upper=math.exp(2 * estimates[j] - samples[j][tail]),
            per_cell={cell: math.exp(statistics.mean(row[j] for row in logs[cell])) for cell in cells},
            descriptive_ci95=[math.exp(2 * estimates[j] - samples[j][int(.975 * draws)]),
                              math.exp(2 * estimates[j] - samples[j][int(.025 * draws)])])
    result.update(familywise=dict(rule=copy.deepcopy(RULE), metrics=uncertainty,
        metric_references=copy.deepcopy(METRIC_REFERENCES), paired_targets=sum(map(len, logs.values())),
        fixed_cells=cells, monte_carlo_tail_index=tail))
    result['ci95'] = uncertainty['instructions']['descriptive_ci95']
    result['candidate_over_baseline'] = uncertainty['instructions']['ratio']
    result['speedup'] = 1 / result['candidate_over_baseline']
    result['per_cell'] = uncertainty['instructions']['per_cell']
    result['native_wall_candidate_over_baseline'] = uncertainty['cold_ns']['ratio']
    result['native_wall_per_cell'] = uncertainty['cold_ns']['per_cell']
    result['native_wall_ci95'] = uncertainty['cold_ns']['descriptive_ci95']
    result['native_wall_status'] = 'fixed-cell paired-target basic bootstrap; complete cold native interval'
    result['online']['ci95'] = uncertainty['online_ns']['descriptive_ci95']
    result['online']['candidate_over_baseline'] = uncertainty['online_ns']['ratio']
    result['online']['per_cell'] = uncertainty['online_ns']['per_cell']
    return result


def promotion_passes(result, contract):
    validate_contract(contract)
    family = result.get('familywise', {})
    if (not result.get('eligible') or family.get('rule') != RULE
            or family.get('metric_references') != METRIC_REFERENCES
            or result.get('metric_references') != METRIC_REFERENCES):
        return False
    metrics = family.get('metrics', {})
    if set(metrics) != set(METRICS):
        return False
    for metric in METRICS:
        row = metrics[metric]
        values = [row['ratio'], row['upper'], *row['per_cell'].values()]
        if (not all(math.isfinite(v) and v > 0 for v in values)
                or set(row['per_cell']) != set(contract['confirmation_cases_per_cell'])
                or row['upper'] >= 1 or max(row['per_cell'].values()) > 1.1
                or row['ratio'] > (1 if metric == 'online_ns' else .8)):
            return False
    return True
