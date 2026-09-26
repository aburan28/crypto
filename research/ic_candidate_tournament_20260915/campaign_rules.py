"""Frozen rules for the September 24 bounded improvement campaign.

Historical evaluators remain unchanged. These rules bind accepted reference
evidence and define final-stage inference on a fixed, finite curve panel.
"""
import hashlib
import json
import math
import random
import statistics

from oracle import require


PURPOSE = 'bounded-improvement-20260924-v1'
QUALIFICATION_SHA256 = 'a043cf343b9bf78bacc4c5c5676f1fc05bcdabf68665fb0175f9f22a495cbc4e'
REFERENCE_ARCHIVE_SHA256 = 'a61e6b7c8154dfc594365a044dc74346880018a9c57e14d699ab4c55444950ac'
HISTORY_SHA256 = '32d684587edc886d4662f381f94f6d4e6ec8ad1bb78132ee2265711acdde60d7'
CELLS = ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0']
HOLDOUT_CELLS = ['n29a1']
IC_SOURCE = '240b8daa0478aacb1f6fc248de9fd5ed9b0773c59d3ae869ced0bab8be438378'
COLD_RHO_SOURCE = '65df2c4417ec963a2b03e511719df28df51ea266a3c4a9e825307b8d1ac8c053'
CONFIG = dict(solver='pair_table', linear_algebra='tiny_gauss', summands=3,
              batch_trials=1, max_trials=65536)
METRICS = ('instructions', 'cold_ns', 'online_ns')
RULE = dict(schema_version=1, campaign=PURPOSE, maximum_attempts=3,
    final_stages=['confirmation', 'replay'], metrics=list(METRICS),
    family_alpha='0.05', comparisons=18, one_sided_alpha='1/360',
    draws=50000, bootstrap_seed=914223,
    estimator='equal-cell mean log ratio of per-point three-process medians',
    resampling='paired targets within each fixed cell; same sampled indices for every metric',
    interval='one-sided basic bootstrap upper bound on log ratio; Bonferroni across 3 attempts x 2 stages x 3 metrics',
    coverage='nominal familywise 95%; bootstrap coverage is approximate, conditional on the declared target law and fixed cells')


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def qualified_binding(report, incumbent, references):
    require(object_hash(report) == QUALIFICATION_SHA256, 'unaccepted qualification report')
    table = {row['alias']: row for row in report['table']}
    require(len(references) == 2, 'both qualified rho references are required')
    declarations = [('incumbent', report['selected_ic_cold'], incumbent),
                    ('rho', report['selected_rho_cold'], references[0]),
                    ('rho_online', report['selected_rho_online'], references[1])]
    require(report['selected_ic_cold'] == report['selected_ic_online'],
            'this protocol requires the accepted common IC reference')
    bindings = {}
    for alias, selected, arm in declarations:
        expected = table[selected]
        require(arm['id'] == alias and
                arm['source_manifest_sha256'] == expected['source_manifest_sha256'] and
                arm['config'] == expected['config'], 'changed qualified reference '+alias)
        bindings[alias] = dict(selected_alias=selected,
            source_manifest_sha256=arm['source_manifest_sha256'], configuration=arm['config'])
    return dict(qualification_sha256=QUALIFICATION_SHA256,
                archive_sha256=REFERENCE_ARCHIVE_SHA256, bindings=bindings)


def validate_contract(contract):
    require(contract.get('purpose') == PURPOSE, 'not the bounded improvement protocol')
    require(contract.get('familywise_rule') == RULE, 'changed familywise confirmation rule')
    require(type(contract.get('attempt_number')) is int and 1 <= contract['attempt_number'] <= 3,
            'attempt is outside the three-round budget')
    require(contract['seed'] == 2026092550 + contract['attempt_number'], 'changed predeclared round seed')
    validate_environment(system=contract['host']['system'], machine=contract['host']['machine'],
                         compiler=contract['compiler'], profiler=contract['profiler_version'])
    require(contract['confirmation_cases'] == 72 and contract['repetitions'] == 3,
            'incomplete confirmation sample/repetitions')
    require(contract['cells'] == CELLS and contract['holdout_cells'] == HOLDOUT_CELLS and
            contract['confirmation_cases_per_cell'] == {cell:12 for cell in CELLS + HOLDOUT_CELLS},
            'changed predeclared six-cell confirmation allocation')
    require(contract['confirmation_ratio'] == .8 and contract['max_cell_ratio'] == 1.1
            and contract['require_native_progress'] and contract['objective'] == 'incumbent',
            'changed improvement objective')
    require(contract['selection_width'] == 6 and contract['exploration_slots'] == 1,
            'changed diverse portfolio budget')
    require(contract['limits']['max_profiled_jobs'] == 3500 and
            contract['limits']['memory_bytes'] == 8 * 1024**3 and
            contract['limits']['timeout_seconds'] == 180 and
            contract['limits']['worker_threads'] == 1 and contract['target_count'] == 1,
            'changed predeclared resource or target limits')
    binding = contract.get('reference_qualification')
    require(isinstance(binding, dict) and binding.get('qualification_sha256') == QUALIFICATION_SHA256
            and binding.get('archive_sha256') == REFERENCE_ARCHIVE_SHA256
            and set(binding.get('bindings', {})) == {'incumbent', 'rho', 'rho_online'},
            'missing qualified reference binding')
    for name, alias, source, config in (
            ('incumbent', 'pairinv', IC_SOURCE, CONFIG),
            ('rho', 'rho_incumbent_4', COLD_RHO_SOURCE, dict(CONFIG, rho_parallel_walks=4)),
            ('rho_online', 'rho_pairinv_4', IC_SOURCE, dict(CONFIG, rho_parallel_walks=4))):
        require(binding['bindings'][name] == dict(selected_alias=alias,
            source_manifest_sha256=source, configuration=config), 'changed bound reference '+name)


def validate_environment(*, system, machine, compiler, profiler):
    require(system == 'Linux' and machine == 'x86_64' and
            compiler.startswith('rustc 1.94.1 ') and profiler == 'valgrind-3.22.0',
            'bounded campaign requires its registered Linux amd64/compiler/profiler environment')


def validate_preparation(args):
    """Reject a changed budget/panel before compiling or exposing fresh targets."""
    require(1 <= args.attempt_number <= 3 and args.seed == 2026092550 + args.attempt_number,
            'invalid bounded round/seed')
    require(not args.qualification and args.profile == 'pilot' and args.confirmation_cases == '',
            'bounded campaign requires the fixed pilot allocation')
    require(args.cells.split(',') == [c[1:] for c in CELLS] and
            args.holdout_cells.split(',') == [c[1:] for c in HOLDOUT_CELLS],
            'changed predeclared curve panel')
    require(args.selection_width == 6 and args.exploration_slots == 1 and
            args.timeout == 180 and args.max_processes == 3500 and args.targets == 1 and
            args.require_native_progress and args.objective == 'incumbent' and
            args.comparison_kind == 'factor-base-policy', 'changed bounded preparation limits')


def final_comparison(rows, candidate, contract, *, baseline='incumbent'):
    """Final-stage paired inference. A process repetition is not a new target."""
    from tournament import comparison
    validate_contract(contract)
    # Reuse existing fixture/support/completion validation before extracting data.
    result = comparison(rows, candidate, baseline=baseline, draws=20,
                        match_support=contract['comparison_kind'] != 'factor-base-policy')
    if not result.get('eligible'):
        return result
    grouped = {}
    for row in rows:
        if row['arm'] in (candidate, baseline):
            grouped.setdefault((row['cell'], row['case'], row['arm']), []).append(row)
    logs = {}
    for cell, case in sorted({(row['cell'], row['case']) for row in rows}):
        a, b = grouped[(cell, case, baseline)], grouped[(cell, case, candidate)]
        require(len(a) == len(b) == 3, 'final inference requires three process repetitions')
        def values(items):
            return [statistics.median(item['total_operations'] for item in items),
                    statistics.median(item['measurement']['native_timing']['cold']['wall_ns'] for item in items),
                    statistics.median(item['measurement']['native_timing']['online']['wall_ns'] for item in items)]
        av, bv = values(a), values(b)
        require(all(type(v) in (int, float) and math.isfinite(v) and v > 0 for v in av + bv),
                'missing positive complete timing/cost')
        logs.setdefault(cell, []).append([math.log(y / x) for x, y in zip(av, bv)])
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
    # Use the lower order statistic for a conservative Monte Carlo tail choice.
    tail = max(0, math.floor(draws / 360) - 1)
    uncertainty = {}
    for j, metric in enumerate(METRICS):
        samples[j].sort()
        uncertainty[metric] = dict(ratio=math.exp(estimates[j]),
            upper=math.exp(2 * estimates[j] - samples[j][tail]),
            per_cell={cell:math.exp(statistics.mean(row[j] for row in logs[cell])) for cell in cells},
            descriptive_ci95=[math.exp(2 * estimates[j] - samples[j][int(.975 * draws)]),
                              math.exp(2 * estimates[j] - samples[j][int(.025 * draws)])])
    result.update(familywise=dict(rule=RULE, metrics=uncertainty, paired_targets=sum(map(len, logs.values())),
                                 fixed_cells=cells, monte_carlo_tail_index=tail))
    # Replace the small validation-only bootstrap with the reported fixed-panel one.
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
    if not result.get('eligible') or result.get('familywise', {}).get('rule') != RULE:
        return False
    metrics = result['familywise'].get('metrics', {})
    if set(metrics) != set(METRICS):
        return False
    for metric in METRICS:
        row = metrics[metric]
        values = [row['ratio'], row['upper'], *row['per_cell'].values()]
        if (not values or not all(math.isfinite(v) and v > 0 for v in values) or
                set(row['per_cell']) != set(contract['confirmation_cases_per_cell']) or
                row['upper'] >= 1 or max(row['per_cell'].values()) > 1.1 or
                row['ratio'] > (1 if metric == 'online_ns' else .8)):
            return False
    return True


def selection_key(row):
    """Choose one balanced cold-cost challenger, preserving online priority in ties."""
    return (max(row['candidate_over_baseline'], row['native_wall_candidate_over_baseline']),
            row['online']['candidate_over_baseline'], row['candidate'])
