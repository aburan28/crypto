"""Export the registered development panel without executing or selecting new work."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def gm(values):
    return math.exp(statistics.mean(math.log(value) for value in values))


def number(value):
    # Display precision only; the archive retains all original full-precision data.
    return format(value, '.7g')


def paired_rho_intervals(campaign, aliases):
    # Reuse the frozen estimator and paired resampling, with rho as numerator.
    # Ratios of marginal interval endpoints would lose target covariance.
    program = """import json, sys
from pathlib import Path
root = Path(sys.argv[1])
sys.path.insert(0, str(root/'evaluator'))
from tournament import read, comparison
contract, report = read(root/'contract.json'), read(root/'qualification.json')
rows = [read(p) for p in sorted(root.glob('runs/development/**/receipt.json'))]
reference = report['selected_rho_online']
result = {alias: comparison(rows, reference, baseline=alias,
    draws=contract['bootstrap_draws'], match_support=False)
    for alias in json.loads(sys.argv[2])} if reference else {}
print(json.dumps(result))
"""
    result = subprocess.run([sys.executable, '-c', program, str(campaign), json.dumps(aliases)],
                            capture_output=True, text=True, timeout=300)
    require(result.returncode == 0, 'paired rho analysis failed:\n' + result.stdout + result.stderr)
    return json.loads(result.stdout)


def export(root):
    campaign = root / 'tournament'
    panel, observed = read(root / 'summary.json'), read(root / 'observer/summary.json')
    qualification = read(campaign / 'qualification.json')
    contract = read(campaign / 'contract.json')
    aa = read(campaign / 'summaries/aa.json')
    require(aa['passed'] and aa['verified_runs'] == 30 and not aa['failures'],
            'registered A/A gate did not pass')
    fixtures = read(campaign / 'fixtures.json')['development']
    require(panel['trial_slots'] == 1350 and panel['observer_pairs'] == 360
            and panel['accepted_reference_binding_changed'] is False
            and panel['promotion_eligible'] is False, 'changed panel or acceptance scope')
    require(observed['pairs'] == 360 and observed['failed_pairs'] == 0
            and observed['legacy_scientific_admission'] is None
            and observed['promotion_eligible'] is False, 'observer did not pass its registered gate')
    require(qualification['eligible_for_improvement'] is False
            and qualification['observer_qualification'] is None
            and qualification['promotion_eligible'] is False, 'unearned competitive qualification')
    require(len(fixtures) == 15 and len(qualification['table']) == 22
            and qualification['repetitions'] == 3, 'changed development panel')
    selected = {key: qualification[key] for key in (
        'selected_ic_cold', 'selected_ic_online', 'selected_rho_cold', 'selected_rho_online')}
    table = {row['alias']: row for row in qualification['table']}
    rho_comparisons = paired_rho_intervals(campaign, [a for a, row in table.items() if row['mode'] == 'ic'])
    online_reference = table.get(selected['selected_rho_online'])
    cold_reference = table.get(selected['selected_rho_cold'])
    online_ref_ratio = (online_reference['comparison_to_archived_incumbent']['online']['candidate_over_baseline']
                        if online_reference else None)
    cold_ref_ratio = (cold_reference['comparison_to_archived_incumbent']['candidate_over_baseline']
                     if cold_reference else None)
    summaries, identities, stage_rows = [], {}, []
    for alias, row in table.items():
        costs, cells = row['comparison_to_archived_incumbent'], row['per_cell']
        floors, ids, inventory = [], {}, {}
        for case in fixtures:
            receipts = [read(campaign / 'runs/development' / case['id'] / alias /
                             f'rep-{rep}/receipt.json') for rep in range(3)]
            identity_key = 'candidate_id' if row['mode'] == 'ic' else 'reference_id'
            ids[case['cell']] = sorted(set(ids.get(case['cell'], [])) |
                                      {item['measurement'][identity_key] for item in receipts})
            if row['mode'] == 'ic':
                admission = read(campaign / 'admissions/development' / case['id'] / alias / 'admission.json')
                base = admission['candidate']['record']['factor_base']['inventory']
                observed_base = dict(usable_points=base['usable_point_count'], folded_columns=base['effective_columns'])
                require(inventory.setdefault(case['cell'], observed_base) == observed_base,
                        'factor base changed within a cell')
                if all(item['status'] == 'VERIFIED' for item in receipts):
                    floors.append(statistics.median(item['ratio_to_floor'] for item in receipts))
            for item in receipts:
                measurement = item['measurement']
                stage_rows.append(dict(alias=alias, case=case['id'], cell=case['cell'],
                    repetition=item['repetition'], mode=row['mode'], identity=measurement[identity_key],
                    workload_id=measurement['workload_id'], run_id=measurement['run_id'],
                    status=item['status'], reason=item.get('reason'), measurement=measurement))
        identities[alias] = ids
        # A failed arm keeps its rows and null aggregate costs; no successful-subset estimate.
        complete = row['qualified'] and costs.get('eligible') and len(cells) == 5
        values = dict.fromkeys(('online_ms', 'online_over_incumbent', 'online_ratio_ci95',
            'selected_rho_over_online', 'selected_rho_over_online_ci95',
            'cold_ms', 'cold_time_over_incumbent', 'cold_time_ratio_ci95',
            'cold_Ir', 'S_Ir_per_sqrt_r', 'Ir_over_incumbent', 'Ir_ratio_ci95',
            'Ir_over_selected_rho', 'Ir_over_K_floor'))
        phase_shares = None
        if complete:
            values.update(online_ms=number(gm([v['online_ns'] for v in cells.values()]) / 1e6),
                online_over_incumbent=number(costs['online']['candidate_over_baseline']),
                online_ratio_ci95=[number(v) for v in costs['online']['ci95']],
                selected_rho_over_online=(number(online_ref_ratio / costs['online']['candidate_over_baseline'])
                                          if online_ref_ratio else None),
                cold_ms=number(gm([v['cold_ns'] for v in cells.values()]) / 1e6),
                cold_time_over_incumbent=number(costs['native_wall_candidate_over_baseline']),
                cold_time_ratio_ci95=[number(v) for v in costs['native_wall_ci95']],
                cold_Ir=number(gm([v['instructions'] for v in cells.values()])),
                S_Ir_per_sqrt_r=number(gm([v['normalized_S'] for v in cells.values()])),
                Ir_over_incumbent=number(costs['candidate_over_baseline']),
                Ir_ratio_ci95=[number(v) for v in costs['ci95']],
                Ir_over_selected_rho=(number(costs['candidate_over_baseline'] / cold_ref_ratio)
                                      if cold_ref_ratio else None),
                Ir_over_K_floor=number(gm(floors)) if len(floors) == len(fixtures) else None)
            if row['mode'] == 'ic' and online_reference:
                paired_cost = rho_comparisons[alias]
                require(paired_cost['eligible'], 'complete arm lost its paired rho observations')
                require(math.isclose(paired_cost['online']['candidate_over_baseline'],
                        online_ref_ratio / costs['online']['candidate_over_baseline'], rel_tol=1e-12),
                        'paired rho point estimate differs')
                values['selected_rho_over_online'] = number(paired_cost['online']['candidate_over_baseline'])
                values['selected_rho_over_online_ci95'] = [number(v) for v in paired_cost['online']['ci95']]
            if row['mode'] == 'ic':
                measured = [r['measurement'] for r in stage_rows if r['alias'] == alias]
                require(len(measured) == 45 and all(m['phase_ledger']['complete'] for m in measured),
                        'phase diagnostics need the full verified workload')
                phase_shares = {phase: number(statistics.mean(
                    m['phase_ledger']['operations'][phase] / m['total_operations'] for m in measured))
                    for phase in measured[0]['phase_ledger']['operations']}
        summaries.append(dict(alias=alias, mode=row['mode'], qualified=row['qualified'],
            verified_runs=row['verified_runs'], scheduled_runs=row['scheduled_runs'],
            smoke_failures=row['smoke_failures'], failures=row['failures'], base_inventory=inventory,
            cold_instruction_phase_shares=phase_shares,
            **values))
    online = read(campaign / 'summaries/development.json')['single_target_online']
    paired = [row for row in online if row['rho_alias'] == selected['selected_rho_online']]
    require(len(stage_rows) == 990 and len(paired) == (60 if online_reference else 0), 'changed export census')
    require(len({(r['identity'], r['workload_id'], r['run_id']) for r in stage_rows}) == 990,
            'colliding development run keys')
    inputs = ['summary.json', 'comparison-summary.json', 'tournament/contract.json',
              'tournament/fixtures.json', 'tournament/qualification.json',
              'tournament/summaries/aa.json', 'tournament/summaries/development.json',
              'observer/contract.json', 'observer/summary.json']
    return dict(schema_version=1, classification='accounting/reference development comparison; no promotion',
        full_trial_slots=1350, development_slots=990, observer_pairs=360, improvement_rounds_used=0,
        promotion_eligible=False, accepted_reference_binding_changed=False,
        input_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in inputs},
        panel=panel, observer=observed, selected=selected, summary=summaries, identities=identities,
        aa=aa, environment={key: contract[key] for key in
            ('host', 'limits', 'compiler', 'unit', 'native_timing_protocol')},
        paired_rho_interval_method='Frozen tournament.comparison with rho as numerator; paired curve/target resampling, original draw count and seed; display precision only.',
        aggregation='Three-process medians per supplied point, then equal-cell geometric means; no target amortization.',
        phase_share_aggregation='Mean per-run fraction of complete cold Ir; the fixed equal target/repetition counts give equal cell and target weight. These are stage diagnostics, not speedups.',
        uncertainty='Descriptive development intervals only; selection does not provide a held-out or familywise claim.',
        observer_scope='Whole-mode effects on outer intervals; no subtraction, low-overhead claim or legacy scientific admission.',
        qualification=qualification, single_target_online=paired, development_runs=stage_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = export(args.bundle)
    if args.verify:
        require(read(args.out) == result, 'export differs from frozen evidence')
    else:
        with args.out.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    print(json.dumps(dict(status='VERIFIED' if args.verify else 'EXPORTED', variants=22,
                         development_runs=990, observer_pairs=360, promotion_eligible=False)))


if __name__ == '__main__':
    main()
