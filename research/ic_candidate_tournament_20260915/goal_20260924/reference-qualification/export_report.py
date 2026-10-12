"""Export review tables from frozen qualification receipts, without new measurements."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


def read(path):
    return json.loads(path.read_text())


def gm(values):
    return math.exp(statistics.mean(math.log(value) for value in values))


def number(value):
    # Display precision only; original full-precision summaries remain archived.
    return format(value, '.7g')


def export(root):
    qualification = read(root / 'qualification.json')
    fixtures = read(root / 'fixtures.json')['development']
    table = qualification['table']
    selected = {key: qualification[key] for key in (
        'selected_ic_cold', 'selected_ic_online', 'selected_rho_cold', 'selected_rho_online')}
    cold_reference = next(row for row in table if row['alias'] == selected['selected_rho_cold'])
    online_reference = next(row for row in table if row['alias'] == selected['selected_rho_online'])
    online_ref_ratio = online_reference['comparison_to_archived_incumbent']['online']['candidate_over_baseline']
    cold_ref_ratio = cold_reference['comparison_to_archived_incumbent']['candidate_over_baseline']
    summaries = []
    identities = {}
    stage_rows = []
    for row in table:
        alias = row['alias']
        costs = row['comparison_to_archived_incumbent']
        cells = row['per_cell']
        floor_ratios = []
        ids = {}
        for case in fixtures:
            receipts = [read(root / 'runs/development' / case['id'] / alias /
                             f'rep-{rep}/receipt.json') for rep in range(qualification['repetitions'])]
            assert all(item['status'] == 'VERIFIED' for item in receipts)
            identity_key = 'candidate_id' if row['mode'] == 'ic' else 'reference_id'
            ids[case['cell']] = sorted(set(ids.get(case['cell'], [])) |
                                     {item['measurement'][identity_key] for item in receipts})
            if row['mode'] == 'ic':
                floor_ratios.append(statistics.median(item['ratio_to_floor'] for item in receipts))
            for item in receipts:
                measurement = item['measurement']
                stage_rows.append(dict(alias=alias, case=case['id'], cell=case['cell'],
                    repetition=item['repetition'], mode=row['mode'],
                    identity=measurement[identity_key],
                    workload_id=measurement['workload_id'], run_id=measurement['run_id'],
                    status=item['status'], measurement=measurement))
        identities[alias] = ids
        summaries.append(dict(alias=alias, mode=row['mode'], verified_runs=row['verified_runs'],
            online_ms=number(gm([value['online_ns'] for value in cells.values()]) / 1e6),
            online_over_old_ic=number(costs['online']['candidate_over_baseline']),
            online_ratio_ci95=[number(value) for value in costs['online']['ci95']],
            selected_rho_over_online=number(online_ref_ratio / costs['online']['candidate_over_baseline']),
            cold_ms=number(gm([value['cold_ns'] for value in cells.values()]) / 1e6),
            cold_time_over_old_ic=number(costs['native_wall_candidate_over_baseline']),
            cold_time_ratio_ci95=[number(value) for value in costs['native_wall_ci95']],
            cold_Ir=number(gm([value['instructions'] for value in cells.values()])),
            S_Ir_per_sqrt_r=number(gm([value['normalized_S'] for value in cells.values()])),
            Ir_over_old_ic=number(costs['candidate_over_baseline']),
            Ir_ratio_ci95=[number(value) for value in costs['ci95']],
            Ir_over_selected_rho=number(costs['candidate_over_baseline'] / cold_ref_ratio),
            Ir_over_K_floor=number(gm(floor_ratios)) if floor_ratios else None))
    online = read(root / 'summaries/development.json')['single_target_online']
    paired = [row for row in online if row['rho_alias'] == selected['selected_rho_online']]
    assert len(stage_rows) == 945 and len(paired) == 45
    assert len({(row['identity'], row['workload_id'], row['run_id']) for row in stage_rows}) == 945
    inputs = ['contract.json', 'fixtures.json', 'qualification.json', 'summaries/development.json']
    return dict(schema_version=1, classification='accounting/reference qualification; no promotion',
        full_qualification_pairs=1290, development_pairs=945, improvement_rounds_used=0,
        promotion_eligible=False, input_sha256={name:hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in inputs}, selected=selected, summary=summaries, identities=identities,
        aggregation='Median of three repetitions per point, then equal-cell geometric mean; no target amortization.',
        uncertainty='Descriptive development bootstrap only; no familywise or held-out significance claim.',
        qualification=qualification, single_target_online=paired, development_runs=stage_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--round', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = export(args.round)
    if args.verify:
        if read(args.out) != data:
            raise ValueError('Export differs from the frozen receipts')
    else:
        with args.out.open('x') as stream:
            json.dump(data, stream, indent=2)
            stream.write('\n')
    print(json.dumps(dict(status='VERIFIED' if args.verify else 'EXPORTED',
        variants=len(data['summary']), development_runs=len(data['development_runs']),
        paired_single_target_rows=len(data['single_target_online']))))


if __name__ == '__main__':
    main()
