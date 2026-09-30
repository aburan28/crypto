"""Read closed F5/SAT evidence; audit counters and mathematical S3 witnesses.

No producer or native solver is invoked. Mathematical chain representation
does not prove the implemented Boolean ANF or matrix criterion is faithful.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys

from generic_stages import projected_columns
from identity import curve_record, sha256, write_immutable
from oracle import Curve, rank, require
from replay_paired_n17_evidence import retained_files

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/'goal_20260924/f5-closed-diagnostics-20260930/protocol.json'
INTEGER_COUNTERS = {'reductions', 'splits', 'propagations', 'infeasible_branches',
                    'max_degree_built', 'oversize', 'eliminated'}
FLAGS = {'exhausted', 'unsupported'}


def s3(curve, x, y, z):
    product = curve.fm(x, y)
    return (curve.fm(curve.fm(x ^ y, x ^ y), curve.fm(z, z))
            ^ curve.fm(product, z) ^ curve.fm(product, product) ^ 1)


def chain_control(curve, base, indices, target):
    require(target is not None and len(indices) == 3
            and all(type(index) is int and 0 <= index < len(base) for index in indices),
            'malformed disclosed finite three-point witness')
    points = [base[index] for index in indices]
    require(all(point is not None for point in points), 'identity summand is not a geometric base point')
    total = None
    for point in points:
        total = curve.add(total, point)
    require(total == target, 'disclosed points do not re-add to the public query')
    orderings = []
    for order in itertools.permutations(range(3)):
        first, second, third = [points[index] for index in order]
        intermediate = curve.add(first, second)
        if intermediate is None:
            orderings.append(dict(order=list(order), intermediate=None,
                                  status='IDENTITY_INTERMEDIATE', s3_values=None,
                                  represented_by_finite_chain=False))
            continue
        values = [s3(curve, first[0], second[0], intermediate[0]),
                  s3(curve, intermediate[0], third[0], target[0])]
        require(values == [0, 0], 'independent mathematical S3 links reject a full-point witness')
        orderings.append(dict(order=list(order), intermediate=list(intermediate),
                              status='FINITE_CHAIN', s3_values=values,
                              represented_by_finite_chain=True))
    return dict(indices=list(indices), points=[list(point) for point in points],
                target=list(target), full_point_replay=True, orderings=orderings,
                finite_chain_exists=any(row['represented_by_finite_chain'] for row in orderings),
                implemented_boolean_anf_validated=False)


def counters(attempts, config, expected_fields):
    ledger = []
    by_outcome = defaultdict(list)
    require(set(expected_fields) == INTEGER_COUNTERS | FLAGS,
            'counter analysis protocol changed')
    for trial, attempt in enumerate(attempts):
        require(attempt['trial'] == trial and attempt['b'] == 0,
                'closed ordinary query chronology changed')
        pdp = attempt['pdp']
        observed = pdp['stats']
        require(set(observed) == {'family', 'engine', 'stats'}
                and observed['family'] == 'groebner'
                and observed['engine'] == {'MatrixF5': {'max_degree': config['groebner_degree']}},
                'closed F5 dispatch changed')
        stats = observed['stats']
        require(set(stats) == set(expected_fields)
                and all(type(stats[key]) is int and stats[key] >= 0 for key in INTEGER_COUNTERS)
                and all(type(stats[key]) is bool for key in FLAGS),
                'retained F5 counters missing or malformed')
        outcome = pdp['outcome']
        require(outcome in {'witness', 'incomplete', 'unsupported', 'proved_unsat'}
                and (pdp['points'] is not None) == (outcome == 'witness')
                and (outcome != 'incomplete' or (stats['exhausted'] and not stats['unsupported']))
                and (outcome != 'unsupported' or (stats['unsupported'] and stats['exhausted']))
                and (outcome != 'proved_unsat' or (not stats['unsupported'] and not stats['exhausted'])),
                'closed F5 verdict contradicts retained completion flags')
        row = dict(trial=trial, scalar=attempt['a'], outcome=outcome, stats=stats,
                   registered_reduction_budget=config['node_budget'],
                   at_reduction_budget=stats['reductions'] == config['node_budget'])
        ledger.append(row)
        by_outcome[outcome].append(row)
    summaries = {}
    for outcome, rows in sorted(by_outcome.items()):
        summaries[outcome] = dict(
            attempts=len(rows), at_reduction_budget=sum(row['at_reduction_budget'] for row in rows),
            # Built degree is a parameter, not an additive work counter.
            totals={key: sum(row['stats'][key] for row in rows)
                    for key in sorted(INTEGER_COUNTERS - {'max_degree_built'})},
            ranges={key: [min(row['stats'][key] for row in rows), max(row['stats'][key] for row in rows)]
                    for key in sorted(INTEGER_COUNTERS)},
            built_degree_histogram=dict(sorted(Counter(str(row['stats']['max_degree_built'])
                                                       for row in rows).items())))
    return dict(ledger=ledger, by_outcome=summaries,
                counter_scope='retained native counters; no conversion into a common operation unit',
                unresolved_instrumentation=['cache_hits_and_misses', 'per_query_matrix_word_counts',
                                           'observed_resolved_variable_and_value_order'])


def diagnose(f5, sat, job, protocol):
    require(job['mode'] == 'ic' and job['degree'] == 17 and job['curve_a'] == 1
            and job['factor_base'] == dict(kind='standard_subspace', dimension=6)
            and job['config']['solver'] == 'f5' and job['config']['summands'] == 3,
            'diagnosis is only the closed n17a1 standard-subspace F5 control')
    config = f5['effective_config']
    require(all(config[key] == value for key, value in job['config'].items())
            and config['node_budget'] == 4096 and config['groebner_degree'] == 3,
            'closed F5 limits changed')
    attempts = [attempt for batch in f5['collection_reports'] for attempt in batch['attempts']]
    require(len(attempts) == protocol['expected_f5_attempts']
            and len(sat['collection']) == protocol['expected_sat_prefix']
            and not f5.get('solutions') and f5['descent_dispatch'] is None,
            'closed streams or unstarted descent changed')
    report = counters(attempts, config, protocol['counter_fields'])
    curve = Curve(f5['fixture'])
    base = [curve.decode(point) for point in f5['factor_base']]
    require(len(base) == 63 and len(set(base)) == len(base), 'closed geometric base changed')
    columns, mapping = projected_columns(curve, base)
    direction = protocol['matrix_direction']
    require(len(columns) == f5['columns'] == 29
            and list(columns[direction['column']]) == direction['representative'],
            'closed orbit direction changed')
    matrix = f5['relation_matrix']
    require(matrix['modulus'] == str(curve.r)
            and [list(map(int, point)) for point in matrix['column_points']] == [list(p) for p in columns]
            and [list(map(int, point)) for point in sat['matrix']['column_points']] == [list(p) for p in columns],
            'closed matrices change subgroup or orbit columns')
    rows = []
    for item in matrix['rows']:
        row = [0]*len(columns)
        for column, value in item['entries']:
            row[column] = int(value) % curve.r
        rows.append(row)
    require(rank(rows, len(columns), curve.r) == 28
            and all(row[direction['column']] == 0 for row in rows),
            'closed F5 missing direction changed')
    witness_controls = []
    for attempt in attempts:
        if attempt['pdp']['points'] is not None:
            witness_controls.append(dict(arm='f5', trial=attempt['trial'],
                control=chain_control(curve, base, attempt['pdp']['points'],
                                      curve.mul(curve.g, attempt['a']))))
    selected = []
    for trial, row in enumerate(sat['collection']):
        attempt = attempts[trial]
        require(row['trial'] == trial and row['probe_scalar'] == attempt['a']
                and curve.decode(row['public_point']) == curve.mul(curve.g, attempt['a']),
                'retrospective SAT/F5 ordinary inputs differ')
        if row['point_witness'] is None:
            continue
        require(row['status'] == 'VALID_POINT_WITNESS' and row['source_model_valid'] is True,
                'retained SAT witness lacks its original validity verdict')
        witness = row['point_witness']
        indices = witness['point_indices']
        control = chain_control(curve, base, indices, curve.mul(curve.g, attempt['a']))
        require([curve.decode(point) for point in witness['points']] == [base[index] for index in indices],
                'SAT/F5 geometric point indices differ')
        witness_controls.append(dict(arm='sat-v2', trial=trial, control=control))
        if trial in protocol['selected_exposed_trials']:
            coefficient = sum(mapping[index][1] for index in indices
                              if mapping[index] is not None and mapping[index][0] == direction['column']) % curve.r
            require(coefficient != 0 and attempt['pdp']['outcome'] == 'incomplete',
                    'selected disclosed witness no longer fills a lost direction')
            selected.append(dict(trial=trial, scalar=attempt['a'],
                                 f5_outcome=attempt['pdp']['outcome'],
                                 f5_stats=attempt['pdp']['stats']['stats'],
                                 missing_column_coefficient=str(coefficient), control=control))
    require([row['trial'] for row in selected] == protocol['selected_exposed_trials'],
            'selected disclosed controls missing')
    usable = {curve.mul(point, curve.h) for point in base} - {None}
    require(len(usable) == 62, 'closed actual distinct usable base size changed')
    report.update(geometric_points=len(base), actual_usable_points=len(usable),
                  folded_columns=len(columns), independently_rebuilt_rank=28,
                  witnessed_controls=witness_controls, selected_exposed_controls=selected,
                  witness_control_count=len(witness_controls),
                  all_witnesses_have_a_finite_mathematical_chain=all(
                      row['control']['finite_chain_exists'] for row in witness_controls),
                  identity_intermediate_orderings=sum(
                      ordering['status'] == 'IDENTITY_INTERMEDIATE'
                      for row in witness_controls for ordering in row['control']['orderings']))
    return report


def analyze(bundle, protocol_file=PROTOCOL):
    bundle = Path(bundle)
    protocol = json.loads(Path(protocol_file).read_text())
    receipt = json.loads((bundle/'receipt.json').read_text())
    require(receipt['archive_sha256'] == protocol['archive_sha256']
            and receipt['archive_bytes'] == protocol['archive_bytes'],
            'diagnosis differs from externally fixed closed archive')
    files = retained_files(bundle)
    roles = ('f5/stdout.json', 'f5/f5-job.json', 'f5/candidate.json',
             'f5/f5-workload.json', 'f5/result.json', 'f5/build-record.json', 'sat-v2/summary.json')
    records = {role: json.loads(files[role]) for role in roles}
    candidate = records['f5/candidate.json']
    workload = records['f5/f5-workload.json']
    original = records['f5/result.json']
    require(sha256(candidate['record']) == candidate['record_sha256']
            and sha256(workload['record']) == workload['record_sha256']
            and candidate['candidate_id'] == original['candidate_id']
            and workload['workload_id'] == original['workload_id']
            and original['status'] == 'AUDITED_BOUNDED_INCOMPLETE'
            and original['online_wall_ns'] is None and original['paired_rho_speedup'] is None
            and original['promotion_eligible'] is False,
            'closed identity or failure classification changed')
    f5 = records['f5/stdout.json']
    require(candidate['record']['curve']['curve_id'] == curve_record(f5['fixture'])['curve']['curve_id'],
            'closed candidate curve differs from diagnostics')
    report = diagnose(f5, records['sat-v2/summary.json'], records['f5/f5-job.json'], protocol)
    loaded_sources = {}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            path = Path(path).resolve()
            if path.suffix == '.py' and path.is_relative_to(HERE) and path.is_file():
                loaded_sources[path.relative_to(HERE).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(schema_version=1, presentation_schema_version=2,
                status='AUDITED_CLOSED_DIAGNOSTICS', protocol_sha256=sha256(protocol),
                protocol_document_sha256=hashlib.sha256(
                    Path(protocol_file).with_name('PROTOCOL.md').read_bytes()).hexdigest(),
                analysis_invocation=list(sys.argv),
                archive_sha256=protocol['archive_sha256'],
                input_sha256={role: hashlib.sha256(files[role]).hexdigest() for role in roles},
                candidate_id=original['candidate_id'], workload_id=original['workload_id'], run_id=original['run_id'],
                native_source_manifest_sha256=records['f5/build-record.json']['source_manifest_sha256'],
                native_worker_sha256=records['f5/build-record.json']['worker_sha256'],
                analysis_sources=dict(sorted(loaded_sources.items())),
                analysis_python=dict(version=platform.python_version(),
                                     executable_sha256=hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest()),
                analysis_scope='post-execution source/import inventory; no new measured candidate or complete-solver admission',
                diagnostics=report, new_candidate_id=None, measured_costs=None,
                implemented_boolean_anf_validated=False, fresh_qualification=False,
                headline_online_admissible=False, promotion_eligible=False, online_speedup=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, default=PROTOCOL)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists(), 'diagnosis output exists; never overwrite evidence')
    result = analyze(args.bundle, args.protocol)
    write_immutable(args.out, result)
    print(json.dumps({key: value for key, value in result.items() if key != 'diagnostics'}, sort_keys=True))
