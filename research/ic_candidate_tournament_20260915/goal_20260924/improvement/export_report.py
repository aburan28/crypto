"""Export a verified bounded round's complete tables; never run a new benchmark."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import random
import subprocess
import sys


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def gm(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def aggregate(rows, getter):
    groups = {}
    for row in rows:
        value = getter(row)
        if value is None:
            return None
        if not math.isfinite(value) or value <= 0:
            raise ValueError('nonpositive/nonfinite complete cost')
        groups.setdefault((row['cell'], row['case']), []).append(value)
    cells = {}
    for (cell, _), values in groups.items():
        cells.setdefault(cell, []).append(statistics.median(values))
    return gm([gm(v) for v in cells.values()])


def rate_interval(points, numerator, denominator, seed):
    """Descriptive target-cluster interval, never an independent-query binomial CI."""
    def ratio(sample):
        bottom = sum(p[denominator] for p in sample)
        return sum(p[numerator] for p in sample)/bottom if bottom else None
    estimate = ratio(points)
    if len(points) < 2 or estimate is None:
        return dict(value=estimate, ci95=None, targets=len(points),
                    scope='insufficient distinct targets for a descriptive interval')
    rng = random.Random(seed)
    values = [ratio(rng.choices(points, k=len(points))) for _ in range(2000)]
    if any(v is None for v in values):
        interval = None
    else:
        values.sort(); interval = [values[50], values[1950]]
    return dict(value=estimate, ci95=interval, targets=len(points),
                scope='descriptive percentile bootstrap of whole target/walk clusters; approximate coverage, especially with few targets')


def stage_diagnostics(rows):
    """Keep stage costs and stopped-walk yield distinct from complete-solve claims."""
    cells = {}
    for row in rows:
        cells.setdefault(row['cell'], []).append(row)
    result = {}
    for cell, selected in sorted(cells.items()):
        complete = all(r['status'] == 'VERIFIED' for r in selected)
        grouped = {}
        for row in selected:
            grouped.setdefault(row['case'], []).append(row)
        data = []
        if complete:
            for case, reps in sorted(grouped.items()):
                counts = [r['measurement']['diagnostics'] for r in reps]
                # Repetitions repeat a walk, and must not triple the query sample.
                if not all(v == counts[0] for v in counts):
                    raise ValueError('process repetitions changed deterministic walk diagnostics')
                d = counts[0]
                collection_cost = statistics.median(sum(r['phase_costs'][p] for p in
                    ('queries', 'pdp', 'relation_check', 'matrix_build', 'relation_la')) for r in reps)
                data.append(dict(case=case, attempts=d['attempts'], verified=d['verified_relations'],
                    novel=d['novel_rows'], rank=d['final_rank'], outcomes=d['outcomes'],
                    collection_Ir=collection_cost))
        phase_names = ('setup', 'isogeny', 'factor_base', 'precompute', 'queries', 'pdp', 'relation_check',
                       'matrix_build', 'relation_la', 'target_descent', 'recovery_check')
        costs = {}
        if complete:
            for phase in phase_names:
                # Arithmetic means of exclusive ledgers remain additive. They
                # are separate from the headline geometric paired estimand.
                costs[phase] = statistics.mean(r['measurement']['phase_ledger']['operations'][phase]
                                               for r in selected)
            if not math.isclose(sum(costs.values()), statistics.mean(r['total_operations'] for r in selected),
                                rel_tol=1e-12):
                raise ValueError('mean exclusive phase ledger does not close')
        result[cell] = dict(complete=complete, distinct_targets=len(grouped),
            per_point_walks=data if complete else None,
            verified_relations_per_attempt=rate_interval(data, 'verified', 'attempts', 8131) if complete else None,
            novel_rows_per_attempt=rate_interval(data, 'novel', 'attempts', 8132) if complete else None,
            collection_Ir_per_novel_row=rate_interval(data, 'collection_Ir', 'novel', 8133) if complete else None,
            arithmetic_mean_phase_Ir=costs if complete else None,
            arithmetic_mean_complete_cold_Ir=statistics.mean(r['total_operations'] for r in selected) if complete else None,
            base_only_peak_bytes=None,
            whole_process_peak_rss_bytes=max(((r.get('native_process') or {}).get('peak_rss_bytes') or 0 for r in selected), default=0) or None,
            scope='Stopped ordinary-query collection walks, not planted decompositions or iid query-yield estimates. Collection cost includes query/PDP/check/matrix/LA; no setup or target work. Base-only memory was not instrumented; whole-process RSS is separate.')
    return result


def export(root):
    c = read(root/'contract.json')
    decision = read(root/'decision.json')
    if c['purpose'] != 'bounded-improvement-20260924-v1':
        raise ValueError('not a bounded improvement round')
    fixtures = read(root/'fixtures.json')
    stages, runs = {}, []
    for stage in c['stages']:
        summary = read(root/'summaries'/f'{stage}.json')
        rows = [read(p) for p in sorted((root/'runs'/stage).glob('*/*/rep-*/receipt.json'))]
        if len(rows) != summary['runs']:
            raise ValueError('receipt count differs from audited summary')
        comparisons = {v['candidate']:v for v in summary['comparisons']}
        comparisons.update(summary.get('rho_comparisons', {}))
        aliases = sorted({r['arm'] for r in rows}, key=lambda a:(a!='incumbent',a))
        tables = []
        for alias in aliases:
            selected = [r for r in rows if r['arm'] == alias]
            if len(selected) != len(fixtures[stage])*c['repetitions']:
                raise ValueError('incomplete scheduled arm')
            complete = all(r['status']=='VERIFIED' and r['total_operations'] is not None for r in selected)
            mode = selected[0]['mode']
            timing = lambda r, part: r['measurement']['native_timing'][part]['wall_ns']
            get = lambda fn: aggregate(selected,fn) if complete else None
            identities = sorted({r['measurement'].get('candidate_id',r['measurement'].get('reference_id')) for r in selected})
            comp = comparisons.get(alias)
            online_ms = get(lambda r:timing(r,'online')/1e6)
            cold_ms = get(lambda r:timing(r,'cold')/1e6)
            cost = get(lambda r:r['total_operations'])
            normalized = get(lambda r:r['normalized_S'])
            floor = get(lambda r:r.get('ratio_to_floor')) if mode=='ic' else None
            online_ratio = 1 if alias=='incumbent' and complete else (comp or {}).get('online',{}).get('candidate_over_baseline')
            ir_ratio = 1 if alias=='incumbent' and complete else (comp or {}).get('candidate_over_baseline')
            cold_ratio = 1 if alias=='incumbent' and complete else (comp or {}).get('native_wall_candidate_over_baseline')
            shapes = {}
            for row in selected:
                m = row['measurement']; cert=m.get('certificate') or {}
                diag=m.get('diagnostics') or {}
                admission=read(root/'admissions'/stage/row['case']/alias/'admission.json')
                inventory=admission.get('candidate',{}).get('record',{}).get('factor_base',{}).get('inventory',{})
                base_size=inventory.get('usable_point_count')
                columns=inventory.get('effective_columns')
                shapes.setdefault(row['cell'],set()).add((base_size,columns,cert.get('rank')))
                runs.append(dict(stage=stage,alias=alias,cell=row['cell'],case=row['case'],repetition=row['repetition'],
                    identity=m.get('candidate_id',m.get('reference_id')),workload_id=m['workload_id'],run_id=m['run_id'],
                    status=row['status'],reason=row.get('reason'),
                    source_manifest_sha256=m['provenance']['source_manifest_sha256'],
                    B=base_size,columns=columns,rank=cert.get('rank'),
                    cold_Ir=row['total_operations'],cold_ns=timing(row,'cold') if row['status']=='VERIFIED' else None,
                    online_ns=timing(row,'online') if row['status']=='VERIFIED' else None,
                    S_Ir_per_sqrt_r=row.get('normalized_S'),Ir_over_K_floor=row.get('ratio_to_floor'),
                    peak_rss_bytes=(row.get('native_process') or {}).get('peak_rss_bytes'),
                    attempts=diag.get('attempts'),novel_rows=diag.get('novel_rows'),verified_relations=diag.get('verified_relations'),
                    phase_Ir=json.dumps((m.get('phase_ledger') or {}).get('operations')
                        if mode == 'ic' else row.get('phase_costs'), sort_keys=True,separators=(',',':')),
                    raw_instrumented_phase_Ir=json.dumps(row.get('phase_costs'),sort_keys=True,separators=(',',':')),
                    online_phase_ns=json.dumps(m['native_timing']['online']['phase_wall_ns'],sort_keys=True,separators=(',',':')) if m.get('native_timing') else None))
            tables.append(dict(alias=alias,mode=mode,complete=complete,verified=sum(r['status']=='VERIFIED' for r in selected),
                scheduled=len(selected),identities=identities,online_ms=online_ms,cold_ms=cold_ms,cold_Ir=cost,
                S_Ir_per_sqrt_r=normalized,Ir_over_K_floor=floor,online_over_ic=online_ratio,
                cold_Ir_over_ic=ir_ratio,cold_time_over_ic=cold_ratio,
                actual_B_columns_rank={k:[list(v) for v in sorted(values,key=repr)] for k,values in sorted(shapes.items())},
                stage_diagnostics=stage_diagnostics(selected) if mode == 'ic' else None,
                comparison=comp))
        by_name={v['alias']:v for v in tables}
        # Ratios of equal-cell geometric means preserve the matched point law.
        for row in tables:
            online_ref=by_name.get('rho_online',{}).get('online_over_ic')
            cold_ref=by_name.get('rho',{}).get('cold_Ir_over_ic')
            row['rho_online_over_IC_online'] = online_ref/row['online_over_ic'] if online_ref and row['online_over_ic'] else None
            row['cold_Ir_over_cold_rho'] = row['cold_Ir_over_ic']/cold_ref if cold_ref and row['cold_Ir_over_ic'] else None
        stages[stage]=dict(table=tables,failures=summary['failures'],paired_online=summary.get('single_target_online'),
            retained_portfolio=summary.get('retained_portfolio'),provisional_challenger=summary.get('provisional_challenger'),
            runs=summary['runs'],verified_runs=summary['verified_runs'])
    keys={(r['identity'],r['workload_id'],r['run_id']) for r in runs}
    if len(keys)!=len(runs):raise ValueError('duplicate canonical run key')
    files=['contract.json','fixtures.json','candidates.json','decision.json']+[f'summaries/{stage}.json' for stage in c['stages']]
    report=dict(schema_version=1,round=c['attempt_number'],scope=c['evidence_scope'],decision=decision,stages=stages,
        primary_metric='single-target native online wall time; median of three processes per point, equal-cell geometric mean',
        cold_unit='complete user-space guest Ir; supplementary complete cold native process wall time',
        floor='K guest instructions for this full-rank collector only; inapplicable to rho',
        uncertainty='Final selected-candidate intervals use the frozen nominal familywise rule. Development and rho comparisons are descriptive.',
        input_sha256={name:sha(root/name) for name in files},total_runs=len(runs),verified_runs=sum(r['status']=='VERIFIED' for r in runs))
    return report,runs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--round', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    # Tables are not an independent evidence audit. Always run the actual sealed
    # evaluator first; do not let an edited summary become a published result.
    subprocess.run([sys.executable, str(a.round/'evaluator/tournament.py'),
                    'verify', '--round', str(a.round)], check=True)
    report, runs = export(a.round)
    if a.verify:
        if read(a.out/'RESULTS.json') != report:
            raise ValueError('exported table differs from audited receipts')
        import io
        expected = io.StringIO(newline='')
        writer = csv.DictWriter(expected, fieldnames=list(runs[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(runs)
        if (a.out/'RUNS.csv').read_bytes() != expected.getvalue().encode():
            raise ValueError('exported run rows differ from audited receipts')
    else:
        a.out.mkdir(parents=True, exist_ok=False)
        (a.out/'RESULTS.json').write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
        with (a.out/'RUNS.csv').open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(runs[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(runs)
    print(json.dumps(dict(status='VERIFIED' if a.verify else 'EXPORTED',
                          rows=len(runs), verified=report['verified_runs'])))

if __name__=='__main__':main()
