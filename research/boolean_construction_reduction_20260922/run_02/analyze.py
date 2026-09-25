#!/usr/bin/env python3
"""Verify completeness and summarize the preregistered bounded experiment."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import sys


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bootstrap(values):
    rng = random.Random(20260922)
    medians = sorted(statistics.median(rng.choices(values, k=len(values))) for _ in range(4000))
    return [medians[100], medians[3899]]


def main(root: Path):
    protocol, metadata = read(root / "protocol.json"), read(root / "metadata.json")
    assert metadata["complete"], "incomplete campaign cannot be promoted"
    for name, expected in metadata["source_hashes"].items():
        assert sha(root / name) == expected, f"source mismatch: {name}"
    test = read(root / "test_receipt.json")
    assert test["exit_code"] == 0 and not test["timed_out"]
    for stream in ("stdout", "stderr"):
        assert hashlib.sha256(test[stream].encode()).hexdigest() == test[f"{stream}_sha256"]
    expected_cells = {
        f"n{n}-{split}-{seed}-{family}-b{batch}"
        for n in protocol["variables"]
        for split in ("discovery", "holdout")
        for seed in protocol[f"{split}_seeds"]
        for family in protocol["families"]
        for batch in protocol["batches"]
    }
    receipts = read(root / "receipts.json")
    assert len(receipts) == len(expected_cells)
    receipts = {r["cell"]: r for r in receipts}
    assert set(receipts) == expected_cells
    raw = defaultdict(list)
    for n in protocol["variables"]:
        for line in (root / f"raw-n{n}.jsonl").read_text().splitlines():
            row = json.loads(line)
            decoded = json.loads(row["raw_line"])
            assert decoded == {k: v for k, v in row.items() if k not in ("cell", "split", "raw_line")}
            raw[row["cell"]].append(row)
    assert set(raw) == expected_cells
    summaries, by_group = [], defaultdict(list)
    for cell, rows in sorted(raw.items()):
        receipt = receipts[cell]
        assert receipt["exit_code"] == 0 and not receipt["timed_out"]
        reconstructed = "".join(r["raw_line"] for r in rows).encode()
        assert hashlib.sha256(reconstructed).hexdigest() == receipt["stdout_sha256"], cell
        assert receipt["stderr_sha256"] == hashlib.sha256(b"").hexdigest(), cell
        assert rows[0]["type"] == "fixture"
        fixture = rows[0]
        for field in ("n", "seed", "family", "batch", "split"):
            assert fixture[field] == receipt[field]
        n, family, batch = fixture["n"], fixture["family"], fixture["batch"]
        inputs = fixture["inputs"]
        expected_active = sum(1 << (i*n//8) for i in range(8)) if family == "restricted_cycle" else (1 << n)-1
        assert fixture["degree"] == 3 and fixture["active"] == expected_active
        assert len(inputs) == batch and len({json.dumps(x) for x in inputs}) == batch
        for i, polys in enumerate(inputs):
            assert len(polys) == protocol["fixture_generators"]
            for j, poly in enumerate(polys):
                assert poly == sorted(set(poly)) and len(poly) <= protocol["limits"]["terms_per_generator"]
                assert all(0 <= term < (1 << n) for term in poly)
                d = max((m.bit_count() for m in poly), default=-1)
                expected_degree = 2
                if j == 0 and family != "quadratic":
                    if i % 4 == 1: expected_degree = 1
                    if i % 4 == 3: expected_degree = -1
                    if family == "restricted_cycle" and i % 4 == 2: expected_degree = 0
                assert d == expected_degree
        names=protocol["bridge_variants"] if n==protocol["bridge_n"] else protocol["variants"]
        samples=rows[1:]
        assert len(samples)==protocol["repetitions"]*len(names)
        pairs={(s['rep'],s['variant']):s for s in samples}
        assert len(pairs)==len(samples)
        assert set(pairs)=={(rep,arm) for rep in range(protocol['repetitions']) for arm in names}
        for rep in range(protocol['repetitions']):
            selected=[s for s in samples if s['rep']==rep]
            assert [s['order'] for s in selected]==list(range(len(names)))
            assert [s['variant'] for s in selected]==[names[(rep+i)%len(names)] for i in range(len(names))]
        for sample in samples:
            arm=sample['variant']
            assert sample['verified_outputs']==batch
            if arm=='stream_reduce':
                assert sample['construction_ns'] is None and sample['reduction_ns'] is None
                assert isinstance(sample['fused_ns'],int) and sample['fused_ns']>=0
            else:
                assert sample['fused_ns'] is None
                assert all(isinstance(sample[f],int) and sample[f]>=0 for f in ['construction_ns','reduction_ns'])
            parts=['setup_ns','construction_ns','reduction_ns','fused_ns','validation_ns']
            assert sample['total_ns']>=sum(sample[f] for f in parts if sample[f] is not None)
            assert 0<=sample['total_rank']<=sample['source_rows']<=batch*protocol['limits']['rows']
            assert sample['total_columns']<=batch*protocol['limits']['columns']
            assert 0<=sample['row_xors']<=sample['word_xors']
            if 'auxiliary_max_bytes' in sample:
                assert sample['auxiliary_max_bytes'] >= 0
                if arm != 'incidence_reduce': assert sample['auxiliary_max_bytes'] == 0
            by_group[fixture['split'],n,family,batch,arm].append(sample)
        for field in ['output_bytes','source_rows','total_rank','total_columns','row_xors']:
            assert len({s[field] for s in samples})==1, (cell,field)
        arms={}
        for arm in names:
            values=[s for s in samples if s['variant']==arm]
            fields=['setup_ns','construction_ns','reduction_ns','fused_ns','validation_ns','total_ns','retained_bytes','basis_max_bytes','word_xors']
            arms[arm]={f'median_{f}':None if values[0][f] is None else statistics.median(s[f] for s in values) for f in fields}
            if 'auxiliary_max_bytes' in values[0]:
                arms[arm]['median_auxiliary_max_bytes']=statistics.median(s['auxiliary_max_bytes'] for s in values)
            if arm!='stream_reduce':
                arms[arm]['median_construction_fraction']=statistics.median(s['construction_ns']/s['total_ns'] for s in values)
                arms[arm]['median_reduction_fraction']=statistics.median(s['reduction_ns']/s['total_ns'] for s in values)
                arms[arm]['construction_zero_time_projection']=statistics.median(s['total_ns']/(s['total_ns']-s['construction_ns']) for s in values)
        summaries.append({'cell':cell,'n':n,'family':family,'batch':batch,'seed':fixture['seed'],'split':fixture['split'],
                          'arms':arms,'source_rows':samples[0]['source_rows'],'total_rank':samples[0]['total_rank'],
                          'total_columns':samples[0]['total_columns'],'output_bytes':samples[0]['output_bytes'],
                          'row_xors':samples[0]['row_xors'],'peak_rss_bytes':receipt['peak_rss_bytes'],
                          'fixture_sha256':hashlib.sha256(fixture['raw_line'].encode()).hexdigest()})
    candidate_name=protocol.get("candidate_arm","stream_reduce")
    controls=protocol.get("comparison_controls",["sorted_reduce","ranked_reduce","sparse_reduce"])
    reference_label=protocol.get("reference_label","pointwise fastest staged arm")
    gates,comparisons=[],[]
    for n in protocol['variables']:
        for family in protocol['families']:
            cells=[c for c in summaries if c['split']=='holdout' and c['n']==n and c['family']==family and c['batch']==8]
            paired=[{(s['rep'],s['variant']):s for s in raw[c['cell']][1:]} for c in cells]
            names=protocol['bridge_variants'] if n==protocol['bridge_n'] else protocol['variants']
            for candidate in names[1:]:
                ratios=[p[rep,'sorted_reduce']['total_ns']/p[rep,candidate]['total_ns'] for p in paired for rep in range(protocol['repetitions'])]
                comparisons.append({'n':n,'family':family,'candidate':candidate,'control':'sorted_reduce',
                                    'paired_ratio_median':statistics.median(ratios),'ci95_paired_median':bootstrap(ratios)})
            if n!=protocol['bridge_n']:
                ratios=[min(p[rep,c]['total_ns'] for c in controls)/p[rep,candidate_name]['total_ns']
                        for p in paired for rep in range(protocol['repetitions'])]
                ci=bootstrap(ratios)
                gates.append({'n':n,'family':family,'batch':8,'control':reference_label,'candidate':candidate_name,
                              'paired_ratio_median':statistics.median(ratios),'ci95_paired_median':ci,'threshold':2.0,'pass':ci[0]>2.0})
    result={'schema_version':1,'scope':protocol['scope'],'correctness':'PASS','cells':len(summaries),
            'paired_batch_samples':sum(len(rows)-1 for rows in raw.values()),
            'verified_outputs':sum(s['verified_outputs'] for rows in raw.values() for s in rows[1:]),
            'dramatic_combined_gate':'PASS' if all(g['pass'] for g in gates) else 'REJECTED','gates_passed':sum(g['pass'] for g in gates),
            'gate_details':gates,'combined_comparisons':comparisons,'cells_detail':summaries,
            'dense_not_executed':[{'n':n,'status':'NOT_EXECUTED','total_ns':None} for n in protocol['variables'] if n>protocol['bridge_n']],
            'classification':protocol['classification'],'memory_scope':protocol['memory'],'counter_scope':protocol['operation_counters'],
            'projection_scope':protocol['projection'],'full_polynomial_solver_cost':None,'rho_ratio':None,
            'ci_scope':'paired bootstrap on repeated measurements and two fixed holdout seeds; not a population claim or independent reproduction'}
    if protocol.get('enable_incidence'): result['candidate_arm']=candidate_name
    (root/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    acceptance_reference = 'the fastest prior arm, including streaming' if protocol.get('enable_incidence') else 'the fastest staged arm'
    lines=['# Boolean construction plus canonical linear reduction','',
           f"Correctness **PASS**: {result['cells']} cells, {result['paired_batch_samples']} batch-arm samples, {result['verified_outputs']} oracle-verified RREF outputs.",'',
           f"Dramatic combined-workload gate: **{result['dramatic_combined_gate']}**, {result['gates_passed']} / {len(gates)} comparisons passed.",'',
           f"Cold batch8 milliseconds include setup, construction, complete forward/backward elimination, compaction, exact validation and destruction. Medians pool two holdout seeds and {protocol['repetitions']} balanced repetitions. Ratios below use pooled medians; acceptance uses paired ratios against {acceptance_reference}.",'',
           '| Variables | Family | Variant | Cold batch (ms) | Sorted / arm | Construction (ms) | Reduction (ms) | Fused (ms) | Reduction word XORs |',
           '|---:|---|---|---:|---:|---:|---:|---:|---:|']
    for n in protocol['variables']:
        names=protocol['bridge_variants'] if n==protocol['bridge_n'] else protocol['variants']
        for family in protocol['families']:
            baseline=statistics.median(s['total_ns'] for s in by_group['holdout',n,family,8,'sorted_reduce'])
            for arm in names:
                samples=by_group['holdout',n,family,8,arm]
                total=statistics.median(s['total_ns'] for s in samples)
                phase=lambda key:'null' if samples[0][key] is None else f"{statistics.median(s[key] for s in samples)/1e6:.6f}"
                lines.append(f"| {n} | {family} | {arm} | {total/1e6:.6f} | {baseline/total:.3f} | {phase('construction_ns')} | {phase('reduction_ns')} | {phase('fused_ns')} | {statistics.median(s['word_xors'] for s in samples):.0f} |")
    lines+=['','Fused construction/reduction phases are not separately observable and remain null, not zero. Source rows and logical row-XOR counts agree across all arms; physical word-XOR counts may differ with ambient coordinate width. These counters cover reduction only, not all calibrated operations.','',
            'Retained context bytes and end-of-insertion basis storage are recorded separately. Neither is whole-process peak memory. Worker RSS includes all arms and the common reference corpus. Fixture and independent oracle construction are outside arm timing and inside process receipts.','',
            'This completes the bounded construction-plus-RREF task. It does not compute complete Groebner closure, enumerate roots, solve the original polynomial system, or measure index-calculus performance. Those costs and rho ratios remain null.']
    if protocol.get('enable_incidence'):
        lines+=['','This additive run evaluates incidence_reduce: sparse construction, cached nonzero pivot words and precomputed backward-elimination incidence. All preparation and auxiliary storage are charged; exact canonical output and logical row-XOR counts still match. The primary gate compares it against the pointwise fastest prior pipeline, including the streaming control.']
    (root/'RESULT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ['correctness','cells','paired_batch_samples','verified_outputs','dramatic_combined_gate','gates_passed']},sort_keys=True))


if __name__=='__main__':
    main(Path(sys.argv[1]).resolve())
