#!/usr/bin/env python3
"""Verify completeness and summarize the preregistered bounded experiment."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
import itertools
from functools import lru_cache
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


@lru_cache(maxsize=1024)
def source_support(active, polys):
    """Independent exact support/count reconstruction for the cache contract."""
    bits=[1<<i for i in range(active.bit_length()) if active&(1<<i)]
    schedules={}
    support=set();count=0
    for poly in polys:
        if not poly:continue
        degree=max(m.bit_count() for m in poly)
        if degree>3:continue
        gap=3-degree
        if gap not in schedules:
            schedules[gap]=[sum(c) for k in range(gap+1) for c in itertools.combinations(bits,k)]
        for multiplier in schedules[gap]:
            row=set()
            for m in poly:
                value=m|multiplier
                if value in row:row.remove(value)
                else:row.add(value)
            if row:count+=1;support.update(row)
    return count,frozenset(support)


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
        n,family,batch=fixture['n'],fixture['family'],fixture['batch']
        inputs=fixture['inputs'];tails=fixture['reference_tails']
        assert fixture['degree']==3 and len(inputs)==len(tails)==batch
        assert len({json.dumps(x) for x in inputs})==batch
        expected_rows=expected_columns=expected_hits=0;layout_cache={}
        allowed=sum(1<<(i*n//8) for i in range(8))
        for i,(entry,tail) in enumerate(zip(inputs,tails)):
            polys=entry['polys'];active=entry['active']
            assert len(polys)==protocol['fixture_generators']
            union=0
            for j,poly in enumerate(polys):
                assert poly==sorted(set(poly)) and len(poly)<=protocol['limits']['terms_per_generator']
                assert all(0<=m<(1<<n) for m in poly)
                for m in poly:union|=m
                degree=max((m.bit_count() for m in poly),default=-1);want=2
                if j==0 and family in ('linear_drop','restricted_cycle'):
                    if i%4==1:want=1
                    if i%4==3:want=-1
                    if i%4==2 and family=='restricted_cycle':want=0
                assert degree==want
                if family=='restricted_cycle':assert all(m&~allowed==0 for m in poly)
            assert active==union
            if family=='cross_cancel':
                difference=set(polys[0])^set(polys[1]);want={1<<(i%n)}
                if i%2:want.add(0)
                assert difference==want
            pivots=[]
            for row in tail:
                assert 0<row<(1<<(n+1))
                pivot=(row&-row).bit_length()-1;pivots.append(pivot)
                assert sum(bool(other&(1<<pivot)) for other in tail)==1
            assert pivots==sorted(set(pivots))
            if family=='cross_cancel':
                value=(1<<(i%n))|((1<<n) if i%2 else 0)
                for row,pivot in zip(tail,pivots):
                    if value&(1<<pivot):value^=row
                assert value==0 and tail
            count,columns=source_support(active,tuple(tuple(p) for p in polys))
            expected_rows+=count;expected_columns+=len(columns)
            key=(n,3,active)
            if key in layout_cache and layout_cache[key]==columns:expected_hits+=1
            layout_cache[key]=columns
        names=protocol['variants'];samples=rows[1:]
        assert len(samples)==len(names)*protocol['repetitions']
        pairs={(s['rep'],s['variant']):s for s in samples}
        assert len(pairs)==len(samples) and set(pairs)=={(rep,arm) for rep in range(protocol['repetitions']) for arm in names}
        for rep in range(protocol['repetitions']):
            selected=[s for s in samples if s['rep']==rep]
            assert [s['order'] for s in selected]==list(range(len(names)))
            assert [s['variant'] for s in selected]==[names[(rep+i)%len(names)] for i in range(len(names))]
        for sample in samples:
            arm=sample['variant']
            assert sample['verified_outputs']==batch
            assert sample['total_ns']>=sample['setup_ns']+sample['apply_ns']+sample['validation_ns']
            assert sample['source_rows']==expected_rows and sample['source_columns']==expected_columns
            assert sample['tail_rank']==sum(map(len,tails)) and sample['nonempty_tails']==sum(bool(t) for t in tails)
            assert sample['high_rank']+sample['tail_rank']<=min(expected_rows,expected_columns)
            assert sample['layout_hits']==(expected_hits if arm=='flat_cached' else 0)
            if arm=='flat_cached':
                assert sample['merge_items'] is None and sample['max_sparse_terms'] is None
                assert sample['high_word_xors']>=sample['high_xors']
            else:
                assert sample['high_word_xors'] is None
                assert sample['merge_items']>=0 and sample['max_sparse_terms']>=0
                assert sample['retained_layout_columns']==0
            if arm!='stream_exchange':assert sample['pivot_exchanges']==0
            by_group[fixture['split'],n,family,batch,arm].append(sample)
        assert len({s['high_rank'] for s in samples})==1
        arms={}
        for arm in names:
            values=[s for s in samples if s['variant']==arm]
            fields=['setup_ns','apply_ns','validation_ns','total_ns','high_xors','high_word_xors','merge_items','max_sparse_terms','pivot_exchanges','retained_schedule_masks','retained_layout_columns']
            arms[arm]={f'median_{f}':None if values[0][f] is None else statistics.median(s[f] for s in values) for f in fields}
            arms[arm]['layout_hits']=values[0]['layout_hits']
        summaries.append({'cell':cell,'n':n,'family':family,'batch':batch,'split':fixture['split'],'seed':fixture['seed'],
                          'source_rows':expected_rows,'source_columns':expected_columns,'high_rank':samples[0]['high_rank'],
                          'tail_rank':samples[0]['tail_rank'],'nonempty_tails':samples[0]['nonempty_tails'],'arms':arms,
                          'peak_rss_bytes':receipt['peak_rss_bytes'],'fixture_sha256':hashlib.sha256(fixture['raw_line'].encode()).hexdigest()})
    gates=[];candidates=['stream_plain','stream_exchange']
    for n in protocol['variables']:
        if n==12:continue
        for family in protocol['families']:
            cells=[c for c in summaries if c['n']==n and c['family']==family and c['batch']==8 and c['split']=='holdout']
            for candidate in candidates:
                ratios=[]
                for cell in cells:
                    pairs={(s['rep'],s['variant']):s for s in raw[cell['cell']][1:]}
                    ratios.extend(min(pairs[rep,arm]['total_ns'] for arm in ('flat_cached','sparse_bucket'))/pairs[rep,candidate]['total_ns']
                                  for rep in range(protocol['repetitions']))
                ci=bootstrap(ratios)
                gates.append({'n':n,'family':family,'batch':8,'candidate':candidate,'control':'pointwise faster of retained flat and sparse-bucket controls',
                              'paired_ratio_median':statistics.median(ratios),'ci95_paired_median':ci,'threshold':2.0,'pass':ci[0]>2.0})
    result={'schema_version':1,'scope':protocol['scope'],'correctness':'PASS','cells':len(summaries),
            'paired_batch_samples':sum(len(rows)-1 for rows in raw.values()),
            'verified_outputs':sum(s['verified_outputs'] for rows in raw.values() for s in rows[1:]),
            'dramatic_gate':{c:'PASS' if all(g['pass'] for g in gates if g['candidate']==c) else 'REJECTED' for c in candidates},
            'gates_passed':{c:sum(g['pass'] for g in gates if g['candidate']==c) for c in candidates},'gate_details':gates,'cells_detail':summaries,
            'classification':protocol['classification'],'memory_scope':protocol['memory'],'counter_scope':protocol['counters'],
            'full_polynomial_solver_cost':None,'full_ic_cost':None,'rho_ratio':None,
            'ci_scope':'paired bootstrap over repeated measurements and two fixed holdout seeds; not a population claim or independent reproduction'}
    (root/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    lines=['# Specialized Boolean linear-tail workload','',
           f"Correctness **PASS**: {result['cells']} cells, {result['paired_batch_samples']} batch-arm samples, {result['verified_outputs']} oracle-verified complete affine-tail outputs.",'',
           f"Dramatic gate: **{result['dramatic_gate']}**. Gates passed: **{result['gates_passed']}**, out of twelve per candidate.",'',
           f"Cold batch8 milliseconds include all setup, product generation, cache guards/misses, high elimination, affine canonicalization, validation and destruction. Values are medians over two holdout seeds and {protocol['repetitions']} balanced repetitions. Ratios below use pooled medians against the cached flat control; acceptance uses paired minima over both retained controls.",'',
           '| Variables | Family | Variant | Cold batch (ms) | Flat / arm | Layout hits / 8 | Nonempty outputs / 8 | Max sparse terms |',
           '|---:|---|---|---:|---:|---:|---:|---:|']
    for n in protocol['variables']:
        for family in protocol['families']:
            base=statistics.median(s['total_ns'] for s in by_group['holdout',n,family,8,'flat_cached'])
            for arm in protocol['variants']:
                values=by_group['holdout',n,family,8,arm];total=statistics.median(s['total_ns'] for s in values)
                fill='null' if values[0]['max_sparse_terms'] is None else f"{statistics.median(s['max_sparse_terms'] for s in values):.0f}"
                lines.append(f"| {n} | {family} | {arm} | {total/1e6:.6f} | {base/total:.3f} | {statistics.median(s['layout_hits'] for s in values):.0f} | {statistics.median(s['nonempty_tails'] for s in values):.0f} | {fill} |")
    lines+=['','The active mask is recomputed from each input. Restricted-cycle fixtures use only eight embedded variables; cross-cancel fixtures guarantee an affine consequence from two nonlinear generators. Exact source dimensions, high rank and canonical tail are checked. Empty and nonempty tail workloads are both retained.','',
            'Schedule/layout retention is reported as counts, not byte estimates. Whole-worker RSS includes all arms and references. Packed XORs and sparse merge items are different diagnostics, not a common complete operation unit. Independent oracle preparation is outside arm timing and inside process receipts.','',
            'The controls retain pre-existing kernel strategies; no novelty claim attaches to sparse high-column elimination or linear-tail restriction. These standalone measurements do not execute the production solver, enumerate roots, recover scalars or measure full index calculus. Those costs remain null.']
    (root/'RESULT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ['correctness','cells','paired_batch_samples','verified_outputs','dramatic_gate','gates_passed']},sort_keys=True))


if __name__=='__main__':
    main(Path(sys.argv[1]).resolve())
