#!/usr/bin/env python3
"""Verify exact generated-system evidence and complete cold comparison gates."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import statistics
import sys
sys.dont_write_bytecode = True


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def satisfies(polys, point):
    return all(sum(point & m == m for m in poly) % 2 == 0 for poly in polys)


def draws(state):
    mask = (1 << 64) - 1
    while True:
        state = (state + 0x9e3779b97f4a7c15) & mask
        z = state
        z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & mask
        z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & mask
        yield z ^ (z >> 31)


def fixture(n, seed, family):
    rng = draws(seed)
    witness = next(rng) & ((1 << n) - 1)
    polys = []
    for _ in range(n + 2):
        terms = set()
        while len(terms) < 12:
            a, b = next(rng) % n, next(rng) % n
            if a != b:
                terms.add((1 << a) | (1 << b))
        for _ in range(4):
            terms.add(1 << (next(rng) % n))
        constant = next(rng) & 1 if family == 'unplanted' else sum(witness & m == m for m in terms) % 2
        if constant:
            terms.add(0)
        polys.append(sorted(terms, key=lambda m: (-m.bit_count(), m)))
    if family == 'cross_planted':
        bit = 1 << (seed % n)
        terms = set(polys[0]) ^ {bit}
        if witness & bit:
            terms ^= {0}
        polys[1] = sorted(terms, key=lambda m: (-m.bit_count(), m))
    return polys, None if family == 'unplanted' else witness


def paired_order(count, seed, rep):
    rng = draws(seed ^ (((rep // 2) * 0xd1b54a32d192ed03) & ((1 << 64) - 1)))
    order = list(range(count))
    for i in range(count - 1, 0, -1):
        j = next(rng) % (i + 1)
        order[i], order[j] = order[j], order[i]
    return order[::-1] if rep % 2 else order


def quotient_rank(polys, k):
    # Independent row-space rank implementation, without the worker's projection.
    basis = {}
    for j in range(k):
        for i in range(j):
            m = (1 << i) | (1 << j)
            value = sum((poly.count(m) % 2) << e for e, poly in enumerate(polys))
            while value:
                p = value.bit_length() - 1
                if p not in basis:
                    basis[p] = value
                    break
                value ^= basis[p]
    return len(basis)


def interval(values):
    rng = random.Random(20260923)
    medians = sorted(statistics.median(rng.choices(values, k=len(values))) for _ in range(4000))
    return [medians[100], medians[3899]]


def check_work(work, polys, n, k, outcome):
    assert work['low_variables'] == k
    assert all(isinstance(v, int) and v >= 0 for key, v in work.items() if key != 'consistent_rank_counts')
    assert work['annihilated_rank'] == quotient_rank(polys, k)
    assert work['quotient_dimension'] == len(polys) - work['annihilated_rank']
    assert work['prefixes'] == 16 * work['batches'] <= 1 << (n - k)
    assert work['affine_queries'] + work['screen_rejected'] <= work['prefixes']
    assert work['affine_rejected'] <= work['affine_queries']
    counts = work['consistent_rank_counts']
    assert len(counts) == 7 and all(isinstance(v, int) and v >= 0 for v in counts)
    assert all(v == 0 for v in counts[k+1:])
    assert sum(counts) == work['affine_queries'] - work['affine_rejected']
    assert work['extension_space'] == sum(v * (1 << (k-r)) for r, v in enumerate(counts[:k+1]))
    assert work['extensions_checked'] <= work['extension_space'] <= 1 << n
    assert work['original_rejected'] <= work['extensions_checked']
    assert sum(r*v for r, v in enumerate(counts)) <= work['rank_sum'] <= k * work['affine_queries']
    if outcome == 'UNSAT':
        assert work['prefixes'] == 1 << (n-k)
        assert work['affine_queries'] + work['screen_rejected'] == work['prefixes']
        assert work['extensions_checked'] == work['extension_space'] == work['original_rejected']
    elif outcome == 'SAT':
        assert work['extensions_checked'] == work['original_rejected'] + 1


def validate(root):
    protocol, meta = read(root/'protocol.json'), read(root/'metadata.json')
    assert meta['complete']
    for name, digest in meta['source_hashes'].items():
        assert sha(root/name) == digest, name
    old = read(root/'REFERENCE_PROTOCOL.json')
    assert protocol['retained_arms'] == protocol['reference_arms'] == old['variants']
    assert protocol['new_candidates'] == ['projected4', 'projected5', 'projected6']
    assert protocol['variants'] == old['variants'] + protocol['new_candidates']
    assert protocol['variables'] == [12,16,20,24] and protocol['repetitions'] == 8
    assert protocol['families'] == ['planted','cross_planted','unplanted']
    assert protocol['splits'] == ['discovery','regression','holdout']
    assert protocol['discovery_seeds'] == old['discovery_seeds']
    assert protocol['regression_seeds'] == old['regression_seeds'] + old['holdout_seeds']
    seeds = sum((protocol[s+'_seeds'] for s in protocol['splits']), [])
    assert len(seeds) == len(set(seeds)) == 18
    assert protocol['holdout_seeds'] == [20261024,3145729]
    lineage = read(root/'SOURCE_LINEAGE.json')
    for name, digest in lineage['files'].items():
        data = (root/name).read_bytes()
        if name == 'worker.rs':
            suffix = b'\ninclude!("projected.rs");\ninclude!("projected_main.rs");\n'
            assert data.endswith(suffix)
            data = data[:-len(suffix)]
        assert hashlib.sha256(data).hexdigest() == digest, name
    test = read(root/'test_receipt.json')
    assert test['exit_code'] == 0 and not test['timed_out']
    for stream in ['stdout','stderr']:
        assert hashlib.sha256(test[stream].encode()).hexdigest() == test[stream+'_sha256']
    expected = {(n,s,seed,f): f'n{n}-{s}-{seed}-{f}' for n in protocol['variables'] for s in protocol['splits'] for seed in protocol[s+'_seeds'] for f in protocol['families']}
    receipts = read(root/'receipts.json')
    assert len(receipts) == len(expected) == 216
    receipts = {r['cell']:r for r in receipts}
    assert set(receipts) == set(expected.values())
    previous = {(c['n'],c['seed'],c['family']):c for c in read(root/'REFERENCE_FIXTURES.json')}
    assert len(previous) == 192
    groups = defaultdict(list)
    summaries = []
    total_samples = 0
    all_complete = True
    for (n,split,seed,family),cell in expected.items():
        receipt = receipts[cell]
        assert receipt['exit_code'] == 0 and not receipt['timed_out']
        assert (receipt['n'],receipt['split'],receipt['seed'],receipt['family']) == (n,split,seed,family)
        for suffix, stream in [('.jsonl','stdout'),('.stderr','stderr')]:
            assert sha(root/(cell+suffix)) == receipt[stream+'_sha256']
        assert (root/(cell+'.stderr')).read_bytes() == b''
        order_seed = (protocol['order_seed'] ^ (n<<48) ^ (seed<<8) ^ protocol['families'].index(family)) & ((1<<64)-1)
        assert receipt['order_seed'] == order_seed
        assert receipt['command'][1:] == [str(n),str(seed),family,'8','200000',str(order_seed)]
        rows = [json.loads(line) for line in (root/(cell+'.jsonl')).read_text().splitlines()]
        source, samples = rows[0], rows[1:]
        assert source['type'] == 'fixture' and (source['n'],source['seed'],source['family']) == (n,seed,family)
        polys,witness = fixture(n,seed,family)
        assert source['polys'] == polys and source['planted_witness'] == witness
        if split != 'holdout':
            assert previous[n,seed,family] == {key:source[key] for key in ['n','seed','family','polys','planted_witness']}
        if witness is not None:
            assert satisfies(polys,witness)
        names = protocol['variants']
        assert len(samples) == 8*len(names)
        pairs = {(s['rep'],s['variant']):s for s in samples}
        assert len(pairs) == len(samples)
        stable = {}
        for rep in range(8):
            chunk = samples[rep*len(names):(rep+1)*len(names)]
            assert [s['variant'] for s in chunk] == [names[i] for i in paired_order(len(names),order_seed,rep)]
            assert [s['order'] for s in chunk] == list(range(len(names)))
            assert all(s['rep'] == rep and s['type'] == 'sample' for s in chunk)
            for sample in chunk:
                arm, outcome = sample['variant'],sample['outcome']
                assert all(isinstance(sample[t],int) and sample[t]>=0 for t in ['solve_ns','validation_ns','total_ns'])
                assert 0 < sample['solve_ns'] + sample['validation_ns'] <= sample['total_ns']
                assert outcome in ['SAT','UNSAT','UNKNOWN']
                if outcome == 'SAT':
                    assert isinstance(sample['model'],int) and 0 <= sample['model'] < 1<<n
                    assert satisfies(polys,sample['model']) and sample['verified'] and sample['reason'] is None
                    assert source['search_reference'] != 'UNSAT'
                elif outcome == 'UNSAT':
                    assert sample['model'] is None and sample['reason'] is None and sample['verified']
                    assert source['search_reference'] == 'UNSAT' and witness is None
                else:
                    assert sample['model'] is None and isinstance(sample['reason'],str) and not sample['verified']
                    all_complete = False
                if arm in protocol['new_candidates']:
                    check_work(sample['projected'],polys,n,int(arm[-1]),outcome)
                    assert sample['trace'] is None and not any(sample['logical'].values())
                else:
                    assert sample['projected'] is None
                signature = {k:v for k,v in sample.items() if k not in ['rep','order','solve_ns','validation_ns','total_ns']}
                if arm in stable:
                    assert stable[arm] == signature
                stable[arm] = signature
                groups[split,n,family,arm].append(sample['total_ns'])
            if source['search_reference'] in ['SAT','UNSAT']:
                assert pairs[rep,'search']['outcome'] == source['search_reference']
        total_samples += len(samples)
        summaries.append({'cell':cell,'n':n,'split':split,'seed':seed,'family':family,'status':source['search_reference'],
            'complete':all(s['outcome']!='UNKNOWN' for s in samples),
            'medians_ns':{arm:statistics.median(pairs[r,arm]['total_ns'] for r in range(8)) for arm in names},
            'projected_work':{a:pairs[0,a]['projected'] for a in protocol['new_candidates']}})
    gates = []
    for candidate in protocol['new_candidates']:
        for split in ['regression','holdout']:
            for n in [16,20,24]:
                for family in protocol['families']:
                    values = groups[split,n,family,candidate]
                    best = [min(groups[split,n,family,a][j] for a in protocol['reference_arms']) for j in range(len(values))]
                    ratios = [a/b for a,b in zip(best,values)]
                    ci = interval(ratios)
                    gates.append({'candidate':candidate,'split':split,'n':n,'family':family,'observations':len(ratios),
                        'median':statistics.median(ratios),'ci95':ci,'dramatic_pass':all_complete and ci[0]>2.0,'incremental_pass':all_complete and ci[0]>1.0})
    decisions = {c:{'dramatic_pass':all(g['dramatic_pass'] for g in gates if g['candidate']==c),
                    'incremental_pass':all(g['incremental_pass'] for g in gates if g['candidate']==c),
                    'dramatic_groups':sum(g['dramatic_pass'] for g in gates if g['candidate']==c),
                    'incremental_groups':sum(g['incremental_pass'] for g in gates if g['candidate']==c)} for c in protocol['new_candidates']}
    table = [{'variant':a,**{f:statistics.median(groups['holdout',24,f,a])/1e6 for f in protocol['families']}} for a in names]
    return {'schema_version':1,'all_complete':all_complete,'cells':len(summaries),'samples':total_samples,'summaries':summaries,
        'decisions':decisions,'gates':gates,'n24_holdout_ms':table,'campaign_seconds':meta['campaign_seconds'],
        'peak_worker_rss_bytes':max(r['peak_rss_bytes'] for r in receipts.values()),
        'production_solver_cost':None,'full_ic_cost':None,'rho_ratio':None,'calibrated_operation_ratio':None}


def main(root):
    result=validate(root)
    dump(root/'results.json',result)
    print(json.dumps({k:result[k] for k in ['all_complete','cells','samples','decisions','campaign_seconds','peak_worker_rss_bytes']}))


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve())
