"""Independent F81 incidence replay; imports no code from the benchmark.

Field multiplication uses explicit polynomial reduction; inverses use exhaustive
search. Curve sums use synthetic division of the line-intersection cubic.
This is a finite computation certificate checker, not a formal proof assistant.
"""
import argparse, collections, gzip, itertools, json, re
from pathlib import Path

def require(condition, message):
    if not condition:
        raise ValueError(message)

def verify(raw):
    digits = list(itertools.product(range(3), repeat=4))
    digits.sort(key=lambda v: sum(c * 3**i for i,c in enumerate(v)))
    def enc(v):
        return sum((c % 3) * 3**i for i,c in enumerate(v))
    add = [[enc([a+b for a,b in zip(x,y)]) for y in digits] for x in digits]
    neg = [enc([-a for a in x]) for x in digits]
    mul = []
    for x in digits:
        row = []
        for y in digits:
            c = [0]*7
            for i in range(4):
                for j in range(4):
                    c[i+j] += x[i]*y[j]
            # X^4 = X^3 - X - 1. Eliminate from highest power downward.
            for k in range(6,3,-1):
                c[k-1] += c[k]
                c[k-3] -= c[k]
                c[k-4] -= c[k]
                c[k] = 0
            row.append(enc(c[:4]))
        mul.append(row)
    inv = [None] + [next(b for b in range(1,81) if mul[a][b] == 1) for a in range(1,81)]
    def sub(a,b): return add[a][neg[b]]
    def quotient(a,b):
        require(b != 0, 'zero denominator')
        return mul[a][inv[b]]
    def frob(a): return mul[mul[a][a]][a]
    def histogram(values, action):
        unseen=set(values); counts=collections.Counter()
        while unseen:
            first=min(unseen); orbit={first}; value=action(first)
            while value != first:
                require(value in unseen and value not in orbit, 'invalid orbit')
                orbit.add(value); value=action(value)
            unseen-=orbit; counts[str(len(orbit))]+=1
        return dict(counts)
    linear={enc([a,b,0,0]) for a in range(3) for b in range(3)}
    fractions={quotient(a,b) for a in linear for b in linear if b}
    require(len(fractions)==27 and {frob(x) for x in fractions}==fractions,'fraction construction')
    cases=[raw['tiny']]+raw['other_five_curves']; require(len(cases)==6,'six curves required')
    seen=set(); rows=0; pair_checks=0
    for case in cases:
        match=re.fullmatch(r'y\^2=x\^3\+([12])x\+([012])',case['curve'])
        require(match is not None,'curve outside frozen family')
        a,b=map(int,match.groups()); require((a,b) not in seen,'duplicate curve'); seen.add((a,b))
        points=[None]+[(x,y) for x in range(81) for y in range(81)
                     if mul[y][y]==add[add[mul[mul[x][x]][x]][mul[a][x]]][b]]
        require([None if p is None else list(p) for p in points]==case['points'],'point enumeration differs')
        index={p:i for i,p in enumerate(points)}
        def divide_linear(poly,root):
            out=[0]*(len(poly)-1); out[-1]=poly[-1]
            for k in range(len(out)-2,-1,-1):
                out[k]=add[poly[k+1]][mul[root][out[k+1]]]
            require(add[poly[0]][mul[root][out[0]]]==0,'line does not intersect given point')
            return out
        def plus(p,q):
            if p is None:return q
            if q is None:return p
            x,y=p; u,v=q
            if x==u and add[y][v]==0:return None
            slope=quotient(sub(v,y),sub(u,x)) if x!=u else quotient(a,mul[2][y])
            intercept=sub(y,mul[slope][x])
            # Substitute line into curve; remove known roots by synthetic division.
            poly=[sub(b,mul[intercept][intercept]),sub(a,mul[2][mul[slope][intercept]]),neg[mul[slope][slope]],1]
            linear_factor=divide_linear(divide_linear(poly,x),u)
            third=quotient(neg[linear_factor[0]],linear_factor[1])
            return third,neg[add[mul[slope][third]][intercept]]
        table=[[index[plus(p,q)] for q in points] for p in points]
        expected_functionals={v for v in digits if any(v) and next(x for x in v if x)==1}
        require(len(case['all_40_hyperplanes'])==40,'missing hyperplane')
        actual_functionals={tuple(v['functional']) for v in case['all_40_hyperplanes']}
        require(actual_functionals==expected_functionals,'hyperplane identities differ')
        for v in [case['fraction']]+case['all_40_hyperplanes']:
            if v is case['fraction']: values=fractions
            else:
                c=v['functional']; values={i for i,x in enumerate(digits) if sum(a*b for a,b in zip(c,x))%3==0}
            base=[i for i,p in enumerate(points) if p is not None and p[0] in values]
            counts=[0]*len(points); distinct=[0]*len(points)
            for i in base:
                for j in base:
                    target=table[i][j]; counts[target]+=1; pair_checks+=1
                    if points[i][0]!=points[j][0]:distinct[target]+=1
            require(counts==v['counts_in_point_order'],'per-target pair counts differ')
            require(distinct==v['distinct_x_counts_in_point_order'],'distinct-x pair counts differ')
            require(len(base)==v['affine_base_points'],'base point count differs')
            require(sum(x>0 for x in counts[1:])==v['covered_affine_targets'],'coverage differs')
            invariant={frob(x) for x in values}==values
            require(invariant==v['frobenius_invariant'],'invariance differs')
            if invariant:
                require(histogram(values,frob)==v['field_orbits'],'field orbits differ')
                require(histogram(base,lambda i:index[(frob(points[i][0]),frob(points[i][1]))])==v['point_orbits'],'point orbits differ')
            rows+=1
    require(seen==set(itertools.product([1,2],range(3))),'curve family incomplete')
    return {'status':'verified','curves':len(seen),'factor_base_rows':rows,'ordered_point_pairs_recomputed':pair_checks,
            'method':'separate polynomial arithmetic; exhaustive inverse tables; synthetic division of line-intersection cubic',
            'scope':'F81 field sets, point enumeration, full pair-count vectors, distinct-x vectors and invariant orbits',
            'not_independently_checked':['F13^7 construction','formal proof of algorithms','full-DLP cost','solver speedups']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('raw');p.add_argument('--output',required=True);args=p.parse_args()
    text=gzip.open(args.raw,'rt').read() if args.raw.endswith('.gz') else Path(args.raw).read_text()
    result=verify(json.loads(text))
    with open(args.output,'x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))
