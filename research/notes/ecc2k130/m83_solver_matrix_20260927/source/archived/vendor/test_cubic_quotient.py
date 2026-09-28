#!/usr/bin/env python3
"""Reproducible exact orbit tests and supplied-orbit point-lifting experiments."""
import argparse
import collections
import json
import random
import statistics
import time
from pathlib import Path

from frobenius_cubic_quotient import FrobeniusQuotient,OrbitKey
from frobenius_curve_arithmetic import Field,Curve,binary_rank


def exhaustive(p):
    start=time.perf_counter();q=FrobeniusQuotient(p)
    representatives={};fibres=collections.Counter();inverse_count=0
    for word in range(1<<p):
        key=q.encode(word)
        if key in representatives:
            assert word in representatives[key]
        else:
            orbit={q.rotate(word,j) for j in range(p)}
            decoded=q.decode(key)
            assert decoded in orbit and q.encode(decoded)==key
            representatives[key]=orbit;inverse_count+=1
        fibres[key]+=1
    expected=2+((1<<p)-2)//p
    assert len(representatives)==expected
    assert dict(collections.Counter(fibres.values()))=={1:2,p:expected-2}
    return {'degree':p,'words':1<<p,'distinct_orbits':len(representatives),
            'decoded_orbits':inverse_count,'fibre_sizes':dict(collections.Counter(fibres.values())),
            'seconds':time.perf_counter()-start}


def derivative(q,base,directions):
    value=0
    for mask in range(1<<len(directions)):
        x=base
        for i,d in enumerate(directions):
            if (mask>>i)&1:x^=d
        value^=q.encode(x).value
    return value


def degree131_tests(count):
    q=FrobeniusQuotient(131);rng=random.Random(20260921)
    fixtures=[0,q.mask,1,2,1<<130,q.mask^1,q.mask^(1<<130)]
    fixtures += [rng.getrandbits(131) for _ in range(count)]
    # Structured small-space words and their complements are included.
    fixtures += [q.rotate(1|(rng.getrandbits(43)<<1),rng.randrange(131)) for _ in range(100)]
    fixtures += [q.mask^x for x in fixtures[-100:]]
    keys=[];enc=[];dec=[];canon=[]
    for word in fixtures:
        start=time.perf_counter();key=q.encode(word);enc.append(time.perf_counter()-start)
        start=time.perf_counter();decoded=q.decode(key);dec.append(time.perf_counter()-start)
        assert q.encode(decoded)==key and q.phase(word,decoded)>=0
        shift=rng.randrange(131)
        assert q.encode(q.rotate(word,shift))==key
        start=time.perf_counter();q.canonical_rotation(word);canon.append(time.perf_counter()-start)
        keys.append(key)
    degree3=derivative(q,0,[1,2,4])
    assert degree3!=0
    for _ in range(100):
        assert derivative(q,rng.getrandbits(131),[rng.getrandbits(131) for _ in range(4)])==0
    invalid=0
    for _ in range(100):
        value=rng.getrandbits(130)
        admissible=value==0 or q.field.pow(value,q.field.image_order)==1
        try:q.decode(OrbitKey(rng.randrange(2),value));assert admissible
        except ValueError:assert not admissible;invalid+=1
    # Identical input orbit keys can have different coordinate-addition results.
    assert q.encode(1)==q.encode(2)
    assert q.encode(1^1)!=q.encode(1^2)
    return {'degree':131,'random_seed':20260921,'tested_words':len(fixtures),
            'random_words':count,'fourth_derivative_tests':100,'nonzero_third_derivative':degree3,
            'random_key_validity_tests':100,'invalid_keys_rejected':invalid,
            'image_group_order':q.field.image_order,'root_exponent':q.field.root_exponent,
            'root_exponent_bitlength':q.field.root_exponent.bit_length(),
            'root_exponent_popcount':q.field.root_exponent.bit_count(),
            'median_encode_seconds':statistics.median(enc),
            'median_decode_validated_seconds':statistics.median(dec),
            'median_minimum_rotation_seconds':statistics.median(canon),
            'example':{'coordinate_word':fixtures[7], 'key':vars(keys[7]),
                       'decoded_word':q.decode(keys[7])}}


class NormalCoordinates:
    def __init__(self,field,seed):
        rng=random.Random(seed);self.f=field
        while True:
            beta=rng.randrange(1,1<<field.n);basis=[beta]
            for _ in range(1,field.n):basis.append(field.sq(basis[-1]))
            if binary_rank(basis)==field.n:break
        self.beta,self.basis=beta,basis;self.pivots={}
        for i,v in enumerate(basis):
            mask=1<<i
            while v:
                j=v.bit_length()-1
                if j in self.pivots:
                    vv,mm=self.pivots[j];v^=vv;mask^=mm
                else:self.pivots[j]=(v,mask);break
        assert len(self.pivots)==field.n
    def to_word(self,x):
        word=0
        while x:
            v,mask=self.pivots[x.bit_length()-1];x^=v;word^=mask
        return word
    def to_field(self,word):
        x=0
        while word:
            low=word&-word;x^=self.basis[low.bit_length()-1];word^=low
        return x


def make_curve(n,modulus):
    curve=Curve.__new__(Curve);curve.f=Field(n,modulus)
    a,b=2,-1
    for _ in range(2,n+1):a,b=b,-b-2*a
    curve.order=(1<<n)+1-b
    assert curve.order%4==0
    curve.cofactor=4;curve.r=curve.order//4
    return curve


def random_point(curve,rng):
    while True:
        lifts=curve.lift(rng.randrange(1,1<<curve.f.n))
        if lifts:
            p=curve.scale(lifts[rng.randrange(len(lifts))],4)
            if p:
                assert curve.valid(p) and curve.scale(p,curve.r) is None
                return p


def point_sum(curve,points):
    value=None
    for p in points:value=curve.add(value,p)
    return value


def lift_orbit(curve,normal,quotient,key):
    word=quotient.decode(key)
    x=normal.to_field(word);points=curve.lift(x)
    if not points:return []
    assert curve.scale(points[0],curve.r) is None
    orbit=set();p=points[0]
    for j in range(curve.f.n):
        orbit.add(p);orbit.add(curve.neg(p));p=curve.frob(p)
    return sorted(orbit)


def compatible_lift(curve,domains,target):
    start=time.perf_counter();pairs={};additions=0
    for p in domains[0]:
        for q in domains[1]:
            pairs.setdefault(curve.add(p,q),(p,q));additions+=1
    setup=time.perf_counter()-start
    start=time.perf_counter();probes=0;found=None
    for p in domains[2]:
        probes+=1;pair=pairs.get(curve.add(target,curve.neg(p)))
        if pair is not None:found=pair+(p,);break
    search=time.perf_counter()-start
    if found is not None:assert point_sum(curve,found)==target
    return found,{'pair_additions':additions,'pair_entries':len(pairs),'third_probes':probes,
                  'pair_setup_seconds':setup,'search_seconds':search,'verified':found is not None}


def curve_tests(n,modulus,trials):
    start=time.perf_counter();curve=make_curve(n,modulus)
    normal=NormalCoordinates(curve.f,20000+n);quotient=FrobeniusQuotient(n)
    setup=time.perf_counter()-start;rng=random.Random(30000+n)
    results=[]
    for trial in range(trials):
        source=[random_point(curve,rng) for _ in range(3)]
        keys=[quotient.encode(normal.to_word(p[0])) for p in source]
        source_phases=[rng.randrange(n) for _ in source]
        signed=[curve.frob(p,a) for p,a in zip(source,source_phases)]
        signed=[curve.neg(p) if rng.randrange(2) else p for p in signed]
        target=point_sum(curve,signed)
        start=time.perf_counter()
        domains=[lift_orbit(curve,normal,quotient,key) for key in keys]
        lifting=time.perf_counter()-start
        assert all(p in domain for p,domain in zip(signed,domains))
        found,metrics=compatible_lift(curve,domains,target)
        assert found is not None
        assert [quotient.encode(normal.to_word(p[0])) for p in found]==keys
        metrics.update(trial=trial,orbit_sizes=list(map(len,domains)),
                       quotient_and_point_lifting_seconds=lifting,target=target,
                       supplied_keys=[vars(key) for key in keys],solution=found)
        results.append(metrics)
        print('LIFT',n,trial,'pairs',metrics['pair_entries'],'seconds',round(metrics['pair_setup_seconds']+lifting,4),flush=True)
    # Demonstrate that independent point-orbit folding does not define group addition.
    for _ in range(100):
        p,q=random_point(curve,rng),random_point(curve,rng)
        sums=[curve.add(p,q),curve.add(p,curve.frob(q))]
        if all(sums):
            outs=[quotient.encode(normal.to_word(r[0])) for r in sums]
            if outs[0]!=outs[1]:
                counterexample={'p':p,'q':q,'sigma_q':curve.frob(q),
                    'sum_keys':[vars(key) for key in outs]};break
    else:raise AssertionError('expected a non-descending addition example')
    return {'degree':n,'curve_modulus':modulus,'normal_generator':normal.beta,
            'curve_order_from_trace_recurrence':curve.order,'annihilator_after_cofactor_four':curve.r,
            'primality_of_annihilator':'not required or certified in this experiment',
            'setup_seconds':setup,'trials':results,'addition_counterexample':counterexample,
            'scope':'Three valid orbit keys supplied. Tests compatible lifting, not finding quotient solutions or DLP.'}


def main(args):
    data={'parameters':vars(args),'exhaustive':[],'curve_lifting':[]}
    def save():Path(args.output).write_text(json.dumps(data,indent=2)+'\n')
    for p in (13,19):
        data['exhaustive'].append(exhaustive(p));print('EXHAUSTIVE',data['exhaustive'][-1],flush=True);save()
    data['degree131']=degree131_tests(args.random_words)
    print('DEGREE131',data['degree131'],flush=True);save()
    for n,mod in [(13,(1<<13)|(1<<4)|(1<<3)|3),
                  (19,(1<<19)|(1<<5)|(1<<2)|3),
                  (131,(1<<131)|(1<<8)|(1<<3)|(1<<2)|1)]:
        data['curve_lifting'].append(curve_tests(n,mod,args.lifting_trials));save()
    print('DONE',args.output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='cubic_quotient_results.json')
    parser.add_argument('--random-words',type=int,default=1000)
    parser.add_argument('--lifting-trials',type=int,default=3)
    main(parser.parse_args())
