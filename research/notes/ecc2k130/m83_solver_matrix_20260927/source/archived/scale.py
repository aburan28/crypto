"""Bounded scaling probe; planted witnesses only, not a relation-rate benchmark."""
import argparse, collections, hashlib, itertools, json, os, random, resource, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'heldout-phase'))
import core as c,experiment as e
from frobenius_curve_arithmetic import Field,Curve
from frobenius_cubic_quotient import CyclotomicFactorSpace
from test_cubic_quotient import NormalCoordinates,make_curve

def record(tag,**values):
    print(json.dumps({'stage':tag,**values}),flush=True)

def irreducible_polynomial(n):
    # Prime n: Rabin test x^(2^n)=x and gcd(x^2-x,f)=1.
    def sqmod(a,poly):
        out=0
        while a:
            bit=a&-a;out^=1<<(2*(bit.bit_length()-1));a^=bit
        while out.bit_length()>n:out^=poly<<(out.bit_length()-n-1)
        return out
    def gcd(a,b):
        while b:
            while a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
            a,b=b,a
        return a
    for k in range(1,n):
        poly=(1<<n)|(1<<k)|1
        if gcd(poly,sqmod(2,poly)^2)!=1:continue
        x=2
        for _ in range(n):x=sqmod(x,poly)
        if x==2:return poly
    for a,b,d in itertools.combinations(range(1,n),3):
        poly=(1<<n)|(1<<a)|(1<<b)|(1<<d)|1
        if gcd(poly,sqmod(2,poly)^2)!=1:continue
        x=2
        for _ in range(n):x=sqmod(x,poly)
        if x==2:return poly
    raise AssertionError('no sparse irreducible polynomial')

class Setup:
    def __init__(self,n,s,seed,arm,modulus):
        self.n,self.s,self.seed,self.arm=n,s,seed,arm
        self.counts=collections.Counter()
        self.curve=make_curve(n,modulus)
        self.normal=NormalCoordinates(self.curve.f,seed)
        self.space=CyclotomicFactorSpace(n,s)
        self.phases=(0,1,2)
        self.phase_maps=[[c.map_polynomial(self.space,self.normal,p,i*s,'linear_payload',self.counts)
                          for p in self.phases] for i in range(3)]
        self.nvars=3*s+n
        self.extra_polys=[];self.extra_equations=[]
        alg=c.Algebra(self.curve.f,self.counts)
        if arm=='ternary_inline':
            self.selector_offset=self.nvars;self.nvars+=6
            self.maps=[]
            for i in range(3):
                a=1<<(self.selector_offset+2*i);b=1<<(self.selector_offset+2*i+1)
                self.extra_polys.append({a|b:1});self.extra_equations.append({a|b})
                base,first,second=self.phase_maps[i]
                self.maps.append(c.addpoly(base,alg.mul({a:1},c.addpoly(first,base)),
                                              alg.mul({b:1},c.addpoly(second,base))))
        elif arm=='implicit':
            offset=self.nvars
            self.maps=[{1<<(offset+i*n+j):1<<j for j in range(n)} for i in range(3)]
            self.selector_offset=offset+3*n;self.nvars=self.selector_offset+9
            for i in range(3):
                selectors=[1<<(self.selector_offset+3*i+j) for j in range(3)]
                eqs=[{0,*selectors}]+[{selectors[j]|selectors[k]} for j in range(3) for k in range(j)]
                self.extra_equations.extend(eqs)
                self.extra_polys.extend({m:1 for m in eq} for eq in eqs)
                selected=c.addpoly(*(alg.mul({selectors[j]:1},self.phase_maps[i][j]) for j in range(3)))
                poly=c.addpoly(self.maps[i],selected)
                self.extra_polys.append(poly);self.extra_equations.extend(c.boolean_equations(poly,n))
        else:raise ValueError(arm)

    def planted(self):
        payloads=[];points=[];used=[]
        for payload in range(1,1<<self.s):
            x=self.normal.to_field(self.space.encode(payload))
            for point in self.curve.lift(x):
                if self.curve.scale(point,self.curve.r) is None:
                    used.append((payload,point))
        assert used,'no factor base points for this deterministic normal basis'
        # Fixed deterministic iteration over all choices; do not cherry-pick a fast solver witness.
        for picks in itertools.product(used,repeat=3):
            payloads=[p for p,_ in picks]
            points=[self.curve.frob(p,i) for i,(_,p) in enumerate(picks)]
            if c.proper(self.curve,points):
                target=c.group_sum(self.curve,points)
                middle=self.curve.add(points[0],points[1])
                if target is not None and middle is not None and target[0] and middle[0]:break
        else:raise AssertionError('no nondegenerate planted triple')
        return {'payloads':payloads,'witness':points,'target':target,'middle':middle,
                'usable_factor_base_lifts':len(used)}

def canonical_hash(polys,eqs,nvars):
    blob=json.dumps([[[list(x) for x in sorted(poly.items())] for poly in polys],
                    [sorted(eq) for eq in eqs],nvars],separators=(',',':')).encode()
    return hashlib.sha256(blob).hexdigest()

def child(n,s,seed,arm,budget):
    resource.setrlimit(resource.RLIMIT_AS,(1<<30,1<<30))
    t=time.perf_counter();poly=irreducible_polynomial(n)
    record('field',n=n,s=s,arm=arm,modulus=hex(poly),seconds=time.perf_counter()-t)
    t=time.perf_counter();setup=Setup(n,s,seed,arm,poly);p=setup.planted()
    record('preparation',seconds=time.perf_counter()-t,normal_generator=setup.normal.beta,
           group_order=setup.curve.order,annihilator=setup.curve.r,
           usable_factor_base_lifts=p['usable_factor_base_lifts'],target=p['target'],
           witness=p['witness'],payloads=p['payloads'],nvars=setup.nvars,
           peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    t=time.perf_counter();deadline=t+budget
    base_polys,base_eq=c.build_system(setup,p['target'],deadline)
    polys=base_polys+setup.extra_polys;eqs=base_eq+setup.extra_equations
    assignment=0
    for i,payload in enumerate(p['payloads']):assignment|=payload<<(i*s)
    assignment|=p['middle'][0]<<(3*s)
    for i in range(3):
        if arm=='implicit':
            assignment|=p['witness'][i][0]<<(3*s+n+i*n)
            assignment|=1<<(setup.selector_offset+3*i+i)
        else:
            if i:assignment|=1<<(setup.selector_offset+2*i+(i==2))
    assert all(c.evaluate(q,assignment)==0 for q in polys), 'planted witness does not satisfy equations'
    digest=canonical_hash(polys,eqs,setup.nvars)
    record('equations',seconds=time.perf_counter()-t,sha256=digest,field_polys=len(polys),
           boolean_equations=len(eqs),monomials=len(set().union(*(poly.keys() for poly in polys))),
           boolean_terms=sum(len(q) for q in eqs),max_input_degree=max(m.bit_count() for q in polys for m in q),
           nvars=setup.nvars,planted_assignment_verified=True,
           peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    try:import z3
    except ImportError:
        record("solver_unavailable",reason="z3-solver not installed in this runtime")
        return
    t=time.perf_counter();deadline=t+budget
    ctx=z3.Context();variables=[z3.Bool('b%d'%i,ctx=ctx) for i in range(setup.nvars)]
    solver=z3.Solver(ctx=ctx);solver.set(random_seed=13)
    monomials=set().union(*(poly.keys() for poly in polys))
    mons={}
    for i,m in enumerate(sorted(monomials)):
        if i%128==0:c.check(deadline)
        mons[m]=z3.And(*[variables[j] for j in range(setup.nvars) if m>>j&1]) if m else z3.BoolVal(True,ctx=ctx)
    for eq in eqs:
        c.check(deadline)
        solver.add(z3.Not(c.xor_tree(z3,[mons[m] for m in sorted(eq)],ctx=ctx)))
    record('sat_setup',seconds=time.perf_counter()-t,z3=z3.get_version_string(),
           peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    t=time.perf_counter();solver.set(timeout=5000);state=solver.check()
    record('solve',status=str(state),reason=solver.reason_unknown() if state==z3.unknown else None,
           seconds=time.perf_counter()-t,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    # SAT is only an algebraic candidate: certify it independently if reached.
    if state==z3.sat:
        model=solver.model();bits=sum(int(z3.is_true(model.eval(v,model_completion=True)))<<i for i,v in enumerate(variables))
        assert all(c.evaluate(q,bits)==0 for q in polys)
        xs=[c.evaluate(q,bits) for q in setup.maps]
        lifted=c.lift(setup.curve,xs,p['target'])
        record('candidate',verified=lifted is not None,witness=lifted)

def parent():
    root=HERE/'results_v2';root.mkdir(exist_ok=False)
    protocol={'degrees':[{'n':53,'s':4},{'n':83,'s':2}],
              'arms':['ternary_inline','implicit'],'seed':260938,
              'equation_build_and_sat_setup_budget_each_seconds':25,
              'z3_check_budget_seconds':5,
              'whole_process_timeout_seconds':70,'per_process_address_space_bytes':1<<30,
              'target_selection':'first nondegenerate planted triple in deterministic payload order',
              'interpretation':'planted satisfiable control only; no relation-rate, rank, or solving-degree claim'}
    (root/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    files=([HERE/'scale.py',HERE.parent/'heldout-phase/core.py',HERE.parent/'heldout-phase/experiment.py'] +
           list(sorted((HERE.parent/'heldout-phase/vendor').glob('*.py'))))
    (root/'source_hashes.json').write_text(json.dumps({str(q.relative_to(HERE.parent)):hashlib.sha256(q.read_bytes()).hexdigest() for q in files},indent=2)+'\n')
    for case in protocol['degrees']:
        for arm in protocol['arms']:
            n,s=case['n'],case['s'];key=f'n{n}-s{s}-{arm}'
            cmd=[sys.executable,str(HERE/'scale.py'),'--child','--n',str(n),'--s',str(s),
                 '--seed',str(protocol['seed']),'--arm',arm,'--budget',str(protocol['equation_build_and_sat_setup_budget_each_seconds'])]
            start=time.perf_counter()
            try:
                proc=subprocess.run(cmd,capture_output=True,text=True,timeout=protocol['whole_process_timeout_seconds'],
                                    env=dict(os.environ,PYTHONHASHSEED='0'))
                out={'key':key,'seconds':time.perf_counter()-start,'returncode':proc.returncode,
                     'stdout':proc.stdout,'stderr':proc.stderr}
            except subprocess.TimeoutExpired as ex:
                out={'key':key,'seconds':time.perf_counter()-start,'returncode':None,
                     'stdout':ex.stdout.decode() if isinstance(ex.stdout,bytes) else ex.stdout,
                     'stderr':ex.stderr.decode() if isinstance(ex.stderr,bytes) else ex.stderr,
                     'timeout':True}
            (root/f'{key}.json').write_text(json.dumps(out,indent=2)+'\n')
            print(key,'return',out['returncode'],'seconds',round(out['seconds'],2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--child',action='store_true')
    parser.add_argument('--n',type=int);parser.add_argument('--s',type=int)
    parser.add_argument('--arm');parser.add_argument('--seed',type=int);parser.add_argument('--budget',type=float)
    a=parser.parse_args();child(a.n,a.s,a.seed,a.arm,a.budget) if a.child else parent()
