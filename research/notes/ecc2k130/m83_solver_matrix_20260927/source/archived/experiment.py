"""Bounded synthetic phase-compatibility experiment; no challenge inputs."""
import collections,itertools,json,math,resource,sys,time
from pathlib import Path
import core as c

CHOICES=(0,1,2)
ORIGINAL_BUILD=c.build_system

class Setup(c.Setup):
    def __init__(self,n,s,seed,arm):
        super().__init__(n,s,seed,'linear_payload')
        self.arm=arm
        self.phase_maps=[[c.map_polynomial(self.space,self.normal,p,i*s,'linear_payload',self.counts)
                          for p in CHOICES] for i in range(3)]
        union=sorted(set(p for d in self.domains for p in d))
        self.domains=[union,union,union]
        self.nvars=3*s+n
        self.extra_polys=[]; self.extra_equations=[]
        if arm=='enumerated':return
        if arm in ('ternary_lifted','ternary_inline'):
            self.selector_offset=self.nvars+(0 if arm=='ternary_inline' else 3*n)
            self.nvars=self.selector_offset+6
            alg=c.Algebra(self.curve.f,self.counts)
            if arm=='ternary_lifted':
                self.maps=[{1<<(3*s+n+i*n+j):1<<j for j in range(n)} for i in range(3)]
            else:self.maps=[]
            for i in range(3):
                a=1<<(self.selector_offset+2*i)
                b=1<<(self.selector_offset+2*i+1)
                self.extra_polys.append({a|b:1})
                self.extra_equations.append({a|b})
                base,first,second=self.phase_maps[i]
                selected=c.addpoly(base,alg.mul({a:1},c.addpoly(first,base)),
                                   alg.mul({b:1},c.addpoly(second,base)))
                if arm=='ternary_inline':self.maps.append(selected)
                else:
                    poly=c.addpoly(self.maps[i],selected)
                    self.extra_polys.append(poly)
                    self.extra_equations.extend(c.boolean_equations(poly,n))
            return
        # x and t are explicit curve-field coordinates. No substitution of
        # bilinear phase maps into S3: keep linking equations instead.
        offset=self.nvars
        self.maps=[{1<<(offset+i*n+j):1<<j for j in range(n)} for i in range(3)]
        offset+=3*n
        self.selector_offset=offset
        self.nvars=offset+9
        alg=c.Algebra(self.curve.f,self.counts)
        for i in range(3):
            selectors=[1<<(offset+3*i+j) for j in range(3)]
            eqs=[{0,*selectors}]+[{selectors[j]|selectors[k]} for j in range(3) for k in range(j)]
            self.extra_equations.extend(eqs)
            self.extra_polys.extend({m:1 for m in eq} for eq in eqs)
            selected=c.addpoly(*(alg.mul({selectors[j]:1},self.phase_maps[i][j]) for j in range(3)))
            poly=c.addpoly(self.maps[i],selected)
            self.extra_polys.append(poly)
            self.extra_equations.extend(c.boolean_equations(poly,n))
        self.key_offset=None
        if arm=='implicit_keys':
            self.key_offset=self.nvars
            basis=[self.space.q.project(word)[1] for word in self.space.basis]
            ka=c.Algebra(self.space.q.field,self.counts)
            for i in range(3):
                z={1<<(i*s+j):basis[j] for j in range(s)}
                k={1<<(self.key_offset+i*s+j):basis[j] for j in range(s)}
                eq=c.addpoly(k,ka.power(z,n))
                self.extra_polys.append(eq)
                self.extra_equations.extend(c.boolean_equations(eq,n-1))
            self.nvars+=3*s

    def certificate(self):
        out=super().certificate()
        out.update(phase_choices=list(CHOICES),variables=self.nvars,
                   selector_offset=getattr(self,'selector_offset',None),key_offset=getattr(self,'key_offset',None))
        return out


def build(setup,target,deadline=math.inf):
    polys,eqs=ORIGINAL_BUILD(setup,target,deadline)
    return polys+setup.extra_polys,eqs+setup.extra_equations

c.build_system=build


def solve(setup,record,seconds=2,context_mode="shared",solver_seed=260926):
    start=time.perf_counter(); results=[]
    phases=list(itertools.product(CHOICES,repeat=3)) if setup.arm=='enumerated' else [None]
    for phase in phases:
        left=seconds-(time.perf_counter()-start)
        if left<=0:break
        if phase is not None:
            setup.maps=[setup.phase_maps[i][p] for i,p in enumerate(phase)]
        # One shared per-target budget includes equation building and SAT.
        r=c.solve_query(setup,record,encoding_seconds=left,query_seconds=0,context_mode=context_mode,solver_seed=solver_seed)
        # Core below derives search budget from the same absolute deadline.
        results.append(r)
        if r['verified']:
            if phase is not None:r['selected_phases']=list(phase)
            elif setup.arm.startswith('ternary'):
                r['selected_phases']=[(1 if r['assignment']>>(setup.selector_offset+2*i)&1
                                      else 2 if r['assignment']>>(setup.selector_offset+2*i+1)&1
                                      else 0) for i in range(3)]
            else:r['selected_phases']=[
                next(p for p in CHOICES if r['assignment']>>(setup.selector_offset+3*i+p)&1) for i in range(3)]
            break
        if r['status'] not in ('no_decomposition',):break
    if not results:
        r=dict(record,status='timeout',verified=False,rank_gain=0)
    else:r=dict(results[-1])
    if not r['verified'] and (len(results)!=len(phases) or any(x['status']!='no_decomposition' for x in results)):
        r['status']='timeout'
    r.update(total_query_seconds=time.perf_counter()-start,phase_systems_attempted=len(results),
             candidate_models=sum(x.get('candidate_models',0) for x in results),
             peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    if results:
        r['anf_degree']=max(x.get('anf_degree',0) for x in results)
        r['anf_boolean_terms_sum']=sum(x.get('anf_boolean_terms',0) for x in results)
    return r


def main():
    import argparse,z3
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--s',type=int,required=True)
    p.add_argument('--seed',type=int,required=True);p.add_argument('--arm',choices=['enumerated','implicit','ternary_lifted','ternary_inline'],required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(1<<30,1<<30))
    t=time.perf_counter();setup=Setup(a.n,a.s,a.seed,a.arm);inputs=setup.inputs(4,2)
    head=dict(type='setup',arm=a.arm,certificate=setup.certificate(),setup_seconds=time.perf_counter()-t,
              python=sys.version,z3=z3.get_version_string())
    with a.out.open('x') as f:
        f.write(json.dumps(head)+'\n');f.flush()
        for record in inputs:
            r=solve(setup,record);f.write(json.dumps(r)+'\n');f.flush()

if __name__=='__main__':main()
