"""Small, exact WDSat reproduction. No DLP or best-solver performance claim."""
from pathlib import Path
import argparse
from datetime import datetime, timezone
from build_pilot import verified_sources
import hashlib, itertools, json, math, platform, re, resource, statistics, subprocess, time

ROOT = Path(__file__).resolve().parent
if not __debug__:
    raise RuntimeError('Run without -O: certificate checks must remain enabled')

def anf_check(text, bits):
    lines = [x for x in text.splitlines() if x.strip()]
    nv, ne = map(int, lines[0].split()[2:])
    assert len(bits) == nv and set(bits) <= {'0', '1'} and len(lines)-1 == ne
    for line in lines[1:]:
        tok=line.split(); assert tok.pop(0)=='x'
        parity=0; i=0
        while tok[i]!='0':
            t=tok[i]; i+=1
            if t=='T': parity ^= 1
            elif t.startswith('.'):
                d=int(t[1:]); term=1
                for v in tok[i:i+d]: term &= int(bits[int(v)-1])
                parity ^= term; i+=d
            else: parity ^= int(bits[int(t)-1])
        assert parity == 1, 'ANF assignment failed'
    return ne

class GF:
    def __init__(self,n,mod): self.n,self.mod=n,mod
    def mul(self,a,b):
        z=0
        while b:
            if b&1:z^=a
            b>>=1;a<<=1
            if a>>self.n:a^=self.mod
        return z
    def sq(self,a):return self.mul(a,a)
    def inv(self,a):
        if not a: raise ZeroDivisionError
        z=1;b=a;k=(1<<self.n)-2
        while k:
            if k&1:z=self.mul(z,b)
            b=self.sq(b);k>>=1
        assert self.mul(a,z)==1
        return z
    def lift(self,x):
        if x==0:return [(0,1)]
        c=x^1^self.sq(self.inv(x));z=c;v=c
        for _ in range((self.n-1)//2):v=self.sq(self.sq(v));z^=v
        if self.sq(z)^z != c:return []
        y=self.mul(x,z)
        return [(x,y),(x,y^x)]
    def on_curve(self,P):
        if P is None:return True
        x,y=P
        return self.sq(y)^self.mul(x,y)==self.mul(self.sq(x),x)^self.sq(x)^1
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u:
            if y!=v or x==0:return None
            lam=x^self.mul(y,self.inv(x));w=self.sq(lam)^lam^1
            ans=(w,self.sq(x)^self.mul(lam^1,w))
        else:
            lam=self.mul(y^v,self.inv(x^u));w=self.sq(lam)^lam^x^u^1
            ans=(w,self.mul(lam,x^w)^w^y)
        assert self.on_curve(ans)
        return ans
    def s4(self,x,y,z,r):
        mul=self.mul;sq=self.sq
        e1=x^y^z;e2=mul(x,y)^mul(x,z)^mul(y,z);e3=mul(mul(x,y),z)
        a=sq(e1);b=sq(e2);c=sq(e3);r2=sq(r);r3=mul(r2,r);r4=sq(r2)
        return (r4^sq(a)^sq(c)^mul(sq(b),r4)^mul(mul(c,e3),r)
                ^mul(mul(e3,b),r3)^mul(mul(e3,a),r)^mul(e3,r3)
                ^mul(mul(a,c),r2)^mul(c,r4)^c^mul(b,r2))

def little(s):return sum(int(x)<<i for i,x in enumerate(s))

def check_sat(field, raw, bits, l, r):
    anf_check(raw, bits)
    xs = [little(bits[i*l:(i+1)*l]) for i in range(3)]
    assert field.s4(*xs, r) == 0, 'S4 witness failed'
    for pts in itertools.product(*(field.lift(x) for x in xs)):
        assert all(field.on_curve(P) for P in pts)
        Q = field.add(field.add(pts[0], pts[1]), pts[2])
        if Q is not None and Q[0] == r:
            return xs, {'summands': pts, 'sum': Q}
    raise AssertionError('SAT assignment has no verified point relation')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', required=True, type=Path)
    parser.add_argument('--solver', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path,
                        help='New directory; existing evidence cannot be overwritten.')
    args = parser.parse_args()
    if args.output.exists(): parser.error('--output must not already exist')
    cache = args.cache_dir.resolve(); solver = args.solver.resolve()
    verified_sources(cache)
    build = json.loads((solver.parent/'build.json').read_text())
    binary_hash = hashlib.sha256(solver.read_bytes()).hexdigest()
    assert build['exit_code'] == 0 and build['binary_sha256'] == binary_hash
    assert build['source_manifest_sha256'] == hashlib.sha256((ROOT/'source_manifest.json').read_bytes()).hexdigest()
    dest=args.output.resolve();dest.mkdir(parents=True, exist_ok=False)
    compiler=build['compiler']
    cpu=next((l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),'unknown')
    output={'status':'running','date':datetime.now(timezone.utc).isoformat(),'cpu':cpu,'platform':platform.platform(),
            'compiler':compiler,'flags':'-O3 -Wall; no OpenMP','repeats':3,
            'build':build,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'solver_commit':'61c6ff3f49445af2729274819edb27b85a89efc1',
            'corpus_commit':'a053971cb155fbab30eaed3e2d11df850b4f3fbb',
            'scope':'Selected polynomial instances; not uniform random targets, full relation collection, or a complete ECDLP.',
            'cases':[]}
    for case in ['n15l5-1-S','n15l5-11-U','n17l6-1-S','n17l6-11-U','n19l6-1-S','n19l6-11-U']:
        src=cache/'instances'/('X'+case+'.anf');raw=src.read_text()
        normalized='\n'.join(x for x in raw.splitlines() if x.strip())+'\n'
        inp=dest/'normalized'/src.name;inp.parent.mkdir(exist_ok=True);inp.write_text(normalized)
        info=(cache/'instances'/('INFO'+case+'.dimacs')).read_text().splitlines()
        n,l=map(int,info[0].split());mod=little(info[1]);r=little(info[2]);field=GF(n,mod)
        cmd=[str(solver),'-i',str(inp),'-n',str(n),'-l',str(l),'-m','3','-b']
        runs=[]
        for rep in range(3):
            before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.perf_counter()
            try:
                p=subprocess.run(cmd,capture_output=True,text=True,timeout=60)
            except subprocess.TimeoutExpired as exc:
                for suffix, data in [('stdout', exc.stdout), ('stderr', exc.stderr)]:
                    (dest/(case+f'-{rep}.'+suffix)).write_bytes(data or b'')
                output.update(status='TIMEOUT', failed_case=case, failed_repeat=rep)
                (dest/'pilot_results.json').write_text(json.dumps(output,indent=2)+'\n')
                raise SystemExit('pilot timed out; incomplete run retained')
            wall=time.perf_counter()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
            (dest/(case+f'-{rep}.stdout')).write_text(p.stdout)
            (dest/(case+f'-{rep}.stderr')).write_text(p.stderr)
            assert p.returncode==0,(case,p.stderr)
            nv=int(normalized.splitlines()[0].split()[2])
            assignment=next((x for x in p.stdout.splitlines() if len(x)==nv and set(x)<=set('01')),None)
            status='SAT' if assignment is not None else 'UNSAT' if re.search(r'^UNSAT$',p.stdout,re.M) else 'ERROR'
            assert status!='ERROR',p.stdout
            conflicts=int(p.stdout.strip().splitlines()[-1])
            if assignment:anf_check(normalized,assignment)
            runs.append({'wall_s':wall,'child_cpu_s':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                         'status':status,'conflicts':conflicts,'assignment':assignment})
        assert len({x['status'] for x in runs})==1
        check_start=time.perf_counter();poly_tests=0;witness=None;point_witness=None
        if runs[0]['status']=='SAT':
            for run in runs:
                witness, point_witness = check_sat(field, normalized, run['assignment'], l, r)
        else:
            for xs in itertools.combinations_with_replacement(range(1<<l),3):
                poly_tests+=1
                if field.s4(*xs,r)==0:raise AssertionError(('false UNSAT',case,xs))
        row={'case':case,'n':n,'l':l,'m':3,'anf_variables':int(normalized.splitlines()[0].split()[2]),
             'status':runs[0]['status'],'median_wall_s':statistics.median(x['wall_s'] for x in runs),
             'median_child_cpu_s':statistics.median(x['child_cpu_s'] for x in runs),
             'conflicts':runs[0]['conflicts'],'runs':runs,'all_anf_equations_verified':witness is not None,
             'point_witness':point_witness,'exhaustive_s4_checks':poly_tests,'verification_s':time.perf_counter()-check_start,
             'normalized_sha256':hashlib.sha256(normalized.encode()).hexdigest(),
             'config_sha256':hashlib.sha256((cache/'vendor/WDSat/src/config.h').read_bytes()).hexdigest()}
        output['cases'].append(row)
        (dest/'pilot_results.json').write_text(json.dumps(output,indent=2)+'\n')
        print(json.dumps({k:row[k] for k in ['case','status','median_wall_s','median_child_cpu_s','conflicts','exhaustive_s4_checks']})+f' point_verified={point_witness is not None}',flush=True)
    output['status']='measured_local_small_pilot'
    (dest/'pilot_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print('complete',flush=True)

if __name__=='__main__':main()
