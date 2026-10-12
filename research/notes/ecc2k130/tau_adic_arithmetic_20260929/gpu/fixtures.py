"""Independent integer-polynomial oracle and immutable PR #913 fixtures."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
spec = importlib.util.spec_from_file_location("tau_reference", PARENT / "benchmark.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)  # Recoding helpers do not require Sage.
sys.path.insert(0, str(HERE.parents[4] / "tools"))
from curve_identity import curve_identity


def sha(data):
    return hashlib.sha256(data).hexdigest()


class Field:
    def __init__(self,m):
        self.m=m
        self.poly=(1<<m)|sum(1<<i for i in ref.TAPS[m])
    def mul(self,a,b):
        r=0
        while b:
            if b&1:r^=a
            b>>=1;a<<=1
            if a>>self.m:a^=self.poly
        return r
    def inv(self,a):
        if a==0:raise ZeroDivisionError()
        u,v,g,h=a,self.poly,1,0
        while u!=1:
            j=u.bit_length()-v.bit_length()
            if j<0:u,v=v,u;g,h=h,g;j=-j
            u^=v<<j;g^=h<<j
        return g
    def add(self,p,q):
        if p is None:return q
        if q is None:return p
        x,y=p;u,v=q
        if x==u:
            if y!=v or not x:return None
            lam=x^self.mul(y,self.inv(x))
            xx=self.mul(lam,lam)^lam
            return xx,self.mul(x,x)^self.mul(lam^1,xx)
        lam=self.mul(y^v,self.inv(x^u))
        xx=self.mul(lam,lam)^lam^x^u
        return xx,self.mul(lam,x^xx)^xx^y
    def multiply(self,p,k):
        r=None
        if k<0:p=(p[0],p[0]^p[1]) if p else None;k=-k
        for bit in bin(k)[2:]:
            r=self.add(r,r)
            if bit=='1':r=self.add(r,p)
        return r


def encoded(p):return None if p is None else [hex(p[0]),hex(p[1])]


def panels():
    raw=(PARENT/'results/run-02.json').read_bytes()
    assert sha(raw)=='f108dfa628f141b7da5c2885472aa960a25403a8f157e84a74cac2e58be8b797', 'Archived receipt changed'
    receipt=json.loads(raw)
    selected=[p for p in receipt['panels'] if p['m'] in (83,131)]
    for p in selected:
        f=Field(p['m']);start=time.perf_counter()
        cases=[(tuple(int(v,16) for v in c['point']),int(c['scalar'])) for c in p['cases']]
        expected=[f.multiply(point,k) for point,k in cases]
        digest=sha(json.dumps([encoded(q) for q in expected],sort_keys=True).encode())
        assert digest==p['output_sha256'], 'Independent oracle disagrees with archived Sage'
        generator=next(q for q in selected if q['m']==p['m'] and not q['holdout'])['cases'][0]['point']
        field={'characteristic':2,'degree':p['m'],'representation':'polynomial',
               'modulus_exponents':[p['m'],*ref.TAPS[p['m']]],'element_encoding':'hex polynomial coefficient bitset, least significant bit is constant'}
        curve={'model':'binary Weierstrass','coefficients':[1,0,0,0,1],
               'subgroup_order':str(int(p['group_order'])//4),'cofactor':4,
               'generator':generator,'target_group':'prime-order subgroup'}
        identity=curve_identity(field,curve,'e0')
        yield {'m':p['m'],'holdout':p['holdout'],'seed':p['seed'],'cases':cases,'expected':expected,
               'sage_output_sha256':digest,'source_receipt_sha256':sha(raw),
               'input_sha256':p['input_sha256'],'identity':identity,'field':field,'curve':curve,
               'oracle_seconds':time.perf_counter()-start}


def pack(points,m):
    w=(m+63)//64
    out=np.zeros((2*w+1,len(points)),dtype=np.uint64)
    for i,p in enumerate(points):
        if p is None:out[2*w,i]=1;continue
        for j in range(w):
            out[2*j,i]=(p[0]>>(64*j))&((1<<64)-1)
            out[2*j+1,i]=(p[1]>>(64*j))&((1<<64)-1)
    return out


def columns(m):
    f=Field(m);w=(m+63)//64
    return np.array([[(f.mul(1<<i,1<<i)>>(64*j))&((1<<64)-1) for j in range(w)] for i in range(m)],dtype=np.uint64)


def prepare(panel,method,repeat=1):
    start=time.perf_counter();m=panel['m'];cases=panel['cases']
    if method=='binary_naf':digits=[ref.binary_digits(k,True) for _,k in cases]
    elif method=='reduced_tau_naf':
        delta=ref.annihilator(m,-1)
        digits=[ref.tau_digits(ref.reduced(k,delta,-1),-1) for _,k in cases]
    else:raise ValueError(method)
    recode=time.perf_counter()-start
    idx=np.tile(np.arange(len(cases)),repeat)
    np.random.default_rng(panel['seed']).shuffle(idx)
    n=len(idx);ds=np.zeros((max(map(len,digits),default=0),n),dtype=np.int8)
    lengths=np.array([len(digits[i]) for i in idx],dtype=np.int32)
    for j,i in enumerate(idx):ds[:len(digits[i]),j]=digits[i]
    return (pack([cases[i][0] for i in idx],m),ds,lengths,
            pack([panel['expected'][i] for i in idx],m),
            {'unique_recode_seconds':recode,'recoding_reused':repeat>1,
             'prepare_seconds':time.perf_counter()-start,'cases':n,'unique_cases':len(cases),
             'repetitions_of_fixtures':repeat,'packed_input_sha256':sha(pack([cases[i][0] for i in idx],m).tobytes()+ds.tobytes()+lengths.tobytes())})
