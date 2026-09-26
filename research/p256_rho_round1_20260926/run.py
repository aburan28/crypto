#!/usr/bin/env python3
"""Bounded, reproducible P-256 rho arithmetic experiments (round 1).

This harness deliberately separates:
  * full-width P-256 field microbenchmarks (E1-E5);
  * small prime-order short-Weierstrass toy groups for complete rho solves
    (E6-E11/E15), so collision quality can be measured without pretending a
    full P-256 DLP is feasible;
  * architecture gates (E12-E14), which report only what this host can measure.

No result is promoted from wall time alone.  JSONL contains every seed.
"""
from __future__ import annotations
import argparse, json, math, os, platform, random, statistics, time
from dataclasses import dataclass
from typing import Optional

P = (1<<256) - (1<<224) + (1<<192) + (1<<96) - 1
MASK=(1<<256)-1
SEED=0x50323536

def mont_constants(p=P):
    R=1<<256
    return R, pow(R,-1,p), (R*R)%p
R,RINV,R2=mont_constants()

def mul_python(a,b): return (a*b)%P
def mul_mont_roundtrip(a,b):
    # Correctness/control backend, not a claim about CIOS instruction count.
    am=(a*R)%P; bm=(b*R)%P
    cm=(am*bm*RINV)%P
    return (cm*RINV)%P

def solinas_fold(x):
    """Reference generalized-Mersenne fold using 2^256 == 2^224-2^192-2^96+1."""
    while x.bit_length()>256:
        lo=x&MASK; hi=x>>256
        x=lo + hi + (hi<<224) - (hi<<192) - (hi<<96)
        if x<0:
            # Keep the next fold non-negative without changing the residue.
            x += ((-x)//P+1)*P
    return x%P

def mul_solinas(a,b): return solinas_fold(a*b)

def weak_add(a,b,k=2):
    # A bounded redundant representative; reduction is deliberately by kp.
    m=k*P
    x=a+b
    return x-m if x>=m else x

def field_correctness(samples=10000):
    rng=random.Random(SEED)
    edge=[0,1,2,P-2,P-1,(1<<255),(1<<224)-1]
    pairs=[(a,b) for a in edge for b in edge]
    pairs += [(rng.randrange(P),rng.randrange(P)) for _ in range(samples)]
    for i,(a,b) in enumerate(pairs):
        ref=mul_python(a,b)
        assert mul_solinas(a,b)==ref, ("solinas",i)
        assert mul_mont_roundtrip(a,b)==ref, ("mont",i)
    return len(pairs)

def bench_field(samples=30000, rounds=3):
    rng=random.Random(SEED)
    pairs=[(rng.randrange(P),rng.randrange(P)) for _ in range(samples)]
    out=[]
    for name,fn in [("python_mod",mul_python),("mont_roundtrip_control",mul_mont_roundtrip),("solinas_fold",mul_solinas)]:
        ts=[]; chk=0
        for _ in range(rounds):
            t=time.perf_counter_ns()
            z=0
            for a,b in pairs: z ^= fn(a,b)&0xffffffff
            ts.append(time.perf_counter_ns()-t); chk ^= z
        out.append({"experiment":"E1_E2_field","backend":name,"samples":samples,
                    "median_ns_per_mul":statistics.median(ts)/samples,"checksum":chk,
                    "stage_diagnostic":True})
    return out

@dataclass(frozen=True)
class Curve:
    p:int; a:int; b:int; n:int; G:tuple[int,int]; name:str

def inv(x,p): return pow(x,p-2,p)
def add(E,A,B):
    if A is None:return B
    if B is None:return A
    x1,y1=A;x2,y2=B;p=E.p
    if x1==x2 and (y1+y2)%p==0:return None
    if A==B:
        if y1==0:return None
        m=(3*x1*x1+E.a)*inv(2*y1%p,p)%p
    else:m=(y2-y1)*inv((x2-x1)%p,p)%p
    x3=(m*m-x1-x2)%p
    return x3,(m*(x1-x3)-y1)%p
def mul(E,k,A):
    Q=None
    while k:
        if k&1:Q=add(E,Q,A)
        A=add(E,A,A);k>>=1
    return Q

def isprime(n):
    if n<2:return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%q==0:return n==q
    d=n-1;s=0
    while d%2==0:s+=1;d//=2
    for a in (2,3,5,7,11):
        if a>=n:continue
        x=pow(a,d,n)
        if x in (1,n-1):continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1:break
        else:return False
    return True

def toy_curve(bits=19):
    # Deterministically search a small prime field and curve with prime order.
    # Counting is intentionally bounded to toy sizes.
    rng=random.Random(SEED+bits)
    p=(1<<bits)-1
    while not isprime(p):p-=2
    for trial in range(200):
        a=rng.randrange(p);b=rng.randrange(1,p)
        if (4*a*a*a+27*b*b)%p==0:continue
        pts=[]
        count=1
        for x in range(p):
            rhs=(x*x*x+a*x+b)%p
            ls=pow(rhs,(p-1)//2,p) if rhs else 0
            count += 1 if rhs==0 else (2 if ls==1 else 0)
        if not isprime(count):continue
        for x in range(p):
            rhs=(x*x*x+a*x+b)%p
            if rhs and pow(rhs,(p-1)//2,p)==1:
                # p chosen 3 mod 4 when possible; fallback brute sqrt at toy scale.
                if p%4==3:y=pow(rhs,(p+1)//4,p)
                else:
                    y=next((y for y in range(p) if y*y%p==rhs),None)
                    if y is None:continue
                E=Curve(p,a,b,count,(x,y),f"toy{bits}")
                if mul(E,count,E.G) is None:return E
    raise RuntimeError("no prime-order toy curve found")

def table(E,Q,Rn,seed):
    rng=random.Random(seed); out=[]
    for _ in range(Rn):
        c=rng.randrange(1,E.n);d=rng.randrange(1,E.n)
        out.append((add(E,mul(E,c,E.G),mul(E,d,Q)),c,d))
    return out

def part(Pt,Rn,kind):
    x=Pt[0]
    if kind=="xlow":h=x
    elif kind=="xor2":h=(x^(x>>16))
    else:h=(x+0x9e3779b9*(x>>8))
    return h&(Rn-1)

def rho_solve(E,secret,Rn=32,kind="xlow",dp_bits=5,seed=0,max_steps=2_000_000):
    Q=mul(E,secret,E.G); T=table(E,Q,Rn,seed)
    rng=random.Random(seed^0xA5A5)
    seen={}; steps=0; dps=0
    # Multiple independent walks; affine is canonical here.
    for restart in range(10000):
        aa=rng.randrange(E.n);bb=rng.randrange(E.n)
        X=add(E,mul(E,aa,E.G),mul(E,bb,Q))
        if X is None:continue
        for _ in range(1<<(dp_bits+3)):
            j=part(X,Rn,kind); M,c,d=T[j]
            X=add(E,X,M);aa=(aa+c)%E.n;bb=(bb+d)%E.n;steps+=1
            if X is None:break
            # negation folding: canonical y and coefficient sign.
            if X[1]>E.p//2:
                X=(X[0],(-X[1])%E.p);aa=(-aa)%E.n;bb=(-bb)%E.n
            if (X[0]&((1<<dp_bits)-1))==0:
                dps+=1; key=X
                if key in seen:
                    a0,b0=seen[key]
                    den=(bb-b0)%E.n
                    if den:
                        k=(a0-aa)*pow(den,-1,E.n)%E.n
                        if mul(E,k,E.G)==Q:
                            return {"verified":True,"steps":steps,"dps":dps,"secret":secret,"recovered":k}
                else:seen[key]=(aa,bb)
                break
            if steps>=max_steps:return {"verified":False,"steps":steps,"dps":dps}
    return {"verified":False,"steps":steps,"dps":dps}

def partition_quality(E,Rn,kind,samples=20000):
    rng=random.Random(SEED); hist=[0]*Rn
    X=E.G
    same=0;prev=None
    for _ in range(samples):
        X=mul(E,rng.randrange(1,E.n),E.G)
        j=part(X,Rn,kind);hist[j]+=1
        same += (prev==j);prev=j
    exp=samples/Rn
    chi=sum((v-exp)**2/exp for v in hist)
    return {"experiment":"E5_E7_partition","partition":kind,"R":Rn,"samples":samples,
            "chi2":chi,"same_transition_rate":same/max(1,samples-1),
            "min_bucket":min(hist),"max_bucket":max(hist)}

def weak_reduction_panel(samples=100000):
    rng=random.Random(SEED); out=[]
    vals=[rng.randrange(P) for _ in range(samples)]
    for k in (2,4,8):
        x=0;mx=0;t=time.perf_counter_ns()
        for v in vals:
            x=weak_add(x,v,k);mx=max(mx,x)
        dt=time.perf_counter_ns()-t
        out.append({"experiment":"E3_E4_weak","k":k,"samples":samples,
                    "ns_per_add":dt/samples,"max_over_p":mx/P,
                    "within_bound":0<=x<k*P,"stage_diagnostic":True})
    return out

def run(out_path,bits=17,seeds=16):
    rows=[]
    checked=field_correctness()
    rows.append({"experiment":"correctness","field_pairs":checked,"verified":True})
    rows += bench_field()
    rows += weak_reduction_panel()
    E=toy_curve(bits)
    rows.append({"experiment":"toy_curve","name":E.name,"p":E.p,"n":E.n,"a":E.a,"b":E.b,"G":E.G})
    for Rn in (8,16,32,64):
        for kind in ("xlow","xor2","mix"):
            rows.append(partition_quality(E,Rn,kind,5000))
    for Rn in (8,16,32,64):
        for kind in ("xlow","xor2","mix"):
            vals=[]
            t0=time.perf_counter_ns()
            for s in range(seeds):
                secret=2+((SEED+s*0x9e37)%(E.n-3))
                z=rho_solve(E,secret,Rn,kind,5,SEED+s)
                z.update({"experiment":"E6_E8_E15_rho","seed":s,"R":Rn,"partition":kind,
                          "n":E.n,"fold":2})
                vals.append(z);rows.append(z)
            good=[x for x in vals if x["verified"]]
            elapsed=time.perf_counter_ns()-t0
            rows.append({"experiment":"E15_summary","R":Rn,"partition":kind,
                         "verified":len(good),"runs":seeds,
                         "mean_steps":statistics.mean(x["steps"] for x in good) if good else None,
                         "median_steps":statistics.median(x["steps"] for x in good) if good else None,
                         "wall_ns":elapsed,
                         "classification":"engineering diagnostic"})
    rows.append({"experiment":"E12_host","platform":platform.platform(),"python":platform.python_version(),
                 "cpu_count":os.cpu_count(),"measured":True})
    rows.append({"experiment":"E13_gpu","measured":False,"reason":"no GPU backend in this Python round"})
    rows.append({"experiment":"E14_fpga","measured":False,"reason":"requires synthesis/board toolchain"})
    with open(out_path,"w") as f:
        for r in rows:f.write(json.dumps(r,sort_keys=True)+"\n")
    return rows

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",default="results.jsonl")
    ap.add_argument("--toy-bits",type=int,default=17)
    ap.add_argument("--seeds",type=int,default=16)
    a=ap.parse_args()
    rows=run(a.out,a.toy_bits,a.seeds)
    print(json.dumps({"rows":len(rows),"out":a.out,
                      "verified_rho":sum(r.get("verified") is True for r in rows if r.get("experiment")=="E6_E8_E15_rho")},indent=2))
