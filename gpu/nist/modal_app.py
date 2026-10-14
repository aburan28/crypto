"""RTX PRO 6000 Blackwell runner for gpu/nist.

modal run modal_app.py::validate
modal run modal_app.py::tune --repeats 3
"""
import os,pathlib,re,statistics,subprocess,modal
CUDA=os.environ.get("NIST_CUDA_VERSION","13.3.1")
GPU=os.environ.get("NIST_GPU","RTX-PRO-6000")
ROOT="/root/nist";LOCAL=pathlib.Path(__file__).parent
image=(modal.Image.from_registry(f"nvidia/cuda:{CUDA}-devel-ubuntu24.04",add_python="3.12")
 .entrypoint([]).apt_install("build-essential","libboost-dev")
  .add_local_dir(LOCAL,remote_path=ROOT,copy=True,ignore=["bench","test_cpu","test_rho","*.ptx","__pycache__"])
 .run_commands(f"cd {ROOT} && make test"))
app=modal.App("nist-prime-gpu-arithmetic")
def sh(cmd,timeout=3600):
 r=subprocess.run(cmd,shell=True,cwd=ROOT,text=True,capture_output=True,timeout=timeout)
 return r.returncode,r.stdout+r.stderr
def build(ptx,maxreg=0):
 return sh(f"make -B bench ARCH=sm_120 NIST_PTX={int(ptx)} MAXREG={maxreg}")
@app.function(image=image,gpu=GPU,timeout=2*3600)
def validate():
 rc,out=sh("make test");print(out)
 if rc: raise RuntimeError("CPU oracle failed")
 for ptx in (0,1):
  rc,out=build(ptx,0);print(out)
  if rc: raise RuntimeError(f"device build failed ptx={ptx}")
  rc,out=sh("./bench 128");print(out)
  if rc: raise RuntimeError(f"device run failed ptx={ptx}")
 return {"ok":True}
RATE=re.compile(r"^(P-(?:256|384) (?:Montgomery|Solinas)) (mul chain1|mul chain2|sqr|point_double|mixed_add) threads=(\d+): ([0-9.]+) ([GM])",re.M)
REG=re.compile(r"Used (\d+) registers")
@app.function(image=image,gpu=GPU,timeout=2*3600)
def sass():
  rc,out=sh("make -B bench ARCH=sm_120 NIST_PTX=1 NVCCFLAGS='-O3 -std=c++17 -arch=sm_120 -Xptxas=-v -ptx'");print(out[-2000:] if rc==0 else out)
  if rc: raise RuntimeError("ptx build failed")
  ptx=open("bench.ptx").read()
  import re
  m=re.search(r"\.visible \.entry _Z5k_rhoI7EcCurveI9NistFieldI7P256Mod.*?\{(.*?)\}\s*$",ptx,re.S)
  print("k_rho P256 PTX body (first 3000 chars):")
  print(m.group(1)[:3000] if m else "NOT FOUND")
  if rc: raise RuntimeError("build failed")
  rc,out=sh("cuobjdump -sass bench | awk '/k_rhoI7EcCurve/{p=1} p{print} /Function : .*/ && !/k_rhoI/ && p && !first {exit} {if(/Function :/) first=1}' | head -80");print(out)
  rc,out=sh("cuobjdump -sass bench | grep -c k_rhoI7EcCurve");print("instantiations:",out)
  rc,out=sh("cuobjdump -sass bench");open("/tmp/nist_sass.txt","w").write(out)
  import re
  sass=open("/tmp/nist_sass.txt").read()
  for m in re.finditer(r"Function : (_Z5k_rho\S+)",sass):
    name=m.group(1); start=m.end(); nxt=sass.find("\n\t\tFunction :",start)
    body=sass[start:nxt if nxt>0 else len(sass)]
    print(f"=== {name}: {len(body.splitlines())} lines ===")
    print(body[:2000])
  return {"ok":True}
@app.function(image=image,gpu=GPU,timeout=2*3600)
def run_ptx(ptx:int=1,threads:int=128):
  rc,log=build(ptx,0)
  rc2,out=sh(f"./bench {threads}")
  print(out)
  return {"rc":rc,"log_tail":log[-1500:],"out":out}
@app.function(image=image,gpu=GPU,timeout=4*3600)
def sweep_rho(reps:int=1):
  """threads/block x walks-per-thread sweep of the rho walk, NIST_PTX=1.
  Short iteration counts keep each cell fast while DPs stay statistically
  checkable against the 2^-24 dp rate."""
  rc,log=build(1,0);print(log[-600:])
  if rc: raise RuntimeError("build failed")
  import re
  regs=re.findall(r"Compiling entry function '(\S*k_rhoI\S*)'.*?Used (\d+) registers.*?(\d+) bytes cumulative stack",log,re.S)
  print("k_rho instantiations (name, regs, stack):",regs)
  rows=[]
  for th in (64,128,256):
    for W in (2,4,8,16):
      iters=8192
      for rep in range(reps):
        rc,out=sh(f"./bench {th} rho {W} {iters}")
        if rc: print("FAILED",th,W,rep,out[-500:]); continue
        for curve,r,thr,blk,ww,gstep,dps in re.findall(
          r"^(P-(?:256|384)) rho: walks=(\d+) threads/blk=(\d+) blocks=(\d+) W=(\d+)\s+([0-9.]+) Gstep/s\s+\((\d+) DPs",out,re.M):
          rows.append({"threads":th,"W":int(W),"repeat":rep,"curve":curve,"Gstep_s":float(gstep),
                      "walks":int(r),"blocks":int(blk),"DPs":int(dps),"expected_DPs":int(r)*iters/2**24})
        print(out)
  print("\n=== best per curve ===")
  best={}
  for r_ in rows:
    k=r_["curve"]
    if k not in best or r_["Gstep_s"]>best[k]["Gstep_s"]: best[k]=r_
  for k,v in best.items(): print(k,v)
  return {"rows":rows}
@app.function(image=image,gpu=GPU,timeout=4*3600)
def sweep_rho2():
  """Round 2: register-resident batch inversion, W up to 32, maxrregcount.
  Cells: (threads, W, maxreg)."""
  import re
  results=[]
  for maxreg in (0,96,128,160):
    rc,log=build(1,maxreg)
    if rc: raise RuntimeError(f"build maxreg={maxreg} failed")
    regs=re.findall(r"Compiling entry function '(\S*k_rhoI\S*)'.*?Used (\d+) registers.*?(\d+) bytes cumulative stack",log,re.S)
    print("maxreg",maxreg,"k_rho:",regs)
    cells=[(128,8),(128,16),(128,32),(256,16),(256,32)] if maxreg==0 else [(256,32),(256,16),(128,32)]
    for th,W in cells:
      rc,out=sh(f"./bench {th} rho {W} 8192")
      if rc: print("FAILED",maxreg,th,W,out[-400:]); continue
      print(out.strip())
      for curve,r,thr,blk,ww,gstep,dps in re.findall(
        r"^(P-(?:256|384)) rho: walks=(\d+) threads/blk=(\d+) blocks=(\d+) W=(\d+)\s+([0-9.]+) Gstep/s\s+\((\d+) DPs",out,re.M):
        results.append({"maxreg":maxreg,"threads":th,"W":W,"curve":curve,"Gstep_s":float(gstep),"DPs":int(dps),
                        "expected_DPs":int(r)*8192/2**24})
  print("\n=== best per curve ===")
  best={}
  for r_ in results:
    k=r_["curve"]
    if k not in best or r_["Gstep_s"]>best[k]["Gstep_s"]: best[k]=r_
  for k,v in best.items(): print(k,v)
  return {"rows":results}
@app.function(image=image,gpu=GPU,timeout=4*3600)
def tune(repeats:int=3):
 rows=[]
 rc,info=sh("nvidia-smi --query-gpu=name,compute_cap,clocks.max.sm,memory.total,memory.used,power.limit --format=csv,noheader")
 print(info)
 configs=[(0,0),(1,0),(1,96),(1,128),(1,160)]
 for ptx,maxreg in configs:
  rc,log=build(ptx,maxreg);print(log)
  if rc: continue
  regs=[int(x) for x in REG.findall(log)]
  for th in (64,96,128,160,192,256):
   for rep in range(repeats):
    rc,out=sh(f"./bench {th}");print(out)
    if rc: continue
    for curve,op,t,v,u in RATE.findall(out):
     rate=float(v)*(1000.0 if u=="G" else 1.0)
     rows.append({"ptx":ptx,"maxreg":maxreg,"threads":int(t),"repeat":rep,
      "curve":curve,"op":op,"Mop_s":rate,"ptxas_registers":regs})
 print("\n=== median-ranked ===")
 for curve in ("P-256 Montgomery","P-256 Solinas","P-384 Montgomery","P-384 Solinas"):
  for op in ("mul chain1","mul chain2","sqr","point_double","mixed_add"):
   groups={}
   for r in rows:
    if r["curve"]==curve and r["op"]==op:
     k=(r["ptx"],r["maxreg"],r["threads"]);groups.setdefault(k,[]).append(r["Mop_s"])
   ranked=sorted(((statistics.median(v),k,v) for k,v in groups.items()),reverse=True)
   if ranked:print(curve,op,{"median_Mop_s":ranked[0][0],"config":ranked[0][1],"samples":ranked[0][2]})
 return {"gpu":info.strip(),"rows":rows}
