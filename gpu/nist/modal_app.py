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
 .add_local_dir(LOCAL,remote_path=ROOT,copy=True,ignore=["bench","test_cpu","*.ptx","__pycache__"])
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
