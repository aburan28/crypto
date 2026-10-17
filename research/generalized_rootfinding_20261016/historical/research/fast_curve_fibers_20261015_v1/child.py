import argparse,json,time,random
from .panel import cases
from .core import FastFiber
from research.advance.general_chain import GeneralS3Fiber
p=argparse.ArgumentParser();p.add_argument('--variant',choices=['old','new'],required=True);p.add_argument('--repeats',type=int,required=True);a=p.parse_args();started=time.monotonic();rows=[]
for name,f,curve in cases():
 rng=random.Random(2026101500+len(rows));pairs=[(0,0),(0,1),(1,1)]+[(rng.randrange(f.order),rng.randrange(f.order)) for _ in range(5)]
 t=time.monotonic();fiber=(GeneralS3Fiber if a.variant=='old' else FastFiber)(f,curve);roots=[]
 for rep in range(a.repeats):
  for x,y in pairs:roots.append(list(fiber.roots(x,y)))
 rows.append(dict(case=name,characteristic=str(f.characteristic),degree=f.degree,order=str(f.order),seconds=time.monotonic()-t,canonical_roots=roots,cache_hits=getattr(getattr(fiber,'context',None),'hits',0),cache_misses=getattr(getattr(fiber,'context',None),'misses',0)))
print(json.dumps(dict(rows=rows,seconds_including_field_construction=time.monotonic()-started,variant=a.variant,repeats=a.repeats),sort_keys=True))
