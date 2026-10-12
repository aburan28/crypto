#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[2];OUT=pathlib.Path(__file__).resolve().parent;SRC=ROOT/'build/g7-cache-six-src'
builds=json.loads((OUT/'builds.json').read_text())['rows'];rows=[]
for b in builds[1:]:
 cmd=list(b['command']);i=cmd.index(str(SRC/'src/main.cu'));cmd[i]=str(SRC/'src/testsharedxcuda.cu');o=cmd.index('-o')+1;target=ROOT/'build'/f"g7-{b['label']}-raw";cmd[o]=str(target)
 t=time.monotonic();p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=600);compile_log=OUT/f"{b['label']}-raw-build.log";compile_log.write_text(p.stdout)
 row={'label':b['label'],'compile_returncode':p.returncode,'compile_seconds':time.monotonic()-t,'command':cmd,'binary':str(target.relative_to(ROOT)),'compile_log':compile_log.name}
 if p.returncode==0:
  row['binary_sha256']=hashlib.sha256(target.read_bytes()).hexdigest();t=time.monotonic();q=subprocess.run([str(target)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300);log=OUT/f"{b['label']}-raw.log";log.write_text(q.stdout);row.update(run_returncode=q.returncode,run_seconds=time.monotonic()-t,run_log=log.name,passed=q.returncode==0 and 'PASS: 136 shared-x-cache storage scenarios' in q.stdout and 'PASS: original/candidate equality' in q.stdout)
 rows.append(row);print(json.dumps(row),flush=True)
 if not row.get('passed'):break
(OUT/'units.json').write_text(json.dumps({'passed':len(rows)==2 and all(r['passed'] for r in rows),'rows':rows},indent=2)+'\n')
raise SystemExit(0 if len(rows)==2 and all(r['passed'] for r in rows) else 1)
