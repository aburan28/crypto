#!/usr/bin/env python3
import hashlib,json,pathlib,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[2];OUT=pathlib.Path(__file__).resolve().parent;rows=[]
for slots in (5,6):
 binary=ROOT/'build'/f'g7-xcache-{slots}-raw'
 for tool,extra in [('memcheck',[]),('initcheck',['--initcheck-address-space','shared'])]:
  label=f'xcache-{slots}-{tool}'+('-shared' if extra else '')
  cmd=['/usr/local/cuda-13.2/bin/compute-sanitizer','--tool',tool,'--error-exitcode','99',*extra,str(binary)]
  t=time.monotonic();p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=600);log=OUT/f'{label}.log';log.write_text(p.stdout)
  passed=p.returncode==0 and 'ERROR SUMMARY: 0 errors' in p.stdout and 'PASS: 136 shared-x-cache storage scenarios' in p.stdout
  row={'label':label,'command':cmd,'returncode':p.returncode,'elapsed_seconds':time.monotonic()-t,'passed':passed,'log':log.name,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()};rows.append(row);print(json.dumps(row),flush=True)
  if not passed:break
 if not rows[-1]['passed']:break
(OUT/'sanitizers.json').write_text(json.dumps({'passed':len(rows)==4 and all(r['passed'] for r in rows),'rows':rows},indent=2)+'\n')
raise SystemExit(0 if len(rows)==4 and all(r['passed'] for r in rows) else 1)
