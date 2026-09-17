from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
rows=[]
for variant in ['indexed','immediate']:
 for tool in ['memcheck','initcheck','synccheck']:
  label=f'{variant}-{tool}'
  cmd=['/usr/local/cuda-13.2/bin/compute-sanitizer','--tool',tool,str(ROOT/'build'/f'g7-two-jump-{variant}'),'--packed','--threads','1024','--steps','2','--launches','1','--verify','0','--bench']
  p=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300)
  log=OUT/f'{label}.log';log.write_text(p.stdout)
  rows.append({'variant':variant,'tool':tool,'returncode':p.returncode,'error_summary_zero':'ERROR SUMMARY: 0 errors' in p.stdout,'log':log.name,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()})
result={'passed':all(r['returncode']==0 and r['error_summary_zero'] for r in rows),'rows':rows}
(OUT/'sanitizers.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
