"""Finite complete-client timing with immutable logs, failure rows and raw telemetry."""
from pathlib import Path
import csv,json,os,re,statistics,subprocess
from experiment import ROOT,OUT,sha
from timed_client import finite
from benchmark_hardware import capture_hardware,hardware_key,check_destination

def bench(binary,label,workers=524288,steps=1024,launches=4,batch=16,monitor=True,collection=False,run_id=62173):
 binary=Path(binary).resolve();hardware=capture_hardware();key=hardware_key(hardware);target=OUT/(label+'.json');check_destination(target,hardware);assert not target.exists()
 args=[str(binary),'--packed','--threads',str(workers),'--steps',str(steps),'--launches',str(launches),'--verify','0','--run-id',str(run_id)]
 if collection:
  corpus=ROOT/'build'/('seven-node-'+label+'.dp');assert not corpus.exists();args+=['--dp-weight','34','--dp-file',str(corpus)]
 else:args+=['--bench']
 monitor_proc=None;telemetry=None
 telemetry_command=['nvidia-smi','-i',hardware['gpu_selector'],'--query-gpu=timestamp,utilization.gpu,utilization.memory,power.draw,power.limit,clocks.sm,clocks.mem,temperature.gpu','--format=csv','-lms','200']
 if monitor:
  telemetry=(OUT/(label+'-gpu.csv')).open('x');monitor_proc=subprocess.Popen(telemetry_command,stdout=telemetry,stderr=subprocess.STDOUT)
 try:process=finite(args,OUT/(label+'.log'),900,env=dict(os.environ,OMP_NUM_THREADS='8',CUDA_DISABLE_PTX_JIT='1'))
 finally:
  if monitor_proc:
   monitor_proc.terminate()
   try:monitor_proc.wait(timeout=5)
   except subprocess.TimeoutExpired:monitor_proc.kill();monitor_proc.wait()
   telemetry.close()
 text=(OUT/process['log']).read_text();expected=workers*batch*steps*launches;matches=re.findall(r'finished: ([0-9.]+) M it/s',text)
 passed=process['gpu_observation']['passed'] and process['returncode']==0 and not process['timed_out'] and 'MISMATCH' not in text and f'{expected} iterations' in text and '0 dropped' in text and len(matches)==1
 row=dict(gpu_observation=process['gpu_observation'],hardware=hardware,hardware_key=key,label=label,command=args,workers=workers,steps=steps,launches=launches,batch=batch,scalar_updates=expected,million_updates_per_second=float(matches[0]) if passed else None,wall_seconds=process['elapsed_seconds'],binary_sha256=sha(binary),passed=passed,returncode=process['returncode'],timed_out=process['timed_out'],timeout_seconds=process['timeout_seconds'],log_sha256=process['log_sha256'],process_wall_scope='All client setup, initialization, execution and process observation; client reported rate excludes setup')
 target.write_text(json.dumps(row,indent=2)+'\n');assert passed,text[-4000:]
 assert hardware_key(capture_hardware())==key
 if collection:
  data=corpus.read_bytes();assert len(data)%32==0;ordered=b''.join(sorted(data[i:i+32] for i in range(0,len(data),32)));row.update(dp_records=len(data)//32,dp_multiset_sha256=__import__('hashlib').sha256(ordered).hexdigest())
 if monitor:
  values=list(csv.reader((OUT/(label+'-gpu.csv')).open()))[1:];samples=values[15:-2];assert samples,'Run too short for steady telemetry'
  row['telemetry_sha256']=sha(OUT/(label+'-gpu.csv'));row['telemetry_command']=telemetry_command
  row['steady_telemetry']=dict(samples=len(samples),excluded_initial_samples=15,excluded_final_samples=2)
  for i,name in [(1,'gpu_busy_percent'),(2,'memory_busy_percent'),(3,'power_watts'),(5,'sm_clock_mhz'),(7,'temperature_c')]:
   v=[float(r[i].strip().split()[0]) for r in samples];row['steady_telemetry'][name]=dict(min=min(v),median=statistics.median(v),max=max(v))
 target.write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True);return row
