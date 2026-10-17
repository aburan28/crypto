import argparse,json,subprocess,sys,time,statistics,math,hashlib
from .panel import cases

def run(started,repeats):
 records=[]
 for rep in range(3):
  for variant in (('old','new') if rep%2==0 else ('new','old')):
   remaining=130-(time.monotonic()-started)
   if remaining<1:raise TimeoutError('wholepanelbudget')
   t=time.monotonic();p=subprocess.run([sys.executable,'-m','research.fast_curve_fibers_20261015_v1.child','--variant',variant,'--repeats',str(repeats)],capture_output=True,text=True,timeout=remaining)
   if p.returncode:raise RuntimeError(p.stderr[-500:])
   r=json.loads(p.stdout);records.append(dict(repetition=rep,variant=variant,cold_seconds=time.monotonic()-t,result=r,stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest()))
 groups={}
 for r in records:
  for row in r['result']['rows']:groups.setdefault(row['case'],[]).append((r['variant'],row))
 rows=[]
 for name,values in groups.items():
  expected=values[0][1]['canonical_roots'];assert all(v['canonical_roots']==expected for _,v in values)
  old=[v['seconds'] for k,v in values if k=='old'];new=[v['seconds'] for k,v in values if k=='new'];ratio=statistics.median(old)/statistics.median(new)
  rows.append(dict(case=name,degree=values[0][1]['degree'],order=values[0][1]['order'],old_seconds=old,new_seconds=new,median_old_over_new=ratio,meets_5x_component_target=ratio>=5,exact_fiber_output_equality=True))
 old=[r['cold_seconds'] for r in records if r['variant']=='old'];new=[r['cold_seconds'] for r in records if r['variant']=='new']
 return dict(field_cases=rows,distinct_fields=len({(row['order'],row['degree']) for row in rows}),curve_contexts=len(rows),calls_per_context=8*repeats,three_cold_process_pairs=True,cold_old_seconds=old,cold_new_seconds=new,cold_panel_median_old_over_new=statistics.median(old)/statistics.median(new),geometric_mean_context_gain=math.exp(sum(math.log(r['median_old_over_new']) for r in rows)/len(rows)),contexts_meeting_5x=sum(r['meets_5x_component_target'] for r in rows),records=records,release_gate_claim=False,full_system_or_Groebner_speed_claim=False,hardware_isolation_qualified=False,A_A_performed=False,calibration_input_known_before_comparison=True)
def main():
 t=time.monotonic();p=argparse.ArgumentParser();p.add_argument('--experiment-id',required=True);p.add_argument('--repeats',type=int,choices=[1,8],required=True);a=p.parse_args()
 try:e=run(t,a.repeats);status='passed';summary='Matched exact finite-fiber panel measured on prime,binary,Koblitz and odd-extension fields; all contexts retained.'
 except Exception as x:e={'release_gate_claim':False,'full_system_or_Groebner_speed_claim':False};status='inconclusive';summary=type(x).__name__+': '+str(x)
 r=dict(schema='scientific_result_v1',experiment_id=a.experiment_id,status=status,summary=summary,evidence=e,next_actions=['Use perfieldresults and worstslowdown,notjustaggregate; integrate into complete source/auxiliary chain and hardmatchedF4/F5/F6 before5xallrootsclaim.']);s=json.dumps(r,sort_keys=True)
 if time.monotonic()-t>=135:r.update(status='inconclusive',summary='Late completepanel serialization.',evidence={'release_gate_claim':False});s=json.dumps(r,sort_keys=True)
 print(s)
if __name__=='__main__':main()
