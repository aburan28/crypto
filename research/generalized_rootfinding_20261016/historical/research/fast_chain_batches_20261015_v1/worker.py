import argparse,json,time,subprocess,sys,statistics

def main():
 t=time.monotonic();p=argparse.ArgumentParser();p.add_argument('--experiment-id',required=True);a=p.parse_args();records=[]
 try:
  for rep in range(3):
   for v in (('old','new') if rep%2==0 else ('new','old')):
    left=130-(time.monotonic()-t)
    if left<1:raise TimeoutError('wholebatchbudget')
    s=time.monotonic();r=subprocess.run([sys.executable,'-m','research.fast_chain_batches_20261015_v1.child','--variant',v],capture_output=True,text=True,timeout=left)
    if r.returncode:raise RuntimeError(r.stderr[-500:])
    records.append(dict(repetition=rep,variant=v,cold_seconds=time.monotonic()-s,result=json.loads(r.stdout)))
  grouped={}
  for r in records:
   for row in r['result']['rows']:grouped.setdefault(row['case'],[]).append((r['variant'],row))
  rows=[]
  for name,values in grouped.items():
   expected=values[0][1]['answers'];assert all(row['answers']==expected for _,row in values)
   old=[row['seconds'] for v,row in values if v=='old'];new=[row['seconds'] for v,row in values if v=='new'];rows.append(dict(case=name,old_seconds=old,new_seconds=new,median_old_over_new=statistics.median(old)/statistics.median(new),all8target_exact_root_sets_identical=True,total_root_count=sum(x['root_count'] for x in expected)))
  old=[r['cold_seconds'] for r in records if r['variant']=='old'];new=[r['cold_seconds'] for r in records if r['variant']=='new'];e=dict(cases=rows,records=records,cold_old_seconds=old,cold_new_seconds=new,cold_panel_median_old_over_new=statistics.median(old)/statistics.median(new),root_contract='complete_original8coordinate_S3chain_roots_D4_8targets',calibration_hardness_not_qualified=True,Groebner_comparison_performed=False,A_A_performed=False,hardware_isolation_qualified=False,release_gate_claim=False);status='passed';summary='Complete small-domain target batches compared old/new on seven wide prime,binary/Koblitz fields with exact output equality.'
 except Exception as x:e={'release_gate_claim':False,'partial_cold_records':records};status='inconclusive';summary=type(x).__name__+': '+str(x)
 out=dict(schema='scientific_result_v1',experiment_id=a.experiment_id,status=status,summary=summary,evidence=e,next_actions=['Retainallfieldregressions and identifysetup/unique-vs-repeated bottlenecks. Hardmatchedallroot5xgating requires genuineF4/F5/F6 onpreregisteredharddomains.']);s=json.dumps(out,sort_keys=True)
 if time.monotonic()-t>=135:out.update(status='inconclusive',summary='Latebatchserialization.',evidence={'release_gate_claim':False});s=json.dumps(out,sort_keys=True)
 print(s)
if __name__=='__main__':main()
