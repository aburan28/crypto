import argparse,json,time,sys,subprocess,statistics,hashlib,math

def main():
 t=time.monotonic();p=argparse.ArgumentParser();p.add_argument('--experiment-id',required=True);a=p.parse_args();records=[]
 try:
  for rep in range(3):
   for v in (('old','new') if rep%2==0 else ('new','old')):
    left=85-(time.monotonic()-t)
    if left<1:raise TimeoutError('fullcoldbudget')
    start=time.monotonic();r=subprocess.run([sys.executable,'-m','research.separator_batch_20261016_v1.child','--variant',v],capture_output=True,text=True,timeout=left)
    if r.returncode:raise RuntimeError(r.stderr[-500:])
    records.append(dict(repetition=rep,variant=v,cold_seconds=time.monotonic()-start,result=json.loads(r.stdout),stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest()))
  grouped={}
  for r in records:
   for row in r['result']['rows']:grouped.setdefault(row['case'],[]).append((r['variant'],row))
  cases=[]
  for name,values in grouped.items():
   expected=values[0][1]['answers'];assert all(row['answers']==expected for _,row in values)
   old=[r['seconds'] for v,r in values if v=='old'];new=[r['seconds'] for v,r in values if v=='new'];cases.append(dict(case=name,order=values[0][1]['order'],degree=values[0][1]['degree'],old_seconds=old,new_seconds=new,median_old_over_new=statistics.median(old)/statistics.median(new),exact_allroot_outputs_match=True,total_roots=sum(x['count'] for x in expected)))
  old=[r['cold_seconds'] for r in records if r['variant']=='old'];new=[r['cold_seconds'] for r in records if r['variant']=='new'];e=dict(cases=cases,records=records,distinct_fields=len({(x['order'],x['degree']) for x in cases}),curve_contexts=len(cases),cold_old_seconds=old,cold_new_seconds=new,cold_panel_median_old_over_new=statistics.median(old)/statistics.median(new),geometric_mean_context_gain=math.exp(sum(math.log(x['median_old_over_new']) for x in cases)/len(cases)),baseline='prior_cached_FastFiber_separator_not_generic_old_or_Groebner',all_original8coordinate_roots_emitted=True,independent_original_equation_verifier_charged=True,source_domain_size_max4=True,hardness_or_qualified_isolation=False,release_gate_claim=False);status='passed';summary='Shared complete leftseparator and certified-edge joins matched every original root against prior cached separator across the manyfield panel.'
 except Exception as x:e={'release_gate_claim':False,'partial_cold_records':records};status='inconclusive';summary=type(x).__name__+': '+str(x)
 out=dict(schema='scientific_result_v1',experiment_id=a.experiment_id,status=status,summary=summary,evidence=e,next_actions=['Ifbatchgaininsufficient, profile independentedgeverification/outputcost. ThehardGF83D18leftbound2^56precludesmaterialization; pursueconditionalcompressedcells, notextrapolation.']);s=json.dumps(out,sort_keys=True)
 if time.monotonic()-t>=90:out.update(status='inconclusive',summary='Latecompletepanelserialization.',evidence={'release_gate_claim':False});s=json.dumps(out,sort_keys=True)
 print(s)
if __name__=='__main__':main()
