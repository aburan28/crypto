#!/usr/bin/env python3
"""Compose current gates with the rejected n59 null-support tail."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-142-current-gate-audit-20260922';F=R/'docs/ic/runs/koblitz-n59-null-support-tail-rejection-20260922.json';P=R/'docs/ic/params/k1n59-cofactor-projected-l15-null-support-control-public-59001.json';G=E/'GATE_STATUS.md';SC='koblitz_stage143_current_gate_audit.v1';SS='koblitz_stage143_current_gate_audit_seal.v1';MARK='Current through Stage 143'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage142_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage142_current_gate_audit.v1','Stage-142 replay changed');return v
def compose():
 p=replay();f=load(F,'null-support evidence');params=load(P,'null-support params');req(f.get('operation')=='koblitz_n59_null_support_tail_rejection','identity');req(params.get('name')=='k1n59-null-support-prefix14-control','params identity');req(params['collection']['units']==14 and params['collection']['targeted_tail_null_support'] is True,'policy');req(f['candidate_params']['sha256']==sha(P),'params pin');patch=R/f['source_patch']['path'];req(f['source_patch']['sha256']==sha(patch),'patch pin');d=f['diagnostic'];req(d['failed_attempts_with_retained_null_support']==0 and d['final_report_has_null_support'] is False and d['relation_stream_exactly_matches_stage142_prefix14'] is True,'diagnostic');r=f['result'];req(r['verified'] is True and r['recovered_scalar']==17861472351607 and r['combined_relations']==29948,'result');req(r['ic_wall_seconds']>116,'wall boundary');a=f['process_accounting'];req(a['retained_processes']==1 and a['retained_total_core_seconds']>644 and a['retained_peak_rss_bytes']==10077847552,'accounting');req(f['decision']['status']=='rejected_no_mechanism_activation','decision');req(MARK in G.read_text(),'gate marker');g=dict(p['gates']);g['1_all_stage_resource_charging']='partial_null_support_rejection_charged_magma_and_one_preliminary_receipt_missing';g['6_full_cost_vs_automorphism_rho']='failed_current_n59_coverage_tail_selected_null_support_unavailable'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target null-support policy rejection; no failed sparse attempt yields a usable null-vector witness, so the Stage-142 selected result is unchanged and remains far behind rho','predecessor':{'stage142_audit_sha256':sha(S/'audit.json'),'stage142_seal_sha256':sha(S/'result-seal.json'),'stage142_status':p['status']},'evidence_pins':{'null_support_evidence_sha256':sha(F),'candidate_params_sha256':sha(P),'candidate_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},**{k:p[k] for k in p if k.startswith('inherited_')},'inherited_current_n59_parallel_column_rejection':p['current_n59_parallel_column_rejection'],'current_n59_null_support_rejection':f,'gates':g,'next_targets':p['next_targets'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=compose();req(a==load(o/'audit.json','audit'),'current changed');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage143-gate-audit: {e}')
if __name__=='__main__':main()
