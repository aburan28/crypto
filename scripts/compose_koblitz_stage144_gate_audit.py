#!/usr/bin/env python3
"""Compose current gates with the rejected n59 prefetch-64 pilot."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-143-current-gate-audit-20260922';F=R/'docs/ic/runs/koblitz-n59-prefetch64-rejection-20260922.json';P=R/'docs/ic/params/k1n59-cofactor-projected-l15-prefetch64-pilot-public-59001.json';G=E/'GATE_STATUS.md';SC='koblitz_stage144_current_gate_audit.v1';SS='koblitz_stage144_current_gate_audit_seal.v1';MARK='Current through Stage 144'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage143_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage143_current_gate_audit.v1','Stage-143 replay changed');return v
def compose():
 p=replay();f=load(F,'prefetch evidence');params=load(P,'prefetch params');req(f.get('operation')=='koblitz_n59_prefetch64_rejection','identity');req(params.get('name')=='k1n59-prefetch-lookahead-pilot','params identity');req(params['collection']['units']==1 and params['collection']['max_units']==1,'pilot volume');req(f['candidate_params']['sha256']==sha(P),'params pin');patch=R/f['source_patch']['path'];req(f['source_patch']['sha256']==sha(patch),'patch pin');panel=f['panel'];req(panel['exact_relation_equality'] is True and len(panel['runs'])==8 and panel['order']==['baseline','candidate','candidate','baseline','baseline','candidate','candidate','baseline'],'panel');c=panel['comparison'];req(c['relation_unit_wall_seconds']['change_fraction']>.32 and c['relation_unit_cpu_seconds']['change_fraction']>.058 and c['whole_wall_seconds']['change_fraction']>.12,'rejection');a=f['process_accounting'];req(a['retained_processes']==8 and a['retained_total_core_seconds']>1291 and a['retained_peak_rss_bytes']==10080190464,'accounting');req(f['decision']['status']=='rejected' and f['decision']['selected_stage143_unchanged'] is True,'decision');req(MARK in G.read_text(),'marker');g=dict(p['gates']);g['1_all_stage_resource_charging']='partial_prefetch64_rejection_charged_magma_and_one_preliminary_receipt_missing';g['6_full_cost_vs_automorphism_rho']='failed_current_n59_coverage_tail_selected_prefetch64_rejected'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target prefetch-64 rejection; the exact relation stream is unchanged while scoped and inclusive costs regress, so the Stage-143 selected result is unchanged and remains far behind rho','predecessor':{'stage143_audit_sha256':sha(S/'audit.json'),'stage143_seal_sha256':sha(S/'result-seal.json'),'stage143_status':p['status']},'evidence_pins':{'prefetch64_evidence_sha256':sha(F),'candidate_params_sha256':sha(P),'candidate_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},**{k:p[k] for k in p if k.startswith('inherited_')},'inherited_current_n59_null_support_rejection':p['current_n59_null_support_rejection'],'current_n59_prefetch64_rejection':f,'gates':g,'next_targets':p['next_targets'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=compose();req(a==load(o/'audit.json','audit'),'current changed');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage144-gate-audit: {e}')
if __name__=='__main__':main()
