#!/usr/bin/env python3
"""Compose current gates with the rejected n59 parallel-column tail."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-141-current-gate-audit-20260922';F=R/'docs/ic/runs/koblitz-n59-parallel-column-tail-rejection-20260922.json';P=R/'docs/ic/params/k1n59-cofactor-projected-l15-parallel-column-control-public-59001.json';G=E/'GATE_STATUS.md';SC='koblitz_stage142_current_gate_audit.v1';SS='koblitz_stage142_current_gate_audit_seal.v1';MARK='Current through Stage 142'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):
 v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage141_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage141_current_gate_audit.v1','Stage-141 replay changed');return v
def compose():
 p=replay();f=load(F,'parallel-column evidence');params=load(P,'parallel-column params');req(f.get('operation')=='koblitz_n59_parallel_column_tail_rejection','identity');req(params.get('name')=='k1n59-parallel-columns-prefix14-control','params identity');req(params['collection']['units']==14 and params['collection']['targeted_tail_rank_columns']==64 and params['collection']['targeted_tail_parallel_columns'] is True,'candidate policy');req(f['candidate_params']['sha256']==sha(P),'params pin');patch=R/f['source_patch']['path'];req(f['source_patch']['sha256']==sha(patch),'patch pin');req(f['relation_stream']['exact_match_between_candidates'] is True and f['relation_stream']['combined_relations']==29948,'relation stream');runs=f['candidate_runs'];req(len(runs)==2 and all(x['verified'] is True and x['recovered_scalar']==17861472351607 for x in runs),'solutions');req(runs[0]['ic_wall_seconds']>118 and runs[1]['ic_wall_seconds']>138,'wall rejection');a=f['process_accounting'];req(a['retained_processes']==2 and a['retained_total_core_seconds']>1309 and a['retained_peak_rss_bytes']==10082500608,'accounting');req(f['decision']['status']=='rejected' and f['decision']['selected_stage141_unchanged'] is True,'decision');req(MARK in G.read_text(),'gate marker');gates=dict(p['gates']);gates['1_all_stage_resource_charging']='partial_parallel_column_rejection_charged_magma_and_one_preliminary_receipt_missing';gates['6_full_cost_vs_automorphism_rho']='failed_current_n59_coverage_tail_selected_parallel_columns_rejected'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target parallel-column scheduling rejection; both nested and single-layer candidates regress full cost, so the Stage-141 selected result is unchanged and remains far behind rho','predecessor':{'stage141_audit_sha256':sha(S/'audit.json'),'stage141_seal_sha256':sha(S/'result-seal.json'),'stage141_status':p['status']},'evidence_pins':{'parallel_column_evidence_sha256':sha(F),'candidate_params_sha256':sha(P),'candidate_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},**{k:p[k] for k in p if k.startswith('inherited_')},'inherited_current_n59_rank_solve_cadence_rejection':p['current_n59_rank_solve_cadence_rejection'],'current_n59_parallel_column_rejection':f,'gates':gates,'next_targets':p['next_targets'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=load(o/'audit.json','audit');req(a.get('schema')==SC,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage142-gate-audit: {e}')
if __name__=='__main__':main()
