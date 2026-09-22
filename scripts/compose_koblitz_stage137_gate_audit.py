#!/usr/bin/env python3
"""Compose current gates with the selected n59 uncovered-column tail."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

R=Path(__file__).resolve().parents[1]
E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates'
S=E/'stage-136-current-gate-audit-20260921'
T=R/'docs/ic/runs/koblitz-n59-targeted-tail-20260922.json'
P=R/'docs/ic/params/k1n59-cofactor-projected-l15-targeted-tail-public-59001.json'
G=E/'GATE_STATUS.md'
SC='koblitz_stage137_current_gate_audit.v1';SS='koblitz_stage137_current_gate_audit_seal.v1';MARK='Current through Stage 137'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):
 v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage136_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage136_current_gate_audit.v1','Stage-136 replay changed');return v
def compose():
 p=replay();t=load(T,'targeted evidence');params=load(P,'targeted params');req(t.get('operation')=='koblitz_n59_uncovered_column_targeted_tail','evidence identity');req(params.get('name')=='k1n59-cofactor-projected-l15-targeted-tail-selected-public-59001-stage137','params identity');req(params.get('targets')==[{'public_hash_seed':59001}],'target');req(t['params']['sha256']==sha(P),'params pin')
 fb=t['factor_base'];req(fb['points']==32934 and fb['projected_columns']==16344,'factor base shape');req(fb['factor_base_logs_known_by_construction'] is False and fb['target_subgroup_enumerated'] is False and fb['scalar_preimages_retained'] is False,'knowledge boundary')
 rows=t['uniform_prefix_frontier']['rows'];req([x['units'] for x in rows]==[22,23,24,25],'prefix widths');req([x['uncovered_columns'] for x in rows]==[5,4,1,1],'prefix frontier')
 policy=t['targeted_policy'];req(policy['selected_missing_column']==15726 and policy['selected_factor_point_index']==21212,'targeted predicate');req(policy['selected_relation']=={'trial':4967,'a':6472497976388,'points':[3990,21212,22115]},'targeted relation');req(policy['selected_relation_sha256']=='c83ba65ae1727ae5aac4b8506a7ffe2d96e81fb06805b3bc32ea971fb2ffaadc','targeted hash')
 stream=t['relation_stream'];req(stream['uniform_units']==25 and stream['uniform_trials']==2500000 and stream['targeted_pair_lookups']==10000 and stream['combined_relations']==52635,'stream shape')
 for name in ('selected_default_thread','selected_one_worker'):
  x=t[name];req(x['verified'] is True and x['recovered_scalar']==17861472351607,f'{name} solution');req(x['uniform_relations_sha256']==stream['uniform_relations_sha256'] and x['targeted_relations_sha256']==stream['targeted_relations_sha256'],f'{name} hashes');req(x['ic_over_rho_wall_ratio']>1,f'{name} rho boundary')
 md=t['matched_ab']['default_thread']['comparison'];mo=t['matched_ab']['one_worker']['comparison'];req(md['ic_wall_reduction_fraction']>.28 and md['whole_cpu_reduction_fraction']>.02,'default matched improvement');req(mo['ic_wall_reduction_fraction']>.009 and mo['whole_cpu_reduction_fraction']>.006,'one-worker matched improvement')
 a=t['process_accounting'];req(a['retained_processes']==23 and a['retained_total_core_seconds']>6492 and a['retained_peak_rss_bytes']==10088611840,'process accounting');req(a['unretained_preliminary']['processes']==1 and a['unretained_preliminary']['core_seconds'] is None,'missing receipt boundary');req(MARK in G.read_text(),'gate marker')
 gates=dict(p['gates']);gates['1_all_stage_resource_charging']='partial_missing_licensed_magma_and_one_preliminary_receipt';gates['3_single_core_core_memory_conflicts_wall']='partial_current_n59_targeted_complete_magma_missing';gates['6_full_cost_vs_automorphism_rho']='failed_current_n59_targeted_matched_improvement_far_behind_rho'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite single-target uncovered-column tail with matched default and one-worker improvements; full cost remains far behind rho, one preliminary receipt is missing, and this is not a SOTA','predecessor':{'stage136_audit_sha256':sha(S/'audit.json'),'stage136_seal_sha256':sha(S/'result-seal.json'),'stage136_status':p['status']},'evidence_pins':{'targeted_evidence_sha256':sha(T),'params_sha256':sha(P),'gate_status_sha256':sha(G)},'inherited_phase_b_same_instance_matrix':p['inherited_phase_b_same_instance_matrix'],'inherited_current_n53':p['inherited_current_n53'],'inherited_current_n41':p['inherited_current_n41'],'inherited_current_n59_standard_cap':p['inherited_current_n59_standard_cap'],'inherited_current_n59_standard_frontier':p['inherited_current_n59_standard_frontier'],'inherited_current_n59_cofactor_projected_ell14':p['inherited_current_n59_cofactor_projected_ell14'],'inherited_current_n59_cofactor_projected_ell15_compact':p['inherited_current_n59_cofactor_projected_ell15_compact'],'inherited_current_n59_collector_tuning':p['inherited_current_n59_collector_tuning'],'inherited_current_n59_witnessed_two_pass':p['inherited_current_n59_witnessed_two_pass'],'inherited_current_n59_cached_witnessed':p['inherited_current_n59_cached_witnessed'],'inherited_current_n59_filter_rejection':p['current_n59_filter_rejection'],'current_n59_targeted_tail':t,'gates':gates,'next_targets':['execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts','reduce the same n59 full cost while charging factor-base discovery, the 10 GB table and targeted-tail policy search','obtain unaffiliated reproduction and a source-pinned novelty/correctness review'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=compose();req(a==load(o/'audit.json','audit'),'current changed');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage137-gate-audit: {e}')
if __name__=='__main__':main()
