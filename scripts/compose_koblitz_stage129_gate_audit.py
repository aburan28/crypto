#!/usr/bin/env python3
"""Compose the current seven-gate audit with the bounded n=59 standard-
subspace end-to-end attempt and same-target rho panel."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
from typing import Any
REPO=Path(__file__).resolve().parents[1];E=REPO/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S128=E/'stage-128-current-gate-audit-20260921';N59=REPO/'docs/ic/runs/koblitz-n59-standard-cap-20260921.json';P59=REPO/'docs/ic/params/k1n59-standard-l9-m3-cap-public-59001.json';G=E/'GATE_STATUS.md';SCHEMA='koblitz_stage129_current_gate_audit.v1';SEAL='koblitz_stage129_current_gate_audit_seal.v1';MARK='Current through Stage 129'
class Error(RuntimeError):pass
def req(v:bool,m:str)->None:
 if not v:raise Error(m)
def load(p:Path,c:str)->dict[str,Any]:
 v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} must be object');return v
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def replay()->dict[str,Any]:
 r=subprocess.run([sys.executable,str(REPO/'scripts/compose_koblitz_stage128_gate_audit.py'),'verify','--output',str(S128)],cwd=REPO,text=True,capture_output=True,check=True);v=json.loads(r.stdout);req(v.get('schema')=='koblitz_stage128_current_gate_audit.v1','Stage-128 replay changed');return v
def compose()->dict[str,Any]:
 pred=replay();n=load(N59,'n59 evidence');p=load(P59,'n59 params');req(n.get('schema_version')==1,'n59 schema changed');req(p.get('name')=='k1n59-standard-l9-m3-eight-unit-cap-public-59001-stage129','n59 params changed');req(p.get('targets')==[{'public_hash_seed':59001}],'n59 target changed')
 fb=n['factor_base'];req(fb['kind']=='standard_subspace' and fb['ell']==9 and fb['points']==483 and fb['projected_columns']==231,'n59 factor base changed');req(fb['target_subgroup_enumerated'] is False and fb['discrete_log_labels_used'] is False,'n59 factor-base boundary changed');req(fb['folded_table_forbidden'] is True and fb['pair_table_tier']=='compact','n59 table boundary changed')
 bc=n['bounded_collection'];req(bc['trial_cap']==1200000 and bc['relations']==0,'n59 cap outcome changed');req(bc['linear_algebra']['attempts']==0,'n59 LA boundary changed');req(bc['status']=='cap_reached_zero_relations_censored','n59 censorship changed')
 rho=n['rho'];req(rho['summary']['runs']==5 and rho['summary']['all_verified'] is True,'n59 rho panel incomplete');req(rho['summary']['recovered_scalars']==[17861472351607],'n59 rho scalar changed');req(n['public_target']['target_scalar_constructed_or_supplied'] is False,'n59 target gained scalar');req(MARK in G.read_text(),'gate marker changed')
 gates=dict(pred['gates']);gates['4_n31_n41_larger_pdp_scaling']='satisfied_finite_coverage_n59_end_to_end_attempt_censored';gates['6_full_cost_vs_automorphism_rho']='partial_n59_ic_incomplete_after_cap_n41_n53_losses_default_n53_wall_pass'
 return {'schema':SCHEMA,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'strong internal engineering and finite public toy-research improvement; not a Koblitz index-calculus SOTA','predecessor':{'stage128_audit_sha256':sha(S128/'audit.json'),'stage128_seal_sha256':sha(S128/'result-seal.json'),'stage128_status':pred['status']},'evidence_pins':{'n59_evidence_sha256':sha(N59),'n59_params_sha256':sha(P59),'gate_status_sha256':sha(G)},'inherited_phase_b_same_instance_matrix':pred['inherited_phase_b_same_instance_matrix'],'inherited_current_n53':pred['inherited_current_n53'],'inherited_current_n41':pred['current_n41'],'current_n59':n,'gates':gates,'next_targets':['execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts','find an algebraic larger-PDP base with non-negligible natural relation yield or prove the materialized-base route infeasible under a stated bound','reduce current n41 and n53 one-core IC/rho ratios below one without losing the n53 default-thread wall crossover','obtain unaffiliated reproduction and a source-pinned novelty/correctness review'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p:Path,v:dict[str,Any])->None:req(not p.exists(),f'refusing overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o:Path)->dict[str,Any]:req(not o.exists(),f'refusing overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SEAL,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o:Path)->dict[str,Any]:s=load(o/'result-seal.json','seal');req(s.get('schema')==SEAL,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'audit seal');a=load(o/'audit.json','audit');req(a.get('schema')==SCHEMA,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main()->None:
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True);b=sp.add_parser('build');b.add_argument('--output',type=Path,required=True);v=sp.add_parser('verify');v.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage129-gate-audit: {e}')
if __name__=='__main__':main()
