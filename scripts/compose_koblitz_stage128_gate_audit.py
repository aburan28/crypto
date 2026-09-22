#!/usr/bin/env python3
"""Compose the current seven-gate audit with the current unified n=41
single-target full-cost IC/rho series."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

REPO=Path(__file__).resolve().parents[1]
EVIDENCE=REPO/'research/sat_factor_base_review_20260908/continuation-05-sota-gates'
STAGE127=EVIDENCE/'stage-127-current-gate-audit-20260921'
N41=REPO/'docs/ic/runs/koblitz-n41-single-target-20260921.json'
N41_PARAMS=REPO/'docs/ic/params/k0n41-fixed-algebraic-single-public-41301.json'
GATE_STATUS=EVIDENCE/'GATE_STATUS.md'
SCHEMA='koblitz_stage128_current_gate_audit.v1'
SEAL_SCHEMA='koblitz_stage128_current_gate_audit_seal.v1'
MARKER='Current through Stage 128'
class Stage128Error(RuntimeError): pass
def require(v:bool,m:str)->None:
 if not v: raise Stage128Error(m)
def load(p:Path,c:str)->dict[str,Any]:
 v=json.loads(p.read_text());require(isinstance(v,dict),f'{c} must be an object');return v
def sha256(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def replay127()->dict[str,Any]:
 r=subprocess.run([sys.executable,str(REPO/'scripts/compose_koblitz_stage127_gate_audit.py'),'verify','--output',str(STAGE127)],cwd=REPO,text=True,capture_output=True,check=True)
 v=json.loads(r.stdout);require(v.get('schema')=='koblitz_stage127_current_gate_audit.v1','Stage-127 replay changed');return v

def compose()->dict[str,Any]:
 pred=replay127();n41=load(N41,'n41 evidence');params=load(N41_PARAMS,'n41 params')
 require(n41.get('schema_version')==1,'n41 schema changed')
 require(params.get('name')=='k0n41-fixed-algebraic-single-public-41301-stage128','n41 parameter identity changed')
 require(params.get('targets')==[{'public_hash_seed':41301}],'n41 target changed')
 inst=n41['instance'];require(inst['factor_base_logs_known_by_construction'] is False,'n41 gained factor-base logs');require(inst['target_subgroup_enumerated_for_factor_base'] is False,'n41 factor base enumerated target subgroup')
 require(inst['factor_base_points']==4759 and inst['projected_columns']==29 and inst['ell']==6,'n41 base shape changed')
 require(inst['target']['target_scalar_constructed'] is False,'n41 target gained scalar')
 rel=n41['relations'];require(rel['relations']==35 and rel['trials']==114688,'n41 relation stream changed');require(rel['canonical_relations_sha256']=='8e51f7adab62f4bc351bcb35225fd6522d7cd56eeeac4f162bc00ef6c9b8058a','n41 relation hash changed')
 sol=n41['solution'];require(sol['recovered_scalar']==281099696942 and sol['verified'] is True,'n41 solution changed')
 panel=n41['panel'];require(panel['run_count']==10 and panel['all_verified'] is True,'n41 panel incomplete')
 med=panel['medians'];require(med['default']['ic_over_rho']>1 and med['onecore']['ic_over_rho']>1,'n41 full-cost loss boundary changed')
 acct=n41['science_process_accounting'];require(acct['processes']==12,'n41 accounting changed')
 require(MARKER in GATE_STATUS.read_text(),'gate-status marker changed')
 gates=dict(pred['gates']);gates['4_n31_n41_larger_pdp_scaling']='satisfied_finite_execution_coverage_current_n41_unified';gates['6_full_cost_vs_automorphism_rho']='partial_n41_and_n53_full_cost_losses_default_n53_wall_pass'
 return {
  'schema':SCHEMA,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,
  'claim_boundary':'strong internal engineering and finite public toy-research improvement; not a Koblitz index-calculus SOTA',
  'predecessor':{'stage127_audit_sha256':sha256(STAGE127/'audit.json'),'stage127_seal_sha256':sha256(STAGE127/'result-seal.json'),'stage127_status':pred['status']},
  'evidence_pins':{'n41_evidence_sha256':sha256(N41),'n41_params_sha256':sha256(N41_PARAMS),'gate_status_sha256':sha256(GATE_STATUS)},
  'inherited_phase_b_same_instance_matrix':pred['phase_b_same_instance_matrix'],
  'inherited_current_n53':pred['current_n53'],
  'current_n41':{
   'parameter_file':str(N41_PARAMS.relative_to(REPO)),'factor_base':inst,'relation_stream':rel,'solution':sol,'projected_predicate_native_cost':n41['algebraic_factor_base']['projected_predicate_cost'],
   'default_threads':med['default'],'one_thread':med['onecore'],'science_process_accounting':acct,
   'factor_base_logs_known_by_construction':False,'target_scalar_constructed_or_supplied':False,
  },
  'gates':gates,
  'next_targets':[
   'execute the frozen 160-input packet under licensed Magma F4 with one-core, total-core, wall, RSS and terminal receipts',
   'run one unified end-to-end algebraic-base IC/rho cost series in the n59 larger PDP regime',
   'reduce current n41 and n53 one-core IC/rho ratios below one without losing the n53 default-thread wall crossover',
   'obtain unaffiliated reproduction and a source-pinned novelty/correctness review',
  ],
  'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False,
 }
def write_new(p:Path,v:dict[str,Any])->None:require(not p.exists(),f'refusing to overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(out:Path)->dict[str,Any]:
 require(not out.exists(),f'refusing to overwrite {out}');out.mkdir(parents=True);a=compose();write_new(out/'audit.json',a);write_new(out/'result-seal.json',{'schema':SEAL_SCHEMA,'status':'audit_frozen','audit_sha256':sha256(out/'audit.json')});return a
def verify(out:Path)->dict[str,Any]:
 s=load(out/'result-seal.json','seal');require(s.get('schema')==SEAL_SCHEMA,'seal schema changed');require(sha256(out/'audit.json')==s.get('audit_sha256'),'audit seal changed');a=load(out/'audit.json','audit');require(a.get('schema')==SCHEMA,'audit schema changed');require(a.get('status')=='current_seven_gate_audit_verified','audit status changed');return a
def main()->None:
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='cmd',required=True);b=sub.add_parser('build');b.add_argument('--output',type=Path,required=True);v=sub.add_parser('verify');v.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Stage128Error) as e:raise SystemExit(f'stage128-gate-audit: {e}')
if __name__=='__main__':main()
