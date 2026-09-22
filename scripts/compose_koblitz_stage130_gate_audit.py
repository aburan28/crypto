#!/usr/bin/env python3
"""Compose current Koblitz gates with the exact n59 standard-subspace
cofactor-class and memory frontier."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
from typing import Any
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-129-current-gate-audit-20260921';F=R/'docs/ic/runs/koblitz-n59-standard-class-frontier-20260921.json';G=E/'GATE_STATUS.md';SC='koblitz_stage130_current_gate_audit.v1';SS='koblitz_stage130_current_gate_audit_seal.v1';MARK='Current through Stage 130'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):
 v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage129_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage129_current_gate_audit.v1','stage129 replay');return v
def compose():
 p=replay();f=load(F,'frontier');req(f.get('operation')=='koblitz_n59_standard_subspace_class_aware_frontier','frontier identity');rows=f['class_counts'];req([x['ell'] for x in rows]==list(range(9,17)),'width coverage');fr=f['frontier'];req(fr['widest_fitting_ell']==15 and fr['first_nonfitting_ell']==16,'memory frontier');req(fr['best_fitting_probe_lower_bound']>1.5e9,'probe lower bound');req(fr['widest_fitting_compact_bytes']<4<<30 and fr['first_nonfitting_compact_bytes']>4<<30,'budget crossing');req(f['process_accounting']['processes']==16,'accounting');req(MARK in G.read_text(),'gate marker')
 gates=dict(p['gates']);gates['4_n31_n41_larger_pdp_scaling']='satisfied_finite_n59_attempt_and_standard_family_frontier';gates['6_full_cost_vs_automorphism_rho']='partial_n59_standard_materialized_frontier_no_go_other_families_open'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite no-go frontier for materialized standard n59 bases under 4 GiB; not a generic or asymptotic lower bound','predecessor':{'stage129_audit_sha256':sha(S/'audit.json'),'stage129_seal_sha256':sha(S/'result-seal.json'),'stage129_status':p['status']},'evidence_pins':{'frontier_sha256':sha(F),'gate_status_sha256':sha(G)},'inherited_phase_b_same_instance_matrix':p['inherited_phase_b_same_instance_matrix'],'inherited_current_n53':p['inherited_current_n53'],'inherited_current_n41':p['inherited_current_n41'],'inherited_current_n59':p['current_n59'],'current_n59_standard_frontier':f,'gates':gates,'next_targets':['execute the frozen 160-input packet under licensed Magma F4 with complete resources and terminals','search a different algebraic or implicit n59 factor-base family that avoids the cofactor-yield/materialized-table frontier','reduce current n41 and n53 one-core IC/rho ratios below one without losing the n53 default-thread wall crossover','obtain unaffiliated reproduction and a source-pinned novelty/correctness review'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=load(o/'audit.json','audit');req(a.get('schema')==SC,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True);b=sp.add_parser('build');b.add_argument('--output',type=Path,required=True);v=sp.add_parser('verify');v.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage130-gate-audit: {e}')
if __name__=='__main__':main()
