#!/usr/bin/env python3
"""Compose current gates with the rejected n59 witnessed-filter width."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-135-current-gate-audit-20260921';F=R/'docs/ic/runs/koblitz-n59-witness-filter-width-20260921.json';G=E/'GATE_STATUS.md';SC='koblitz_stage136_current_gate_audit.v1';SS='koblitz_stage136_current_gate_audit_seal.v1';MARK='Current through Stage 136'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):
 v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage135_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage135_current_gate_audit.v1','Stage-135 replay changed');return v
def compose():
 p=replay();f=load(F,'filter evidence');req(f.get('operation')=='koblitz_n59_witness_filter_width_rejection','identity');req(f['instance']['public_target_seed']==59001,'target');patch=R/f['source_patch']['path'];req(sha(patch)==f['source_patch']['sha256'],'patch pin');c=f['candidate_filter'];req(c['nominal_bits_per_pair']==8 and c['allocated_filter_bytes']==1073741824,'candidate shape');req(c['relation_hash_preserved'] is True,'relation hash');req(all(x>1 for x in c['relation_unit_wall_ratios_vs_selected']),'wall boundary');req(c['decision']=='rejected_wall_regression','decision');a=f['process_accounting'];req(a['processes']==2 and a['total_core_seconds']>328 and a['peak_rss_bytes']==10076028928,'accounting');req(MARK in G.read_text(),'gate marker')
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target witnessed-filter rejection; selected Stage135 result unchanged and not a SOTA','predecessor':{'stage135_audit_sha256':sha(S/'audit.json'),'stage135_seal_sha256':sha(S/'result-seal.json'),'stage135_status':p['status']},'evidence_pins':{'filter_evidence_sha256':sha(F),'rejected_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},'inherited_phase_b_same_instance_matrix':p['inherited_phase_b_same_instance_matrix'],'inherited_current_n53':p['inherited_current_n53'],'inherited_current_n41':p['inherited_current_n41'],'inherited_current_n59_standard_cap':p['inherited_current_n59_standard_cap'],'inherited_current_n59_standard_frontier':p['inherited_current_n59_standard_frontier'],'inherited_current_n59_cofactor_projected_ell14':p['inherited_current_n59_cofactor_projected_ell14'],'inherited_current_n59_cofactor_projected_ell15_compact':p['inherited_current_n59_cofactor_projected_ell15_compact'],'inherited_current_n59_collector_tuning':p['inherited_current_n59_collector_tuning'],'inherited_current_n59_witnessed_two_pass':p['inherited_current_n59_witnessed_two_pass'],'inherited_current_n59_cached_witnessed':p['current_n59_cached_witnessed'],'current_n59_filter_rejection':f,'gates':p['gates'],'next_targets':['execute the frozen 160-input packet under licensed Magma F4 with complete resource and terminal receipts','reduce the same n59 cached witnessed full cost without hiding its 10 GB construction peak','obtain unaffiliated reproduction and a source-pinned novelty/correctness review'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=load(o/'audit.json','audit');req(a.get('schema')==SC,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage136-gate-audit: {e}')
if __name__=='__main__':main()
