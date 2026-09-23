#!/usr/bin/env python3
"""Compose current gates with the rejected selected FIFO admission fast path."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-154-current-gate-audit-20260922';F=R/'docs/ic/runs/koblitz-n59-fifo-admitted-fastpath-rejection-20260922.json';G=E/'GATE_STATUS.md';SC='koblitz_stage155_current_gate_audit.v1';SS='koblitz_stage155_current_gate_audit_seal.v1';MARK='Current through Stage 155'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage154_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage154_current_gate_audit.v1','Stage-154 replay changed');return v
def compose():
 p=replay();f=load(F,'FIFO admitted evidence');req(f.get('operation')=='koblitz_n59_fifo_admitted_fastpath_rejection','identity');patch=R/f['source_patch']['path'];req(f['source_patch']['sha256']==sha(patch),'patch pin');panel=f['panel'];req(panel['exact_relation_equality'] is True and len(panel['runs'])==4,'panel');c=panel['comparison'];req(c['relation_unit_wall_seconds']['change_fraction']>.019 and c['relation_unit_cpu_seconds']['change_fraction']>.055,'rejection');a=f['process_accounting'];req(a['retained_processes']==4 and a['retained_total_core_seconds']>686 and a['retained_peak_rss_bytes']==10085679104,'accounting');req(f['decision']['status']=='rejected_wall_and_cpu_regression' and f['decision']['selected_stage149_unchanged'] is True,'decision');req(MARK in G.read_text(),'marker');g=dict(p['gates']);g['1_all_stage_resource_charging']='partial_fifo_admitted_rejection_charged_magma_and_one_preliminary_receipt_missing';g['6_full_cost_vs_automorphism_rho']='failed_current_n59_bucket_prefetch_selected_fifo_admitted_rejected'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target selected-FIFO admitted fast-path rejection; both wall and CPU regress, so Stage 149 remains selected and far behind rho','predecessor':{'stage154_audit_sha256':sha(S/'audit.json'),'stage154_seal_sha256':sha(S/'result-seal.json'),'stage154_status':p['status']},'evidence_pins':{'fifo_admitted_evidence_sha256':sha(F),'candidate_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},**{k:p[k] for k in p if k.startswith('inherited_')},'inherited_current_n59_fifo_lookahead16_rejection':p['current_n59_fifo_lookahead16_rejection'],'current_n59_fifo_admitted_rejection':f,'gates':g,'next_targets':p['next_targets'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=load(o/'audit.json','audit');req(a.get('schema')==SC,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage155-gate-audit: {e}')
if __name__=='__main__':main()
