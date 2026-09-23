#!/usr/bin/env python3
"""Compose current gates with the rejected n59 admitted-key fast path."""
from __future__ import annotations
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-145-current-gate-audit-20260922';F=R/'docs/ic/runs/koblitz-n59-admitted-key-fastpath-rejection-20260922.json';G=E/'GATE_STATUS.md';SC='koblitz_stage146_current_gate_audit.v1';SS='koblitz_stage146_current_gate_audit_seal.v1';MARK='Current through Stage 146'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage145_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage145_current_gate_audit.v1','Stage-145 replay changed');return v
def compose():
 p=replay();f=load(F,'admitted-key evidence');req(f.get('operation')=='koblitz_n59_admitted_key_fastpath_rejection','identity');patch=R/f['source_patch']['path'];req(f['source_patch']['sha256']==sha(patch),'patch pin');one=f['one_unit_panel'];helper=f['helper_full_panel'];full=f['minimal_full_panel'];req(one['exact_relation_equality'] is True and len(one['runs'])==8,'pilot');req(helper['exact_relation_equality'] is True and len(helper['runs'])==4,'helper');req(full['exact_relation_equality'] is True and len(full['runs'])==12,'minimal');req(one['comparison']['relation_unit_wall_seconds']['change_fraction']<-.14,'pilot improvement');req(helper['comparison']['logs_wall_seconds']['change_fraction']>.63,'helper rejection');req(full['comparison']['relation_units_wall_seconds']['change_fraction']>.012 and full['comparison']['relation_units_cpu_seconds']['change_fraction']<-.004 and full['paired_relation_unit_wall_ratio_median']>1,'minimal rejection');req(all(x['verified'] is True and x['recovered_scalar']==17861472351607 for x in helper['runs']+full['runs']),'solutions');a=f['process_accounting'];req(a['retained_processes']==24 and a['retained_total_core_seconds']>9477 and a['retained_peak_rss_bytes']==10084401152,'accounting');req(f['decision']['status']=='rejected_wall_nonimprovement' and f['decision']['selected_stage145_unchanged'] is True,'decision');req(MARK in G.read_text(),'marker');g=dict(p['gates']);g['1_all_stage_resource_charging']='partial_admitted_key_rejection_charged_magma_and_one_preliminary_receipt_missing';g['6_full_cost_vs_automorphism_rho']='failed_current_n59_coverage_tail_selected_admitted_key_rejected'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target admitted-key fast-path rejection; CPU improves but six-run relation-unit wall and paired ratios do not, so the Stage-145 selected result is unchanged and remains far behind rho','predecessor':{'stage145_audit_sha256':sha(S/'audit.json'),'stage145_seal_sha256':sha(S/'result-seal.json'),'stage145_status':p['status']},'evidence_pins':{'admitted_key_evidence_sha256':sha(F),'candidate_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},**{k:p[k] for k in p if k.startswith('inherited_')},'inherited_current_n59_prefetch16_rejection':p['current_n59_prefetch16_rejection'],'current_n59_admitted_key_rejection':f,'gates':g,'next_targets':p['next_targets'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=load(o/'audit.json','audit');req(a.get('schema')==SC,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 try:r=build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage146-gate-audit: {e}')
if __name__=='__main__':main()
