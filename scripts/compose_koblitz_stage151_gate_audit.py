#!/usr/bin/env python3
"""Compose current gates with the rejected n59 cached filter hashes."""
from __future__ import annotations
import argparse,gzip,hashlib,json,subprocess,sys,tempfile
from pathlib import Path
R=Path(__file__).resolve().parents[1];E=R/'research/sat_factor_base_review_20260908/continuation-05-sota-gates';S=E/'stage-150-current-gate-audit-20260922';F=R/'docs/ic/runs/koblitz-n59-cached-filter-hash-rejection-20260922.json';G=E/'GATE_STATUS.md';SC='koblitz_stage151_current_gate_audit.v1';SS='koblitz_stage151_current_gate_audit_seal.v1';MARK='Current through Stage 151'
SRC='src/cryptanalysis/koblitz_index_calculus.rs';PRE='cd730af687555ffff5d3b7f53b59405dce05c765';POST='f751e356d16dca674b7e495563cf85efb01046d1';SNAP=R/f'docs/ic/patches/historical_snapshots/koblitz_index_calculus.rs.{PRE}.gz';SNAP_GZ='e5d4567098f539ec470e54332fde9b22eab524d9f56a164ae44675366d3b1424';SNAP_SRC='175aa5b966bc176e281098e536f5611829fc1b01c4f51f795b8762008eefaf56';PATCH=R/'docs/ic/patches/koblitz-n59-cached-filter-hash-rejected-20260922.patch'
class Error(RuntimeError):pass
def req(v,m):
 if not v:raise Error(m)
def load(p,c):v=json.loads(p.read_text());req(isinstance(v,dict),f'{c} object');return v
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
def check_patch():
 req(sha(SNAP)==SNAP_GZ,'snapshot pin');src=gzip.decompress(SNAP.read_bytes());req(hashlib.sha256(src).hexdigest()==SNAP_SRC and blob(src)==PRE,'snapshot source');req(f'index {PRE[:8]}..{POST[:8]} ' in PATCH.read_text(),'patch preimage')
 with tempfile.TemporaryDirectory() as t:
  d=Path(t);(d/SRC).parent.mkdir(parents=True);(d/SRC).write_bytes(src);run=lambda *c:subprocess.run(['git',*c],cwd=d,text=True,capture_output=True,check=True);run('init','-q');run('apply','--check',str(PATCH));run('apply',str(PATCH));req(blob((d/SRC).read_bytes())==POST,'patch postimage')
 return {'path':str(SNAP.relative_to(R)),'sha256':SNAP_GZ,'source_path':SRC,'source_sha256':SNAP_SRC,'preimage_git_blob':PRE,'postimage_git_blob':POST,'patch_applies':True}
def replay():
 x=subprocess.run([sys.executable,str(R/'scripts/compose_koblitz_stage150_gate_audit.py'),'verify','--output',str(S)],cwd=R,text=True,capture_output=True,check=True);v=json.loads(x.stdout);req(v.get('schema')=='koblitz_stage150_current_gate_audit.v1','Stage-150 replay changed');return v
def compose():
 p=replay();f=load(F,'cached-hash evidence');req(f.get('operation')=='koblitz_n59_cached_filter_hash_rejection','identity');patch=R/f['source_patch']['path'];req(f['source_patch']['sha256']==sha(patch),'patch pin');panel=f['panel'];req(panel['exact_relation_equality'] is True and len(panel['runs'])==4,'panel');c=panel['comparison'];req(c['relation_unit_wall_seconds']['change_fraction']>.005 and c['relation_unit_cpu_seconds']['change_fraction']>.003,'rejection');a=f['process_accounting'];req(a['retained_processes']==4 and a['retained_total_core_seconds']>662 and a['retained_peak_rss_bytes']==10072932352,'accounting');req(f['decision']['status']=='rejected' and f['decision']['selected_stage147_unchanged'] is True,'decision');req(MARK in G.read_text(),'marker');g=dict(p['gates']);g['1_all_stage_resource_charging']='partial_cached_filter_hash_rejection_charged_magma_and_one_preliminary_receipt_missing';g['6_full_cost_vs_automorphism_rho']='failed_current_n59_coverage_tail_selected_cached_hash_rejected'
 return {'schema':SC,'status':'current_seven_gate_audit_verified','all_seven_gates_passed':False,'koblitz_index_calculus_sota':False,'claim_boundary':'finite same-target cached-filter-hash rejection; exact relations are preserved while block scratch traffic regresses the changed path, so the Stage-150 selected result is unchanged and remains far behind rho','predecessor':{'stage150_audit_sha256':sha(S/'audit.json'),'stage150_seal_sha256':sha(S/'result-seal.json'),'stage150_status':p['status']},'evidence_pins':{'cached_hash_evidence_sha256':sha(F),'candidate_patch_sha256':sha(patch),'gate_status_sha256':sha(G)},**{k:p[k] for k in p if k.startswith('inherited_')},'inherited_current_n59_split_filter_rejection':p['current_n59_split_filter_rejection'],'current_n59_cached_hash_rejection':f,'gates':g,'next_targets':p['next_targets'],'licensed_magma_complete':False,'independent_external_reproduction_satisfied':False,'full_cost_gate_passed':False}
def new(p,v):req(not p.exists(),f'overwrite {p}');p.write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def build(o):req(not o.exists(),f'overwrite {o}');o.mkdir(parents=True);a=compose();new(o/'audit.json',a);new(o/'result-seal.json',{'schema':SS,'status':'audit_frozen','audit_sha256':sha(o/'audit.json')});return a
def verify(o):s=load(o/'result-seal.json','seal');req(s.get('schema')==SS,'seal schema');req(sha(o/'audit.json')==s.get('audit_sha256'),'seal');a=load(o/'audit.json','audit');req(a.get('schema')==SC,'audit schema');req(a.get('status')=='current_seven_gate_audit_verified','audit status');return a
def main():
 p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
 for c in ('build','verify'):q=sp.add_parser(c);q.add_argument('--output',type=Path,required=True)
 sp.add_parser('check-patch')
 a=p.parse_args()
 try:r=check_patch() if a.cmd=='check-patch' else build(a.output.resolve()) if a.cmd=='build' else verify(a.output.resolve(strict=True));print(json.dumps(r,indent=2,sort_keys=True))
 except (OSError,ValueError,KeyError,subprocess.CalledProcessError,Error) as e:raise SystemExit(f'stage151-gate-audit: {e}')
if __name__=='__main__':main()
