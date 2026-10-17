"""Package existing evidence with exhaustive content-addressed reconstruction.
This script performs no experimental or statistical computation.
"""
import json,hashlib,zipfile,shutil,stat
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];STUDY=ROOT/'research/n83_basis_setup_20261015';D=STUDY/'delivery';PRIVATE=Path('/workspace/n83_basis_setup')
BUNDLE=D/'n83_basis_setup_evidence.zip';MANIFEST=D/'HASH_MANIFEST.json';RECEIPT=D/'EVIDENCE_ARCHIVE_RECEIPT.json'
excluded={BUNDLE.resolve(),MANIFEST.resolve(),RECEIPT.resolve()}
paths=set()
def add(p):
 p=Path(p).resolve()
 if p.is_dir():
  for f in p.rglob('*'):
   if f.is_file() and '__pycache__' not in f.parts and f.suffix!='.pyc' and f.resolve() not in excluded:paths.add(f.resolve())
 else:
  assert p.is_file(),f'Missing required payload {p}'
  if p not in excluded:paths.add(p)
add(STUDY);add(ROOT/'benchmarks/20261015-n83-basis-setup');add(ROOT/'benchmarks/README.md')
for p in ['AGENTS.md','Cargo.toml','Cargo.lock','src','examples/n83_basis_setup.rs','examples/n83_basis_rho.rs','scripts/build_ic_leaderboard.py','scripts/render_n83_basis_diagnostics.py','scripts/build_lab_browser.py','scripts/site/panel_index.py','docs/index-calculus-scoreboard.html','docs/ic/progress-timeline.json','docs/ic/leaderboard.json','docs/ic/LEADERBOARD.md','docs/ic-leaderboard.html','docs/ic-current-state.html','docs/browser/data.json','docs/ic/boundary_targets.json','problems/n83-eight-target-rho-contract.md']:add(ROOT/p)
# Every originally named frozen input/source/executable remains represented.
for manifest in [STUDY/'provenance/SOURCE_SHA256SUMS',STUDY/'provenance/BINARY_SHA256SUMS',STUDY/'provenance/holdout/FREEZE_SHA256SUMS',STUDY/'checks/cold/POSTPROCESSING_V1_SHA256SUMS',STUDY/'checks/final/AGGREGATE_V2_SHA256SUMS']:
 for line in manifest.read_text().splitlines():
  if not line.strip():continue
  expected,name=line.split(maxsplit=1);p=Path(name.strip());p=p if p.is_absolute() else ROOT/p
  assert hashlib.sha256(p.read_bytes()).hexdigest()==expected,(str(p),'frozen seal mismatch');add(p)
final=json.loads((STUDY/'reports/final-20261016/native_final.json').read_text())
for v in final['bindings']:
 p=Path(v['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256'];add(p)
for p in [PRIVATE/'frozen_binaries',PRIVATE/'postprocessing_binaries/v1',PRIVATE/'final_aggregate_v2',PRIVATE/'test_final_aggregate_v2',PRIVATE/'launch_cold.sh',PRIVATE/'launch_cohort.sh',Path('/workspace/n83_evaluation/frozen_binaries/isolated_bench'),Path('/workspace/report_output/unified_curve_basis_report_matched_ic_rho.pdf')]:add(p)
rows=[];objects={}
for p in sorted(paths):
 data=p.read_bytes();h=hashlib.sha256(data).hexdigest();objects.setdefault(h,p)
 relative=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else 'external/'+str(p).lstrip('/')
 rows.append({'path':relative,'original_path':str(p),'sha256':h,'size':len(data),'mode':stat.S_IMODE(p.stat().st_mode),'object':'objects/'+h})
assert len({v['path'] for v in rows})==len(rows)
manifest={'schema':'n83_basis_evidence_reconstruction/v1','coverage':'Every present file in the study, all55attempted raw directories, old/input/source/binary seals, final reports and canonical diagnostics. No repeated object omitted from original path inventory. Missing answers remain failures, not absent payload files.','planned_processes':55,'actual_completed':50,'actual_failed':5,'canonical_answers':92,'rank_log_audits':30,'AA_pairs':5,'files':rows,'original_paths':len(rows),'unique_content_objects':len(objects),'excluded_circular_metadata':['delivery/n83_basis_setup_evidence.zip','delivery/HASH_MANIFEST.json','delivery/EVIDENCE_ARCHIVE_RECEIPT.json']}
manifest_bytes=(json.dumps(manifest,indent=2)+'\n').encode();MANIFEST.write_bytes(manifest_bytes)
with zipfile.ZipFile(BUNDLE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9,allowZip64=True) as z:
 z.writestr('RECONSTRUCTION_MANIFEST.json',manifest_bytes)
 for h,p in sorted(objects.items()):z.write(p,'objects/'+h)
assert BUNDLE.stat().st_size<100*1024*1024,'Archive exceeds100MiB; retain it and adjust packaging without omitting evidence'
with zipfile.ZipFile(BUNDLE) as z:
 assert z.testzip() is None,'CRC failure'
 assert z.read('RECONSTRUCTION_MANIFEST.json')==manifest_bytes
 checked={}
 for name in z.namelist():
  if not name.startswith('objects/'):continue
  data=z.read(name);h=hashlib.sha256(data).hexdigest();assert name=='objects/'+h;checked[h]=len(data)
 for v in rows:assert checked[v['sha256']]==v['size']
receipt={'status':'PASS_CRC_ALL_OBJECT_SHA256_ALL_RECONSTRUCTION_PATHS','archive_sha256':hashlib.sha256(BUNDLE.read_bytes()).hexdigest(),'archive_bytes':BUNDLE.stat().st_size,'original_paths':len(rows),'checked_unique_objects':len(checked),'every_path_size_and_hash_verified':True,'reconstruction_manifest_sha256':hashlib.sha256(manifest_bytes).hexdigest(),'no_scientific_values_recomputed':True}
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
for src,name in [(BUNDLE,BUNDLE.name),(MANIFEST,'n83_basis_setup_HASH_MANIFEST.json'),(RECEIPT,'n83_basis_setup_archive_receipt.json')]:shutil.copy2(src,Path('/workspace/report_output')/name)
print(json.dumps(receipt),flush=True)
