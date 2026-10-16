from pathlib import Path
import hashlib,json
repo=Path('/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015')
p=repo/'research/isogeny_world_frontier_20261016/replay-optimization-v1/custody-approved-20261016'
inputs=json.loads((p/'inputs.json').read_text());idx=json.loads((p/'store-index.json').read_text())
assert (p/'transfer.exit').read_text().strip()=='0'
# Pull and verify are documented text outputs; seal them as stdout files.
for f in (p/'pull').glob('*.json'):f.rename(f.with_suffix('.stdout'))
if (p/'remote-verification.json').exists(): (p/'remote-verification.json').rename(p/'remote-verification.stdout')
expected={r['store_name']:r for r in inputs['inputs']}
manifest_bytes=(p/'inputs.json').read_bytes()
manifest_name=inputs['logical_prefix']+'inputs.json'
expected[manifest_name]={'id':'input-manifest','path':str(p/'inputs.json'),'sha256':hashlib.sha256(manifest_bytes).hexdigest(),'bytes':len(manifest_bytes),'store_name':manifest_name}
assert set(idx['entries'])==set(expected)
readbacks=[]
for name,row in sorted(expected.items()):
 entry=idx['entries'][name];assert entry['sha256']==row['sha256'] and entry['bytes']==row['bytes']
 pushed=json.loads((p/'push'/f"{row['id']}.json").read_text());assert pushed['sha256']==row['sha256'] and pushed['bytes']==row['bytes'] and not pushed['replaced']
 downloaded=Path((p/'pull'/f"{row['id']}.stdout").read_text().strip());assert downloaded.is_file()
 b=downloaded.read_bytes();assert hashlib.sha256(b).hexdigest()==row['sha256'] and len(b)==row['bytes']
 assert b==Path(row['path']).read_bytes()
 readbacks.append({'name':name,'sha256':row['sha256'],'bytes':len(b),'downloaded_path':str(downloaded),'byte_equal_to_frozen_input':True,'python_hashlib_sha256_matches_native_store_digest':True,'s3_uri':'s3://crypto-autoresearcher/store/sha256/'+row['sha256'][:2]+'/'+row['sha256']})
lines=(p/'remote-verification.stdout').read_text().splitlines()
assert len(lines)==len(expected)
assert set(lines)=={n+': cache ok, remote ok' for n in expected}
progress=[]
for f in p.rglob('*.stderr'):
 text=f.read_text().strip()
 if not text: continue
 if f.parent.name=='pull':
  assert text.startswith('fetched ') and ', sha256 ' in text and text.endswith(')'), (f,text)
  progress.append(f.name)
 else: assert f.name=='remote-verification.stderr' and text=='9 entries verified', (f,text)
assert len(progress)==8
native_uploads=sum(json.loads(f.read_text())['uploaded'] for f in (p/'push').glob('*.json'))
result={'schema':'isogeny-approved-custody-verification/v1','status':'VERIFIED_REMOTE_AND_FULL_READBACK','task_id':'TASK-20261016-de577e','remote':inputs['remote'],'logical_prefix':inputs['logical_prefix'],'native_source_commit':inputs['native_source_commit'],'authorization':'Explicit user approval received before this transfer; earlier rejected attempt executed zero upload operations.','native_transfer_exit':0,'logical_push_commands':9,'native_uploads_reported':native_uploads,'native_fresh_fetches_reported':len(progress),'stderr_classification':'Only expected native fetch/verified progress; complete raw text retained','logical_names':len(expected),'distinct_content_objects':len({r['sha256'] for r in expected.values()}),'full_readback_hashes':'PASS','all_readback_bytes_equal_originals':True,'native_head_size_checksum_verification':'PASS','native_client_sha256':(p/'store-client.sha256').read_text().split()[0],'native_client_build_source':'not independently established; transport client hash bound, remote contents independently rehashed','compiled_scientific_executables_uploaded':4,'scientific_runs':0,'started_at_utc':(p/'transfer-started-utc.txt').read_text().strip(),'ended_at_utc':(p/'transfer-ended-utc.txt').read_text().strip(),'objects':readbacks,'source_gate_complete':False,'independent_source_review_started':False}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','logical_names','distinct_content_objects','compiled_scientific_executables_uploaded','scientific_runs']}))
