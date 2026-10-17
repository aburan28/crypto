from pathlib import Path
import json,hashlib,re,subprocess,os
repo=Path('/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015')
root=Path('/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/clean-retake-v1')
source='015fa7bb504ef63953e0e04e464bb266f4f207bf'
assert (root/'completed.txt').read_text().strip()=='completed'
assert subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD']).decode().strip()==source
assert not subprocess.check_output(['git','-C',str(repo),'status','--porcelain=v1'])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_bytes())
arms={}
for name,policy,features in [('reference','reference-rem-v1',['--no-default-features']),('candidate','cached-fixed-modulus-v1',['--no-default-features','--features','cached-fixed-modulus-reduction'])]:
 a=root/name
 assert (a/'feature-arguments.txt').read_text().splitlines()==features
 for p in ['build.exit','tests.exit','monitor-isogeny.exit','monitor-fixed-modulus-bench.exit','timeout-isogeny.exit','timeout-fixed-modulus-bench.exit']:
  assert (a/p).read_text().strip()=='0',(name,p)
 for p in ['git-status-before.txt','git-status-after.txt']:assert (a/p).read_bytes()==b''
 version=read(a/'isogeny-version.json');bench=read(a/'fixed-modulus-bench-version.json')
 assert version['tool']['git_commit']==source and version['tool']['git_dirty']=='false'
 assert bench['git_commit']==source and bench['git_dirty']=='false'
 assert version['compiled_reduction_policy']==policy and bench['compiled_policy']==policy
 binaries={}
 for binary in ['isogeny','fixed-modulus-bench']:
  installed=a/'install'/binary;compiled=a/'target/release'/binary
  assert installed.read_bytes()==compiled.read_bytes()
  binaries[binary]={'installed_path':str(installed),'compiled_path':str(compiled),'sha256':sha(installed),'bytes':installed.stat().st_size,'mode':oct(installed.stat().st_mode&0o777),'compiled_equals_installed':True}
 log=(a/'tests.stdout').read_text()
 counts=[tuple(map(int,m)) for m in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',log)]
 assert len(counts)==8 and sum(x[0] for x in counts)==128 and sum(x[1] for x in counts)==0
 for binary in ['isogeny','fixed-modulus-bench']:
  monitor=read(a/'controls'/f'monitor-{binary}'/'compact-cm-before-missing/native.receipt.json')
  paths=list((a/'controls'/f'monitor-{binary}').glob('*/native.receipt.json'))
  assert len(paths)==8
  assert all(read(p)['status']=='RESOURCE_MONITOR_UNAVAILABLE' for p in paths)
  for kind in ['monitor','timeout']:
   selected=(a/f'{kind}-{binary}.stdout').read_text()
   assert '1 passed; 0 failed;' in selected,(name,binary,kind,selected)
  timeout=read(a/'controls'/f'timeout-{binary}'/'native.receipt.json')
  assert timeout['status']=='TIMEOUT'
 arms[name]={'features':features,'policy':policy,'version':version,'companion_version':bench,'binaries':binaries,
   'native_counts':counts,'native_tests_passed':128,'monitor_injection_receipts':16,'deadline_receipts':2}
comparisons=[]
for label,relative in [('P192/13','challenge/p192/judgement/verification.json'),('P521/7','challenge/p521/judgement/verification.json'),
                       ('P224/1471','cache/valid-judge/verification.json'),('P192/13 trace off','trace/valid-off/verification.json'),
                       ('P192/13 trace on','trace/valid-on/verification.json'),('P192/13 trace second','trace/valid-second/verification.json')]:
 a=root/'reference/controls'/relative;b=root/'candidate/controls'/relative
 assert a.read_bytes()==b.read_bytes(),label
 v=read(a)
 assert v['status']=='PASS'
 assert all(r['public_scalar_transport_checks']==20 for r in v['fresh_certificate']['records'])
 comparisons.append({'case':label,'reference':str(a.relative_to(root)),'candidate':str(b.relative_to(root)),'sha256':sha(a),'bytes':a.stat().st_size,'literal_equal':True,'excluded_fields':[]})
rejects={}
for name in ['reference','candidate']:
 data=read(root/name/'controls/cache/controls.json')
 rows=[{'case':r['case'],'status':r['status'],'exit':r['exit'],'reason':r['entry'].get('reason'),'score':r['entry'].get('score')} for r in data['controls']]
 trace=[]
 for p in sorted((root/name/'controls/trace').glob('command-*.json')):
  c=read(p)
  try:v=json.loads(c['stdout'])
  except (KeyError,json.JSONDecodeError):continue
  argv=c['args']
  if len(argv)>3 and argv[:2]==['challenge','judge'] and '/bad-' in argv[3]:
   trace.append({k:v.get(k) for k in ['status','reason','score']})
 assert len(trace)==18
 rejects[name]={'cache':rows,'adversarial_trace_pairs':trace}
assert rejects['reference']==rejects['candidate']
paths=['research/p192_p224_large_degree_isogenies_20261009/verification_poly.rs','research/p192_p224_large_degree_isogenies_20261009/verification_kernel.rs',
 'tools/isogeny-cli/src/replay.rs','tools/isogeny-cli/Cargo.toml','tools/isogeny-cli/Cargo.lock','tools/isogeny-cli/src/main.rs','tools/isogeny-cli/src/supervisor.rs',
 'tools/isogeny-cli/build.rs','tools/isogeny-cli/src/bin/fixed_modulus_bench.rs','tools/isogeny-cli/tests/fixed_modulus_replay.rs','docs/identity-certificates/fixed-modulus-newton.json',
 'docs/identity-certificates/fixed-modulus-newton.md','src/cryptanalysis/identity_certificate.rs','docs/curves/registry.json']
files={str(p.relative_to(root)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(root.rglob('*')) if p.is_file() and 'target' not in p.relative_to(root).parts and 'install' not in p.relative_to(root).parts}
record={'schema':'fixed-modulus-clean-build-custody-input/v1','repository':'aburan28/crypto','source_commit':source,
 'source_tree':subprocess.check_output(['git','-C',str(repo),'rev-parse',source+'^{tree}']).decode().strip(),
 'source_files':{p:{'sha256':sha(repo/p),'bytes':(repo/p).stat().st_size} for p in paths},'arms':arms,
 'mathematical_outputs':comparisons,'rejection_comparison':rejects['reference'],'mathematical_fields_excluded':[],
 'raw_files':files,'scientific_runs':0,'performance_claim':False,'independent_review':'pending'}
with (root/'manifest.json').open('x') as f:json.dump(record,f,indent=2)
print(json.dumps({'source':source,'arms':{n:{'sha256':a['binaries']['isogeny']['sha256'],'native_passed':128} for n,a in arms.items()},'math_pairs_equal':len(comparisons),'raw_file_count':len(files)}))

