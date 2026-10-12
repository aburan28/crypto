"""Offline public-command audit checker. Never starts/polls jobs or reads live logs."""
from pathlib import Path
import argparse,csv,hashlib,json,math,re,shlex,statistics
ROOT=Path('/private/tmp/ecc2k-shared-sigma-public-src-20260912')
FREEZE=Path('/private/tmp/ecc2k-shared-sigma-public-freeze-20260912.json')
MEASURED=Path('/private/tmp/ecc2k-shared-sigma-walk-20260912/source')
OUT=Path('/Users/adamburan/Documents/Codex/2026-09-08/pl/outputs')
HEAD='6a48e3ec6dbd51a54e22cace4c24b27e0022b3d0'
FROZEN_FILE_COUNT=251
REMOTE_FILE_COUNT=71
FREEZE_SHA='bbb3029be60601874e31acbdf9d12a7e94b6436da5253ae3ae4a088f6dd74ee1'
SOURCE_AGGREGATE='837296393f447dd524bfa1cfc6066a3d0c1375dd8fa31ee336a1a9eb135f891d'
UPDATES=201863462912;WORKERS=385024;BATCH=16;WALKS=6160384
# Owner-specified known-answer expectation for this fixed workload, not the general artifact schema.
EXPECTED_COLLECTION_RECORDS=5149;EXPECTED_COLLECTION_BYTES=164768
# Explicit measured-runtime bindings; these prove source identity, not linked code identity.
CORE={'include/bigmod.h': 'a71b62074d426b71f255ec664eb8c6921e1235664b26fabe0425315a1eba45c5', 'include/bitslice.h': 'f8cf0d47abe28080fac82a7c26dfe1a03debd7b7239a94142249d72b96278b34', 'include/curveparams.h': 'ee8ffb38e12d546e2fecea63b2678879388994f526d36a420d929a97f8c3085b', 'include/fieldbs.h': '49313d9e3fd4214c3b618ca1fbe8c0c8d02d1c7c2ac119dfe8fe5c0f992a65cf', 'include/fieldpb.h': '9e9fddbefa90c734f88aa97dcfa785a08cc71d539926df3f491bcf004733b919', 'include/kernel.h': '5186e0c44830856bdc6c692faaf81c055419dd3f549c03d0127ce99af5bdcbcf', 'include/packed131.h': '930ffa60f1b42ac14700ec94b42294c228bf449501ab4accef4097d0b9762ca0', 'include/packedcompactstate.cuh': '1f2bb0215807ef3dc0880c417c515185da9200ee4a8fcc983759cb57309249c4', 'include/packeddirectreduce131.h': '8a25b93b3d271b84fbc8720f7ba074a58c0b560eb27ca8925344effafdad07e2', 'include/packedengine.cuh': '054645487d1c39a6daf9ac0e4de537eca420c889b7a7958a6c79b0c8db6a3496', 'include/packedgeneratedproduct131.h': '8561eed5899fd58505e8bbf9ea96aab2fea5cd0d159d56e6c3afa47e622801ea', 'include/packedkernels.cuh': '19bc960f45cba4a6aa45188fb9d921be6682d0b384d590593f5974e6b0d16249', 'include/packedpolyreduce131.h': '0bbef59761d3f95544a4120e843b9d602aef464f8a770ab3d4cd1e31dee92730', 'include/packedsigma131.h': 'a3276bbd3c1cd1467b5dba3afad1dcbd29813b390f5847b032c5b3a91e03a623', 'include/packedtransform131.h': '06e3a2a473d44c5d40393eed703725f9aa16d2d02d024ca3fbbfa0d05d0e71ee', 'include/ref.h': 'e48e9167c817b72e16c5c3f0b2ca094deb7c38a037157c6d95b15796b01531c8', 'include/solver.h': 'a2d23d5a10a588005ee01fb956fff98925630d3452f3cf1e41a1ff5fa611877c', 'include/walk.h': '0af88d29695599f401e6eb4a0cc660282142b158b647e0ec0407e39ecd6b0582', 'src/main.cu': '1aacf9eba08a76d825434d081caa998ca4d615065b21006b1c8d33b940ffd2ab', 'src/profile_probe.cu': '8a86edf891d8bddb778b8014f8adfb33afa5195aff7fd2ffb79d7530cdd211b7', 'src/ptxspike.cu': 'd39ffbf847052af4fe9f0eccb8179f04ace908a470d133317eda70cc8b9f99ad', 'src/testpacked.cpp': '07fad5ccb5214faa4d6d0b3445b6b504c5e34f21282b95adf8d45dd57c59f45f', 'src/testpackedcuda.cu': 'e5d09cd9260e92d5afc8b58328ed98e080b19e122f0a24323164497ef8124598', 'src/testpackedstatecuda.cu': '6c4956263002e12cefcb360fa6dfdbfe6fd5514b4de411358d20f3e982c24eba', 'src/testschedule.cpp': '79cc84736ff42ce55408ebd4a6e0dae9c507241acd9baee8a95039fc740ff99f', 'src/testsharedsigmacuda.cu': '0d1a642bd0ac66414ef6cf524f2366ea9a0848044d006e04bbaf87d0a26d7ff5', 'src/testtiming.cpp': '726aed66d869a11e680075bc54f0bd2d942542a85639d022b3cea96640f1c9a5', 'generated/eccF131.h': '64c4ca165fdfad360800bcffc8825e4e8de01789b6844649a8795cbfb437a9c2', 'generated/eccF23.h': 'e07a0685ee084dd99c250a1161ef6cf7c83092894c4b54e9f73b3b58c54bf5c8', 'generated/eccF41.h': '8eea324fec5fd0ff04267d9a2d38eb9ef2d5d296bcd29b011e4212036fcadb31', 'generated/eccF83.h': '886e1543a9a25d1633a353a7af23bd2b3ac2468c792d0991992dde8a3cdbc3f2', 'generated/eccP13.h': '4c7eb2be0286c4524622b3dec9ca0317e897ee06d6720c2bd641f434fcae4415', 'generated/eccP19.h': '95438543886989db7720daf7e63ed46c7d38d57a8ea6cdebd87879ea9b4bfbd5', 'generated/eccP41.h': 'e83923268b16a821b57a3de11b0c9ec58057424e88e5bb4f54f8a3962b323081', 'generated/eccP97.h': '8a8f9e6267f4d603d760d3a784d7936e52be1e723fd64fefb853cc476d6bd26f', 'generated/rtl/ecc_hamming131.v': '1edcfdace9f086f67e22d955d12982d2403672ade1cf81b00932814838e41454', 'generated/rtl/ecc_mul131.v': 'e39ef0ca0d4ef83efb3c4b5c3053894f8d971338fb373c18896923eca8ac5573', 'generated/rtl/ecc_sigma131.v': '16bdc5b629dd139370ce6e53018284dacf081d3b4c7ccc430d4a8f30f939291e', 'generated/rtl/tb_ecc_mul131.v': '20fcc12ac3e12f87857d11b8491741aa38d30a236dc043704b941a3a0bb48b63', 'generated/rtl/vectors_mul131.txt': 'ff4e79c857cfce97c9c9365f2e7e9c359d43ec08abeb71d9ca78c9ddaddba6a8', 'codegen/gensigma.py': '6360039e5ff0b2c31ff092c4eb2c1e006e5e681eb19cae88f35c6c1be5f3e0a7', 'codegen/genpackedproduct.py': '490760f7345e742b360bbbd3676a61cb9910ce43817c227365ec71e746fe8602'}
PASS_LINES=[
 'PASS: 3120 GPU Frobenius vectors, every field basis vector for all selected powers plus dense cases',
 'PASS: 2526 GPU polynomial reductions against long division, including ignored upper-word bits and canonical outputs',
 'PASS: 18194 GPU polynomial products, including all 17161 basis pairs',
 'PASS: 18194 GPU paired polynomial products against independent multiplication',
 'PASS: 1157 GPU polynomial squares against independent multiplication and long division',
 'PASS: 6240 GPU paired Frobenius vectors, both inputs against independent routing']
STORAGE_PASS='PASS: 128 GPU storage cases, 297344 records, independent physical images and logical reads with canaries'
SHARED_PASSES=[
 'PASS: 21 GPU sigma scenarios, 21036 input pairs, global and selected helpers against independent routing',
 'PASS: 114 complete block mask snapshots, 51072 words, output guards and inactive blocks']
INTEGRATION=[
 'PASS: DP replay across launch boundaries, restart and checkpoint resume',
 'PASS: byte-identical resume, scalar iteration count, incompatible checkpoint preserved',
 'PASS: overdue walks restart without false distinguished-point reports']
MARKERS={'denominator cache':1,'multiply by value':1,'Frobenius network':3,'polynomial chain':1,'polynomial state':1,'unrolled inversion':1,'paired products':1,'direct reduction':1,'generated product':1,'native carryless multiply':1,'weighted prefix':2,'compact state':1,'shared sigma':1,'state tile':256}
MAKE_FLAGS={'BATCH':'16','THREADS':'256','MINBLOCKS':'2','PACKED_SINGLE_PRODUCT':'1','PACKED_CACHE_DENOM':'1','PACKED_BY_VALUE':'1','PACKED_PERM_SIGMA':'3','PACKED_POLY_CHAIN':'1','PACKED_UNROLL_INV':'1','PACKED_PAIR_PRODUCTS':'1','PACKED_POLY_STATE':'1','PACKED_DIRECT_REDUCE':'1','PACKED_GENERATED_PRODUCT':'1','PACKED_CLMAD':'1','PACKED_COMPACT_STATE':'1','PACKED_SHARED_SIGMA':'1','PACKED_WEIGHTED_PREFIX':'2','PACKED_STATE_TILE':'256'}
def need(x,m):
 if not x:raise ValueError(m)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def digest(x):return isinstance(x,str) and re.fullmatch('[0-9a-f]{64}',x) is not None

def prepare():
 need(sha(FREEZE)==FREEZE_SHA,'committed freeze changed');f=read(FREEZE)
 need(f['commit']==HEAD and f['source']==str(ROOT) and f['fileCount']==FROZEN_FILE_COUNT and len(f['files'])==FROZEN_FILE_COUNT and FROZEN_FILE_COUNT>0,'source snapshot identity')
 need(all(not Path(p).is_absolute() and '..' not in Path(p).parts for p in f['files']),'snapshot path outside committed source')
 need(all(sha(ROOT/p)==h for p,h in f['files'].items()),'committed source bytes changed')
 paths=[ROOT/'Makefile',ROOT/'modal_app.py']
 for folder in ('include','src','generated','codegen'):
  paths+=sorted(p for p in (ROOT/folder).rglob('*') if p.suffix in ('.h','.cuh','.cu','.cpp','.py'))
 need(len(paths)==REMOTE_FILE_COUNT and len(set(paths))==REMOTE_FILE_COUNT and REMOTE_FILE_COUNT>0 and all(str(p.relative_to(ROOT)) in f['files'] for p in paths),'exact remote source aggregate coverage')
 h=hashlib.sha256()
 for p in paths:h.update(str(p.relative_to(ROOT)).encode()+b'\0'+p.read_bytes()+b'\0')
 need(h.hexdigest()==SOURCE_AGGREGATE,'local source aggregate changed')
 need(CORE and all(sha(ROOT/p)==sha(MEASURED/p)==v for p,v in CORE.items()),'measured core/test source identity')
 make=(ROOT/'Makefile').read_text();need('RTX_PRO6000_COMPACT_STATE ?= 1' in make and 'RTX_PRO6000_SHARED_SIGMA ?= 1' in make and 'ECC_PACKED_SHARED_SIGMA=$(RTX_PRO6000_SHARED_SIGMA)' in make and 'ECC_PACKED_WEIGHTED_PREFIX=2 ECC_PACKED_COMPACT_STATE=$(RTX_PRO6000_COMPACT_STATE)' in make and '--batch 16 --min-blocks 2 --block-threads 256 --repeats 3 --workers $(RTX_PRO6000_WORKERS)' in make,'public preset/opt-out source')
 wrapper=(ROOT/'packed_audit.py').read_text();need(wrapper.index('storage = result["deviceStorage"] = run')<wrapper.index('probe = result["deviceSharedSigma"] = run')<wrapper.index('result["integration"] = run')<wrapper.index('result["benchmark"] = client.measureBench'),'source storage/probe/integration/timing order')
 return dict(commit=HEAD,frozenSourceFiles=FROZEN_FILE_COUNT,remoteSourceFiles=REMOTE_FILE_COUNT,remoteSourceSha256=h.hexdigest(),aggregatePaths=[str(p.relative_to(ROOT)) for p in paths],measuredSourceMatches=CORE,wrapperOrderingSourceBound=True)

def terminal_context(owner):
 need(owner.get('confirmed') is True,'owner confirmation of terminal result required before reading artifact')
 need(type(owner.get('session')) is int and owner['session']>0 and type(owner.get('exitCode')) is int and owner['exitCode']==0,'expected owner terminal session/exit')
 need(isinstance(owner.get('appId'),str) and re.fullmatch(r'ap-[A-Za-z0-9]+',owner['appId']) is not None,'owner-confirmed app ID required')

def component_command(row,target):
 need(row['returncode']==0 and not row.get('error') and row.get('timedOut') in (None,False),'component did not complete')
 argv=row['command'];need(isinstance(argv,list) and argv[:2]==['make',target],'component make target')
 assignments={}
 for token in argv[2:]:
  need(isinstance(token,str) and '=' in token,'unexpected component argument');k,v=token.split('=',1);need(k not in assignments,'duplicate component setting');assignments[k]=v
 need(assignments==dict(ARCH='-gencode arch=compute_120,code=sm_120',**MAKE_FLAGS),'component architecture/flags differ')

def arithmetic(row):
 component_command(row,'test-packed-cuda');lines=row['output'].splitlines()
 for label,value in [('direct reduction',1),('generated product',1),('native carryless multiply',1),('weighted prefix',2)]:
  need(re.findall(r'^packed arithmetic '+re.escape(label)+r': (.*)$',row['output'],re.M)==[str(value)],'arithmetic marker '+label)
 need([x for x in lines if x.startswith('PASS:')]==PASS_LINES,'all six exact arithmetic suites')
 for k in ['packedGeneratedProduct','expectedPackedGeneratedProduct','packedClmad','expectedPackedClmad']:need(row[k] is True,'arithmetic boolean identity')
 need(row['packedWeightedPrefix']==row['expectedPackedWeightedPrefix']==2,'arithmetic weighted identity')

def storage(row,arithmetic_command):
 component_command(row,'test-packed-storage-cuda')
 wanted=list(arithmetic_command);wanted[1]='test-packed-storage-cuda';need(row['command']==wanted,'same storage build settings')
 need(row['packedCompactState'] is True and row['expectedPackedCompactState'] is True,'storage compact boolean identity')
 need(re.findall(r'^packed storage compact state: (.*)$',row['output'],re.M)==['1'] and re.findall(r'^packed storage batch: (.*)$',row['output'],re.M)==['16'],'storage mode/batch transcript')
 need([x for x in row['output'].splitlines() if x.startswith('PASS:')]==[STORAGE_PASS] and row['cases']==128 and row['records']==297344,'complete storage validation')

def shared_probe(row,arithmetic_command):
 component_command(row,'test-shared-sigma-cuda')
 wanted=list(arithmetic_command);wanted[1]='test-shared-sigma-cuda';need(row['command']==wanted,'same shared probe build settings')
 need(row['packedSharedSigma'] is True and row['expectedPackedSharedSigma'] is True,'probe shared boolean identity')
 need(re.findall(r'^packed shared sigma probe: (.*)$',row['output'],re.M)==['1'],'shared probe mode transcript')
 need([x for x in row['output'].splitlines() if x.startswith('PASS:')]==SHARED_PASSES,'both exact shared probe suites')
 for key,value in [('scenarios',21),('inputPairs',21036),('blockSnapshots',114),('maskWords',51072)]:
  need(type(row[key]) is int and row[key]==value,'shared probe count '+key)

def sample(row,weight):
 need(row['returncode']==0 and row['valid'] is True and not row.get('error') and row.get('timedOut') in (None,False),'sample completion')
 for k in ['packedGeneratedProduct','expectedPackedGeneratedProduct','packedClmad','expectedPackedClmad','packedCompactState','expectedPackedCompactState','packedSharedSigma','expectedPackedSharedSigma','packedDirectReduction','expectedPackedDirectReduction']:need(row[k] is True,'sample boolean identity '+k)
 need(row['packedWeightedPrefix']==row['expectedPackedWeightedPrefix']==2 and row['packedStateTile']==row['expectedPackedStateTile']==256,'sample mode identity')
 exact={'reportedIterations':UPDATES,'expectedIterations':UPDATES,'requestedBatch':16,'actualBatch':16,'requestedWorkers':WORKERS,'actualWorkers':WORKERS,'scalarWalks':WALKS}
 need(all(type(row[k]) is int and row[k]==v for k,v in exact.items()),'exact sample geometry/count metadata')
 argv=shlex.split(row['command']) if isinstance(row['command'],str) else row['command'];need(isinstance(argv,list),'sample command type')
 if weight==0:
  need(argv==['./ecc2k130','--curve','131','--bench','--steps','1024','--launches','32','--verify','0','--packed','--threads',str(WORKERS)],'exact benchmark argv')
 else:
  prefix=['./ecc2k130','--packed','--curve','131','--dp-weight','34','--steps','1024','--launches','32','--verify','0','--dp-file'];need(argv[:len(prefix)]==prefix and len(argv)==len(prefix)+3 and argv[-2:]==['--threads',str(WORKERS)] and isinstance(argv[len(prefix)],str) and argv[len(prefix)],'exact collection argv')
 raw=row['raw'];need(not any(s in raw for s in ['MISMATCH','stopping:','warning: could not write checkpoint']),'interrupted or failed sample')
 for label,value in MARKERS.items():need(re.findall(r'^packed '+re.escape(label)+r': (.*)$',raw,re.M)==[str(value)],'runtime marker '+label)
 need(re.findall(r'^packed kernel: (.*)$',raw,re.M)==['104 registers/thread, 0 local bytes/thread, 1792 shared bytes/block, single-product multiplier'],'runtime kernel resources')
 need(re.findall(r'^packed driver reserved shared bytes/block: (.*)$',raw,re.M)==['1024, device 0'],'separate driver shared reservation')
 need(re.findall(r'^backend cuda-packed131: (.*)$',raw,re.M)==[f'{WORKERS} threads x 16 slots x 1 lanes = {WALKS} walks, dp weight {weight}, 1024 steps per launch'],'runtime scalar geometry')
 counts=list(map(int,re.findall(r'(\d+) iterations\s+\d+ dp',raw)));need(counts and counts[-1]==UPDATES and all(0<x<=UPDATES and x%(WALKS*1024)==0 for x in counts) and all(a<b for a,b in zip(counts,counts[1:])),'monotone complete scalar progress on launch grid')
 progress=[tuple(map(int,x)) for x in re.findall(r'M it/s\s+(\d+) iterations\s+(\d+) dp\s+(\d+) stored\s+(\d+) dropped',raw)]
 need([x[0] for x in progress]==counts and all(x[3]==0 and 0<=x[2]<=x[1] for x in progress),'complete progress accounting and zero drops')
 need(all(a[1]<=b[1] and a[2]<=b[2] for a,b in zip(progress,progress[1:])),'monotone fresh report/storage counters')
 final=re.findall(r'^\s*finished:\s+(\S+) M it/s, (\d+) distinguished points \((\d+) verified against the reference, (\d+) dropped\)$',raw,re.M)
 need(len(final)==1 and tuple(map(int,final[0][1:]))==((EXPECTED_COLLECTION_RECORDS if weight else 0),0,0),'final reports/drops')
 need(progress[-1][1:]==(int(final[0][1]),int(final[0][1]),0),'final fresh progress agrees with reports and stored count')
 rate=float(final[0][0]);need(math.isfinite(rate) and rate>0 and math.isfinite(row['rate']) and row['rate']==rate,'finite final rate')
 if weight:need(row['corpusRecords']==EXPECTED_COLLECTION_RECORDS and row['corpusBytes']==EXPECTED_COLLECTION_BYTES,'public collection counts/sizes')
 return dict(rateM=rate,rateB=rate/1000,updates=UPDATES,reports=int(final[0][1]),dropped=0,compactState=True,sharedSigma=True,weightedPrefix=2)

def evaluate(record,source,owner):
 terminal_context(owner)
 r=record;need(r['valid'] is True and not r.get('error'),'public audit failed or incomplete')
 for k in ['packedGeneratedProduct','expectedPackedGeneratedProduct','packedClmad','expectedPackedClmad','packedCompactState','expectedPackedCompactState','packedSharedSigma','expectedPackedSharedSigma','packedDirectReduction']:need(r[k] is True,'top-level boolean '+k)
 need(r['packedWeightedPrefix']==r['expectedPackedWeightedPrefix']==2 and r['packedStateTile']==r['expectedPackedStateTile']==256,'top-level layout/schedule')
 need(tuple(r[k] for k in ['batch','blockThreads','minBlocks','requestedWorkers','steps','launches','repeats'])==(16,256,2,WORKERS,1024,32,3),'public command geometry')
 identity=r['identity'];need(identity['sourceSha256']==source['remoteSourceSha256'] and digest(identity['binarySha256']),'reported source and binary identity')
 for k in ['packedDenominatorCache','packedByValue','packedPolynomialChain','packedUnrolledInverse','packedPairedProducts','packedPolynomialState','packedDirectReduction','packedGeneratedProduct','packedClmad','packedCompactState','packedSharedSigma']:need(identity[k] is True,'identity '+k)
 need(identity['packedWeightedPrefix']==2 and identity['packedStateTile']==256 and identity['packedFrobeniusMode']==3 and identity['activeBackend']=='packed-poly131' and identity['packedMultiplier']=='single-product','packed backend identity')
 need(identity['compilerReturncode']==0 and 'V13.3.73' in identity['compiler'] and identity['cudaImageVersion']=='13.3.1' and identity['gpuStateReturncode']==0,'compiler/image/inventory status')
 gpu=[[v.strip() for v in line] for line in csv.reader(identity['gpuState'].splitlines()) if line];need(len(gpu)==2 and len(gpu[0])==len(gpu[1])==8 and gpu[0][:3]==['name','uuid','driver_version'],'single GPU snapshot')
 need(gpu[1][0]=='NVIDIA RTX PRO 6000 Blackwell Server Edition' and re.fullmatch(r'GPU-[0-9a-fA-F]{8}(?:-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}',gpu[1][1]),'GPU identity')
 arithmetic(r['deviceArithmetic']);storage(r['deviceStorage'],r['deviceArithmetic']['command'])
 need(r['sharedSigmaProbeApplicable'] is True,'selected preset must execute shared probe')
 shared_probe(r['deviceSharedSigma'],r['deviceArithmetic']['command'])
 integ=r['integration'];need(integ['returncode']==0 and not integ.get('error') and integ.get('timedOut') in (None,False) and integ['command']==['python3','codegen/testpackedclient.py','./ecc2k130'] and integ['output'].splitlines()==INTEGRATION,'complete integration and exact transcript')
 need(r['collection']==r['collectionSummary']['samples'] and len(r['collection'])==3,'collection summary identity')
 summaries={};samples=[]
 for phase,key,weight in [('benchmark','benchmark',0),('collection','collectionSummary',34)]:
  summary=r[key];need(summary['valid'] is True and len(summary['samples'])==3,'three complete repetitions')
  checked=[sample(x,weight) for x in summary['samples']];rates=[x['rateM'] for x in checked]
  need(summary['rate']==statistics.median(rates) and summary['minRate']==min(rates) and summary['maxRate']==max(rates),'recomputed public statistics')
  summaries[phase]=dict(ratesM=rates,medianM=statistics.median(rates),minimumM=min(rates),maximumM=max(rates),medianB=statistics.median(rates)/1000)
  samples += [dict(phase=phase,repeat=i,**x) for i,x in enumerate(checked)]
 need(len(samples)==6,'six exact samples')
 return dict(valid=True,gitHead=HEAD,sourceFreezeSha256=FREEZE_SHA,frozenSourceFiles=FROZEN_FILE_COUNT,remoteDigestFiles=REMOTE_FILE_COUNT,sourceSha256=source['remoteSourceSha256'],sourceAggregatePaths=source['aggregatePaths'],measuredRuntimeSourceMatches=source['measuredSourceMatches'],reportedLinkedBinarySha256=identity['binarySha256'],compiler='13.3.73',cudaImageVersion='13.3.1',gpuSnapshot=dict(zip(gpu[0],gpu[1])),ownerTerminal=owner,compactState=True,weightedPrefix=2,clmad=True,tile=256,batch=16,workers=WORKERS,scalarWalks=WALKS,registers=104,localBytes=0,sharedBytes=1792,driverReservedSharedBytes=1024,sharedSigma=True,arithmeticPassLines=PASS_LINES,storageCases=128,storageRecords=297344,correctnessGatesPassed=True,sourceOrderedStorageBeforeIntegrationAndTiming=True,samples=samples,summaries=summaries,corpusCountsAndSizesOnly=dict(collectionSamples=3,recordsEach=EXPECTED_COLLECTION_RECORDS,bytesEach=EXPECTED_COLLECTION_BYTES,zeroDrops=True,ownerSpecifiedWorkloadRegression=True),goalEvaluation=dict(targetBillionScalarUpdatesPerSecond=26,benchmarkMedianAtLeast26B=summaries['benchmark']['medianM']>=26000,allThreeBenchmarkSamplesAtLeast26B=summaries['benchmark']['minimumM']>=26000,collectionMedianAtLeast26B=summaries['collection']['medianM']>=26000),findings=[],limits=[
 'Separate-allocation public-command audit; this does not establish another paired gain or inherit the private comparison rates.',
 'Owner-supplied terminal session/app/exit are recorded; this offline checker does not independently query the scheduler or read the execution log.',
 'The committed-source snapshot and exact source aggregate are verified. The aggregate excludes root packed_audit.py, whose bytes are retained in the committed snapshot.',
 'The pinned core/test source files match the measured comparison runtime. The public artifact reports a pre-validation linked-binary hash but has no post-run binary/source hashes or full SASS, so whole-binary/native-code identity with the private measurement is not established here.',
 'Storage-before-shared-probe-before-integration-before-timing order is verified in the frozen producer source; the artifact has no independent event chronology.',
 'The checker additionally pins the owner-specified5149-record/164768-byte expectation for this exact workload. The general public schema accepts other positive consistent counts. Content hashes/payloads are not retained, so corpus multiset equality is not independently established here.',
 'GPU telemetry is one pre-validation snapshot, not sustained clocks, achieved occupancy, cache traffic or a hardware ceiling. The public artifact does not retain the JIT environment.',
 'No compiler, client, GPU, remote job, live-log read or numerical rerun is performed by this checker.'])

def audit(path,expected_sha,owner):
 terminal_context(owner)
 need(digest(expected_sha),'terminal artifact SHA required')
 source=prepare();need(sha(path)==expected_sha,'owner terminal artifact SHA mismatch')
 result=evaluate(read(path),source,owner);result.update(rawSha256=expected_sha,auditorSha256=sha(__file__))
 return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepare-only',action='store_true');p.add_argument('--artifact',type=Path);p.add_argument('--sha256');p.add_argument('--owner-confirmed-terminal',action='store_true');p.add_argument('--session',type=int);p.add_argument('--app');p.add_argument('--exit-code',type=int);args=p.parse_args()
 if args.prepare_only:
  source=prepare();print(json.dumps(dict(valid=True,auditorSha256=sha(__file__),**source),indent=2));return
 owner=dict(confirmed=args.owner_confirmed_terminal,session=args.session,appId=args.app,exitCode=args.exit_code)
 result=audit(args.artifact or ROOT/'build/rtx-pro6000-audit.json',args.sha256,owner)
 path=OUT/'shared-sigma-public-native-review.json';path.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(valid=True,receipt=str(path),sha256=sha(path),summaries=result['summaries'],goalEvaluation=result['goalEvaluation']),indent=2))
if __name__=='__main__':main()
