#!/usr/bin/env python3
"""Freeze, run, independently check, and select complete IC candidate experiments."""
import argparse
import copy
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import resource
import shutil
import signal
import statistics
import subprocess
import sys
import time

from oracle import Curve, InvalidEvidence, require, verify

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
UNIT = 'valgrind-3.22-amd64-Ir'
STAGES = ['aa', 'smoke', 'development', 'selection', 'confirmation', 'replay']
BASE_CONFIG = {'solver':'pair_table','linear_algebra':'sparse','batch_trials':64,
               'max_trials':4096,'summands':3}
PHASES = {'startup_and_input','curve_and_targets','factor_base_and_tables',
          'collection_and_decomposition','verify_filter_and_linear_algebra',
          'log_certification','individual_log','final_verification',
          'reporting_and_cleanup','rho_solve'}


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for data in iter(lambda:f.read(1024*1024), b''):
            h.update(data)
    return h.hexdigest()


def objhash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def write(path, value, *, exclusive=False):
    path = Path(path)
    data = json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'
    path.parent.mkdir(parents=True,exist_ok=True)
    if exclusive:
        with path.open('x') as f:
            f.write(data)
    else:
        temporary = path.with_name(path.name+'.tmp')
        temporary.write_text(data)
        temporary.replace(path)


def candidates():
    changes = [
        ('batch16', {'batch_trials':16}, 'Reduce surplus collection and verification before a rank check.'),
        ('batch32', {'batch_trials':32}, 'Balance surplus relations against repeated rank checks.'),
        ('batch128', {'batch_trials':128}, 'Amortize rank checks when relation yield is low.'),
        ('dense', {'linear_algebra':'dense'}, 'Avoid sparse-solver setup on small relation matrices.'),
        ('excess8', {'sparse':{'filter':{'target_excess':8}}}, 'Reduce the filtered core by retaining fewer surplus rows.'),
        ('window8', {'collection_window':8}, 'Trade extra walked probes for fewer pair-table scans.'),
    ]
    out = [{'id':'incumbent','parent':None,'hypothesis':'Frozen complete library pipeline.',
            'config':copy.deepcopy(BASE_CONFIG)}]
    for name, delta, hypothesis in changes:
        out.append({'id':name,'parent':'incumbent','hypothesis':hypothesis,
                    'falsification':'Fails correctness, total-cost or per-cell confirmation gate.',
                    'changed_parameters':delta,'config':dict(BASE_CONFIG,**delta)})
    return out


def proposals_from_previous(root):
    result=read(root/'decision.json')
    require(result.get('status') in ('promoted','retained') and result.get('winner'),
            'finish and audit the previous decision before generating a successor')
    previous=read(root/'candidates.json')
    parent=next(a for a in previous if a['id']==result['winner'])
    config=copy.deepcopy(parent['config'])
    baseline={'id':'incumbent','parent':parent['id'],'parent_round':str(root),
              'hypothesis':'Retain the previous selected complete algorithm.','config':config}
    seen={objhash(a['config']) for a in previous}
    batch=config.get('batch_trials',64)
    changes=[]
    for value in [batch//2,batch//4,batch*2,batch*4,1,2,8,256]:
        if 1<=value<=min(config.get('max_trials',4096),4096):
            changes.append((f'batch{value}',{'batch_trials':value},'Test the next collection/rank-check tradeoff around the selected batch.'))
    mode='dense' if config.get('linear_algebra','sparse')=='sparse' else 'sparse'
    changes.append((mode,{'linear_algebra':mode},'Test linear-algebra setup against the selected collection policy.'))
    for value in [4,16,32]:
        changes.append((f'window{value}',{'collection_window':value},'Test a new walked-probe scan window with all failed attempts charged.'))
    out=[baseline]
    for name,delta,hypothesis in changes:
        cfg=dict(config,**delta)
        key=objhash(cfg)
        if key in seen:
            continue
        seen.add(key)
        out.append({'id':name,'parent':'incumbent','parent_round':str(root),
                    'hypothesis':hypothesis,'changed_parameters':delta,'config':cfg,
                    'falsification':'Fails the unchanged complete-workload confirmation gate.'})
        if len(out)==7:
            break
    require(len(out)>1,'no distinct configured candidates remain')
    return out,root/parent.get('source_directory','source')


def child_env():
    # Numerical/runtime knobs and preloaded libraries must not leak into a run.
    keep = ('PATH','LANG','LC_ALL','TZ','CARGO_HOME','RUSTUP_HOME','HOME')
    env = {k:os.environ[k] for k in keep if k in os.environ}
    env.update(RAYON_NUM_THREADS='1',IC_ARTIFACT_CACHE='off',IC_F2_BACKEND='cpu',
               PYTHONHASHSEED='0')
    return env


def execute(command, job, directory, timeout, memory, cpu):
    directory.mkdir(parents=True,exist_ok=True)
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(memory,memory))
        resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        os.sched_setaffinity(0,{cpu})
    start = time.monotonic()
    status = 'EXITED'
    with (directory/'stdout.json').open('w') as stdout, (directory/'stderr.txt').open('w') as stderr:
        process = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=stdout,stderr=stderr,
            env=child_env(),start_new_session=True,preexec_fn=limits)
        try:
            process.communicate(json.dumps(job,sort_keys=True).encode(),timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL)
            process.communicate()
            status = 'TIMEOUT'
    return {'exit_code':process.returncode,'process_wall_seconds':time.monotonic()-start,
            'process_status':status,'command':command,'memory_cap_bytes':memory,'cpu':cpu}


def parse_profiles(directory, *, compressed=False):
    pattern = 'callgrind.out*.gz' if compressed else 'callgrind.out*'
    paths = sorted(Path(directory).glob(pattern))
    require(bool(paths),'missing instruction profiles')
    phases = {}
    terminated = 0
    parts = set()
    for path in paths:
        if not compressed and path.suffix == '.gz':
            continue
        text = gzip.decompress(path.read_bytes()).decode() if compressed else path.read_text()
        def one(prefix):
            values = [l[len(prefix):].strip() for l in text.splitlines() if l.startswith(prefix)]
            require(len(values)==1, 'missing/duplicate profile field '+prefix)
            return values[0]
        require(one('events:')=='Ir','incompatible instruction unit')
        part = one('part:')
        require(part not in parts,'duplicate profile interval')
        parts.add(part)
        total = int(one('summary:'))
        require(total >= 0 and int(one('totals:'))==total,'invalid profile total')
        trigger = one('desc: Trigger:')
        if trigger == 'Program termination':
            phase = 'reporting_and_cleanup'
            terminated += 1
        else:
            require(trigger.startswith('Client Request: '),'unknown profiling boundary')
            phase = trigger.removeprefix('Client Request: ')
        require(phase in PHASES,'unknown cost phase')
        phases[phase] = phases.get(phase,0)+total
    require(terminated==1,'missing or duplicate final interval')
    stderr = (Path(directory)/'stderr.txt').read_text()
    collected = re.findall(r'Collected\s*:\s*(\d+)',stderr)
    require(len(collected)==1,'missing whole-process instruction checksum')
    require(sum(phases.values())==int(collected[0]),'phase costs do not sum to process cost')
    return phases


def frozen_inputs(round_dir):
    c = read(round_dir/'contract.json')
    require(objhash(c)==read(round_dir/'seal.json')['contract_sha256'],'changed contract')
    for relative, expected in c['pinned_files'].items():
        require(digest(round_dir/relative)==expected,'changed pinned artifact '+relative)
    require(digest(Path(__file__))==c['evaluator_sha256']['tournament.py'],'run the frozen evaluator copy')
    require(digest(HERE/'oracle.py')==c['evaluator_sha256']['oracle.py'],'changed checker')
    return c, read(round_dir/'fixtures.json'), read(round_dir/'candidates.json')


def snapshot_build(source,destination):
    snap = destination/'source'
    snap.mkdir()
    for name in ('Cargo.toml','Cargo.lock','build.rs'):
        if (source/name).exists():
            shutil.copy2(source/name,snap/name)
    shutil.copytree(source/'src',snap/'src')
    # Cache identities include native sources even in CPU-only builds.
    # Follow literal include dependencies so the source seal covers them too.
    for rust in (source/'src').rglob('*.rs'):
        for relative in re.findall(r'include(?:_str|_bytes)?!\s*\(\s*"([^"]+)"',rust.read_text()):
            dependency=(rust.parent/relative).resolve()
            require(dependency.is_relative_to(source),'include escapes source root')
            include_target=snap/dependency.relative_to(source)
            include_target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(dependency,include_target)
    (snap/'examples').mkdir()
    shutil.copy2(source/'examples/ic_tournament_worker.rs',snap/'examples/ic_tournament_worker.rs')
    manifest = {str(p.relative_to(snap)):digest(p) for p in sorted(snap.rglob('*')) if p.is_file()}
    write(destination/'source-manifest.json',manifest,exclusive=True)
    with (destination/'build.log').open('w') as log:
        subprocess.run(['cargo','build','--release','--offline','--locked','--jobs','2',
            '--example','ic_tournament_worker','--target-dir',str(destination/'build')],
            cwd=snap,env=child_env(),stdout=log,stderr=subprocess.STDOUT,check=True)
    binary = destination/'worker'
    shutil.copy2(destination/'build/release/examples/ic_tournament_worker',binary)
    return binary, manifest


def prepare(args):
    out = args.out.resolve()
    out.mkdir(parents=True,exist_ok=False)
    require(platform.machine()=='x86_64','instruction protocol currently supports amd64 only')
    require(shutil.which('valgrind') is not None,'Valgrind required')
    version = subprocess.check_output(['valgrind','--version'],text=True).strip()
    require(version=='valgrind-3.22.0','version requires a new calibrated protocol')
    evaluator = out/'evaluator'
    evaluator.mkdir()
    for name in ('tournament.py','oracle.py'):
        shutil.copy2(HERE/name,evaluator/name)
    source = args.source_root.resolve()
    binary, manifest = snapshot_build(source,out)
    arms = read(args.candidates) if args.candidates else candidates()
    require(isinstance(arms,list) and 2 <= len(arms) <= 16,'candidate count must be 2..16')
    require(arms[0]['id']=='incumbent','first candidate must be incumbent')
    ids, hashes = set(),set()
    for arm in arms:
        name = arm['id']
        require(re.fullmatch(r'[a-z][a-z0-9_-]{0,31}',name) is not None and name not in ids,'invalid candidate id')
        require(name not in {'rho','aa_control'},'reserved candidate id')
        if arm.get('source_root'):
            require(name!='incumbent','use --source-root for the incumbent')
            destination=out/'source_candidates'/name
            destination.mkdir(parents=True)
            arm_binary, arm_manifest=snapshot_build(Path(arm['source_root']).resolve(),destination)
            arm['binary_relative']=str(arm_binary.relative_to(out))
            arm['source_directory']=str((destination/'source').relative_to(out))
            arm['source_manifest_relative']=str((destination/'source-manifest.json').relative_to(out))
            arm.pop('source_root')
        else:
            arm_manifest=manifest
            arm['binary_relative']='worker'
            arm['source_directory']='source'
            arm['source_manifest_relative']='source-manifest.json'
        h = objhash(arm['config'])
        identity=objhash([objhash(arm_manifest),h])
        require(identity not in hashes,'duplicate candidate source/configuration')
        ids.add(name); hashes.add(identity)
        arm['configuration_sha256'] = h
        arm['source_manifest_sha256'] = objhash(arm_manifest)
    write(out/'candidates.json',arms,exclusive=True)
    cpus = sorted(os.sched_getaffinity(0))
    cpu = args.cpu if args.cpu is not None else cpus[-1]
    require(cpu in cpus,'CPU outside permitted affinity')
    # A pilot deliberately limits sample count; every confirmation has 60 new
    # distinct public-target jobs over five cells, three repetitions per arm.
    cells = [(13,0),(17,1),(19,0),(23,0)]
    profile = {'development':3,'selection':3,'confirmation':12}
    if args.profile=='standard':
        profile = {'development':30,'selection':30,'confirmation':100}
    limits = {'timeout_seconds':args.timeout,'memory_bytes':8*1024**3,'cpu':cpu,
              'worker_threads':1,'max_profiled_jobs':args.max_processes}
    fixtures = {}
    rng = random.Random(args.seed)
    for stage in STAGES:
        stage_cells = cells+[(19,1)] if stage in ('confirmation','replay') else cells
        count = profile.get(stage,1)
        if stage=='replay':
            fixtures[stage] = copy.deepcopy(fixtures['confirmation'])
            continue
        cases = []
        for degree,a in stage_cells:
            for index in range(count):
                public_seed = rng.getrandbits(64)
                case = {'id':f'n{degree}a{a}-{index:03d}', 'cell':f'n{degree}a{a}',
                    'job':{'mode':'fixture','degree':degree,'curve_a':a,
                           'target_seeds':[public_seed],'algorithm_seed':rng.getrandbits(64),
                           'factor_base':{'kind':'subgroup_orbits','seed':43,'points':6*degree},
                           'config':BASE_CONFIG}}
                raw = subprocess.run([str(binary)],input=json.dumps(case['job']),text=True,
                    capture_output=True,env=child_env(),timeout=args.timeout,check=True)
                case['fixture'] = json.loads(raw.stdout)['fixture']
                Curve(case['fixture'])
                case['fixture_sha256'] = objhash(case['fixture'])
                cases.append(case)
        fixtures[stage] = cases
    write(out/'fixtures.json',fixtures,exclusive=True)
    write(out/'calibration.json',{
        'unit':UNIT,'conversion':'1 Ir event = 1 user-space guest instruction; identity, no curve-op conversion',
        'all_intervals_sum_to_collected_required':True,
        'boundary':'complete worker user-space execution, including startup, allocations, reporting and all algorithm phases',
        'excluded':'kernel/device work, profiler implementation and external independent audit; reported separately, not silently converted',
        'restrictions':'CPU only; fixed compiler/ISA/profiler. This is implementation cost, not a hardware-independent arithmetic complexity claim.',
        'floor':'For this required full-rank collector, K columns require at least K relation-producing trials and at least K guest instructions. Weak implementation-specific floor; never evidence of a non-generic advance.',
        'paired_aa_gate':'same executable/config; both complete, no >=30% promotion, each cell within 5%'},exclusive=True)
    c = {'schema_version':1,'profile':args.profile,'seed':args.seed,'created_unix':time.time(),
        'unit':UNIT,'evidence_scope':'bounded public-hash ECDLP configuration tournament',
        'metric_class':'implementation_instruction_cost','family_wide_or_scaling_claim':False,
        'equivalent_suite_reason':'Point-base collector/rank/descent API, not WDSat ANF/conflict protocol. Fresh reference/candidates, independent point and scalar-field certificates; final suite has 60 inputs in pilot.',
        'curve_diversity_limit':'One Koblitz curve per development degree; additional holdout curve in confirmation. Does not meet three curves per size for a broad family claim.',
        'limits':limits,'repetitions':3,'confirmation_ratio':0.7,'max_cell_ratio':1.1,
        'ci_level':0.95,'bootstrap_draws':2000,'host':platform.uname()._asdict(),
        'profiler_version':version,'compiler':subprocess.check_output(['rustc','--version'],text=True).strip(),
        'evaluator_sha256':{n:digest(evaluator/n) for n in ('tournament.py','oracle.py')},
        'source_manifest_sha256':objhash(manifest),
        'target_count':1,'native_timings':'paired native reruns retained; no runtime claim from profiled elapsed time',
        'pinned_files':{n:digest(out/n) for n in set(['worker','source-manifest.json','candidates.json','fixtures.json','calibration.json']+[a['binary_relative'] for a in arms]+[a['source_manifest_relative'] for a in arms])}}
    write(out/'contract.json',c,exclusive=True)
    write(out/'seal.json',{'contract_sha256':objhash(c)},exclusive=True)
    print(json.dumps({'status':'prepared','round':str(out),'command':f'python3 {evaluator / "tournament.py"} run --round {out}'}))


def trial_path(root, stage, case, arm, repetition):
    return root/'runs'/stage/case['id']/arm/f'rep-{repetition}'


def verify_receipt(directory, c, case, arm):
    r = read(directory/'receipt.json')
    require(r['case']==case['id'] and r['cell']==case['cell'],'changed case identity')
    require(r['arm']==arm['id'] and r['unit']==c['unit'],'changed arm or unit')
    require(r['repetition']==int(directory.name.removeprefix('rep-')),'changed repetition')
    require(r['case_sha256']==objhash(case),'changed case binding')
    require(r['arm_sha256']==objhash(arm),'changed candidate binding')
    for path, h in r['artifacts'].items():
        require(digest(directory/path)==h,'changed trial artifact '+path)
    if r['status']=='VERIFIED':
        report = read(directory/'profile/stdout.json')
        proof = verify(report,case['fixture'],expected_mode=r['mode'],summands=arm['config']['summands'])
        require(proof==r['certificate'],'certificate summary changed')
        costs = parse_profiles(directory/'profile',compressed=True)
        require(costs==r['phase_costs'] and sum(costs.values())==r['total_operations'],'changed cost summary')
        native = read(directory/'native/stdout.json')
        native_proof = verify(native,case['fixture'],expected_mode=r['mode'],summands=arm['config']['summands'])
        require(native_proof['solutions']==proof['solutions'],'native/profile mismatch')
    return r


def run_trial(root, c, stage, case, arm, repetition):
    directory = trial_path(root,stage,case,arm['id'],repetition)
    if (directory/'receipt.json').exists():
        return verify_receipt(directory,c,case,arm)
    require(not directory.exists(),'interrupted trial retained; start a new campaign rather than overwrite it')
    directory.mkdir(parents=True)
    job = copy.deepcopy(case['job'])
    job.update(mode='rho' if arm['id']=='rho' else 'ic',config=arm['config'])
    write(directory/'job.json',job,exclusive=True)
    binary=root/arm.get('binary_relative','worker')
    command = ['valgrind','--tool=callgrind','--cache-sim=no','--branch-sim=no',
        '--separate-threads=no','--collect-atstart=yes','--instr-atstart=yes',
        '--callgrind-out-file='+str(directory/'profile/callgrind.out'),str(binary)]
    lim = c['limits']
    run = execute(command,job,directory/'profile',lim['timeout_seconds'],lim['memory_bytes'],lim['cpu'])
    r = {'schema_version':1,'stage':stage,'arm':arm['id'],'case':case['id'],'cell':case['cell'],
         'case_sha256':objhash(case),'arm_sha256':objhash(arm),'repetition':repetition,
         'unit':c['unit'],'mode':job['mode'],'status':'ERROR','profile_process':run,
         'total_operations':None,'phase_costs':None,'certificate':None}
    try:
        require(run['process_status']!='TIMEOUT','TIMEOUT')
        require(run['exit_code']==0,'worker failed or incomplete')
        report = read(directory/'profile/stdout.json')
        proof = verify(report,case['fixture'],expected_mode=job['mode'],summands=job['config']['summands'])
        costs = parse_profiles(directory/'profile')
        expected = {'startup_and_input','curve_and_targets','final_verification','reporting_and_cleanup'}
        expected |= {'rho_solve'} if job['mode']=='rho' else PHASES-{'rho_solve'}
        require(set(costs)==expected,'missing exclusive phases')
        native = execute([str(binary)],job,directory/'native',lim['timeout_seconds'],lim['memory_bytes'],lim['cpu'])
        r['native_process'] = native
        require(native['exit_code']==0 and native['process_status']=='EXITED','native worker failed')
        nproof = verify(read(directory/'native/stdout.json'),case['fixture'],expected_mode=job['mode'],summands=job['config']['summands'])
        require(nproof['solutions']==proof['solutions'],'native/profile answers differ')
        total = sum(costs.values())
        r.update(status='VERIFIED',total_operations=total,phase_costs=costs,certificate=proof,
                 normalized_S=total/math.sqrt(int(case['fixture']['subgroup_order'])),
                 floor_operations=proof['rank'],
                 ratio_to_floor=total/proof['rank'] if proof['rank'] else None)
    except (InvalidEvidence,ValueError,KeyError,TypeError) as exc:
        r['reason'] = str(exc)
        if run['process_status']=='TIMEOUT':
            r['status']='TIMEOUT'
        elif 'memory allocation' in (directory/'profile/stderr.txt').read_text():
            r['status']='OOM'
        else:
            r['status']='INVALID_OR_INCOMPLETE'
    for path in sorted((directory/'profile').glob('callgrind.out*')):
        zipped = path.with_name(path.name+'.gz')
        zipped.write_bytes(gzip.compress(path.read_bytes(),mtime=0))
        path.unlink()
    r['artifacts'] = {str(p.relative_to(directory)):digest(p) for p in sorted(directory.rglob('*')) if p.is_file()}
    write(directory/'receipt.json',r,exclusive=True)
    return r


def synthetic_arm(name, base):
    arm = copy.deepcopy(base)
    arm['id']=name
    return arm


def load_stage(root, stage, fixtures, arms, repetitions):
    rows=[]
    for case in fixtures:
        for arm in arms:
            for rep in range(repetitions):
                path=trial_path(root,stage,case,arm['id'],rep)/'receipt.json'
                require(path.exists(),'missing scheduled trial')
                rows.append(read(path))
    return rows


def comparison(rows, candidate_id, *, baseline='incumbent', draws=2000):
    selected=[r for r in rows if r['arm'] in (candidate_id,baseline)]
    cases={r['case'] for r in rows}
    grouped={}
    for r in selected:
        grouped.setdefault((r['case'],r['arm']),[]).append(r)
    logs={}
    native_logs={}
    for case in sorted(cases):
        a=grouped.get((case,baseline),[]); b=grouped.get((case,candidate_id),[])
        if not a or not b or len(a)!=len(b) or len({x['repetition'] for x in a})!=len(a) or len({x['repetition'] for x in b})!=len(b) or {x['repetition'] for x in a}!={x['repetition'] for x in b}:
            return {'candidate':candidate_id,'eligible':False,'reason':'missing paired trials'}
        if any(x['status']!='VERIFIED' or x['total_operations'] is None for x in a+b):
            return {'candidate':candidate_id,'eligible':False,'reason':'unverified or unpriced workload'}
        require(len({x['case_sha256'] for x in a+b})==1,'changed paired fixtures')
        if candidate_id!='rho' and baseline!='rho':
            require(len({x['certificate']['factor_base_sha256'] for x in a+b})==1,'changed paired factor-base support')
        cell=a[0]['cell']
        ratio=statistics.median(x['total_operations'] for x in b)/statistics.median(x['total_operations'] for x in a)
        logs.setdefault(cell,[]).append(math.log(ratio))
        wall_ratio=statistics.median(x['native_process']['process_wall_seconds'] for x in b)/statistics.median(x['native_process']['process_wall_seconds'] for x in a)
        native_logs.setdefault(cell,[]).append(math.log(wall_ratio))
    require(bool(logs),'empty comparison')
    cell_values={cell:statistics.mean(xs) for cell,xs in logs.items()}
    estimate=math.exp(statistics.mean(cell_values.values()))
    rng=random.Random(751203)
    cells=sorted(logs)
    boot=[]
    for _ in range(draws):
        means=[]
        for _ in cells:
            cell=rng.choice(cells)
            means.append(statistics.mean(rng.choices(logs[cell],k=len(logs[cell]))))
        boot.append(math.exp(statistics.mean(means)))
    boot.sort()
    return {'candidate':candidate_id,'eligible':True,'candidate_over_baseline':estimate,
        'speedup':1/estimate,'ci95':[boot[int(.025*draws)],boot[min(draws-1,int(.975*draws))]],
        'per_cell':{cell:math.exp(v) for cell,v in cell_values.items()},
        'native_wall_candidate_over_baseline':math.exp(statistics.mean(statistics.mean(v) for v in native_logs.values())),
        'native_wall_status':'diagnostic; separate timing-confidence gate not applied',
        'paired_cases':len(cases),'independent_curve_blocks':len(cells)}


def gate(result,c):
    return bool(result.get('eligible') and result['candidate_over_baseline']<=c['confirmation_ratio']
        and result['ci95'][1]<1 and max(result['per_cell'].values())<=c['max_cell_ratio'])


def stage_arms(root, stage, arms):
    base=arms[0]
    if stage=='aa':
        return [base,synthetic_arm('aa_control',base)]
    if stage=='smoke':
        return arms+[synthetic_arm('rho',base)]
    if stage=='development':
        smoke=read(root/'summaries/smoke.json')
        failed={r['arm'] for r in smoke['failures']}
        require('incumbent' not in failed,'incumbent failed smoke')
        return [a for a in arms if a['id'] not in failed]+[synthetic_arm('rho',base)]
    if stage=='selection':
        d=read(root/'summaries/development.json')
        eligible=sorted([r for r in d['comparisons'] if r.get('eligible')],key=lambda r:r['candidate_over_baseline'])[:2]
        return [base]+[next(a for a in arms if a['id']==r['candidate']) for r in eligible]+[synthetic_arm('rho',base)]
    d=read(root/'summaries/selection.json')
    provisional=d.get('provisional_challenger')
    return [base]+([next(a for a in arms if a['id']==provisional)] if provisional else [])+[synthetic_arm('rho',base)]


def summarize(root,c,stage,fixtures,arms,*,save=True):
    rows=load_stage(root,stage,fixtures,arms,c['repetitions'])
    comps=[comparison(rows,a['id'],draws=c['bootstrap_draws']) for a in arms if a['id'] not in ('incumbent','rho')]
    rho = comparison(rows,'rho',draws=c['bootstrap_draws']) if any(a['id']=='rho' for a in arms) else None
    result={'stage':stage,'runs':len(rows),'verified_runs':sum(r['status']=='VERIFIED' for r in rows),
            'comparisons':comps,'rho_over_incumbent':rho,
            'process_wall_seconds_including_profiling':sum(r['profile_process']['process_wall_seconds']+r.get('native_process',{}).get('process_wall_seconds',0) for r in rows),
            'failures':[{k:r.get(k) for k in ('case','arm','repetition','status','reason')} for r in rows if r['status']!='VERIFIED']}
    if stage=='aa':
        aa=comps[0]
        result['passed']=bool(aa.get('eligible') and all(.95<=r<=1.05 for r in aa['per_cell'].values()) and not gate(aa,c))
    if stage=='selection':
        eligible=sorted([r for r in comps if r.get('eligible')],key=lambda r:r['candidate_over_baseline'])
        result['provisional_challenger']=eligible[0]['candidate'] if eligible else None
    if save:
        write(root/'summaries'/f'{stage}.json',result,exclusive=True)
    return result


def decision(root,c,fixtures,all_arms,*,save=True):
    select=read(root/'summaries/selection.json')
    confirm=read(root/'summaries/confirmation.json')
    replay=read(root/'summaries/replay.json')
    challenger=select['provisional_challenger']
    conf=next((r for r in confirm['comparisons'] if r['candidate']==challenger),{})
    rep=next((r for r in replay['comparisons'] if r['candidate']==challenger),{})
    passed=gate(conf,c) and gate(rep,c)
    qualified=bool(conf.get('eligible') and rep.get('eligible'))
    base_complete=True
    for stage in ('confirmation','replay'):
        for case in fixtures[stage]:
            for repetition in range(c['repetitions']):
                receipt=read(trial_path(root,stage,case,'incumbent',repetition)/'receipt.json')
                base_complete &= receipt['status']=='VERIFIED' and receipt['total_operations'] is not None
    choice=challenger if passed else ('incumbent' if base_complete and (qualified or not challenger) else None)
    reasons=[]
    if not passed:
        reasons.append('No challenger passed every confirmation and replay threshold.')
    result={'status':'promoted' if passed else ('retained' if choice else 'inconclusive'),
        'winner':choice,'provisional_challenger':challenger,'confirmation':conf,'replay':rep,
        'unit':c['unit'],'classification':'engineering' if passed else 'accounting',
        'scope':c['evidence_scope'],'family_wide_or_complexity_claim':False,
        'arithmetic_operation_speedup':None,'instruction_speedup':conf.get('speedup') if passed else None,
        'reasons':reasons,'independent_replay':'fresh processes, frozen executable, same final cases; independent Python group and rank checker',
        'rho_comparison':confirm['rho_over_incumbent'],
        'beats_rho':False}
    # Rho's cost ratio to the incumbent divided by the winner's ratio gives rho/winner.
    rc=confirm['rho_over_incumbent']
    if choice and rc and rc.get('eligible'):
        ratio = conf['candidate_over_baseline'] if passed else 1
        result['rho_over_winner']=rc['candidate_over_baseline']/ratio
        result['beats_rho']=result['rho_over_winner']>1
        result['rho_verdict_scope']='observed complete instruction cost; no asymptotic crossover claim'
    if save:
        write(root/'decision.json',result,exclusive=True)
    return result


def run_campaign(args):
    root=args.round.resolve()
    with (root/'operation.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        c,fixtures,arms=frozen_inputs(root)
        completed=sum(1 for _ in (root/'runs').glob('**/receipt.json'))
        for stage in STAGES:
            if args.stage!='all' and stage!=args.stage:
                continue
            index=STAGES.index(stage)
            for prior in STAGES[:index]:
                require((root/'summaries'/f'{prior}.json').exists(),'complete earlier stage '+prior)
            if (root/'summaries'/f'{stage}.json').exists():
                if stage=='aa':
                    require(read(root/'summaries/aa.json')['passed'],'A/A gate failed')
                continue
            active=stage_arms(root,stage,arms)
            schedule=[(case,rep) for case in fixtures[stage] for rep in range(c['repetitions'])]
            rng=random.Random(c['seed']+index)
            rng.shuffle(schedule)
            for case,rep in schedule:
                order=list(active);rng.shuffle(order)
                for arm in order:
                    present=(trial_path(root,stage,case,arm['id'],rep)/'receipt.json').exists()
                    require(present or completed<c['limits']['max_profiled_jobs'],'campaign process budget exhausted')
                    r=run_trial(root,c,stage,case,arm,rep)
                    completed+=not present
                    write(root/'state.json',{'stage':stage,'completed_jobs':completed,'last_case':case['id'],
                          'last_arm':arm['id'],'last_status':r['status'],'updated_unix':time.time()})
                    print(json.dumps({'stage':stage,'arm':arm['id'],'case':case['id'],'rep':rep,
                                      'status':r['status'],'operations':r['total_operations']}),flush=True)
            summary=summarize(root,c,stage,fixtures[stage],active)
            if stage=='aa':
                require(summary['passed'],'A/A gate failed; preserve evidence and investigate')
        if (root/'summaries/replay.json').exists() and not (root/'decision.json').exists():
            print(json.dumps(decision(root,c,fixtures,arms)),flush=True)


def audit(args):
    root=args.round.resolve()
    c,fixtures,arms=frozen_inputs(root)
    manifest=read(root/'source-manifest.json')
    source_pairs={(a.get('source_manifest_relative','source-manifest.json'),a.get('source_directory','source')) for a in arms}
    for manifest_path,source_path in source_pairs:
        for name,h in read(root/manifest_path).items():
            require(digest(root/source_path/name)==h,'changed source snapshot '+name)
    count=0
    for stage in STAGES:
        if not (root/'summaries'/f'{stage}.json').exists():
            continue
        active=stage_arms(root,stage,arms)
        for case in fixtures[stage]:
            for arm in active:
                for rep in range(c['repetitions']):
                    verify_receipt(trial_path(root,stage,case,arm['id'],rep),c,case,arm)
                    count+=1
        # Recompute selection, statistics, completion counts and final verdict.
        saved=read(root/'summaries'/f'{stage}.json')
        require(saved==summarize(root,c,stage,fixtures[stage],active,save=False),
                'changed stage summary or provisional selection')
    if (root/'decision.json').exists():
        require(read(root/'decision.json')==decision(root,c,fixtures,arms,save=False),'changed decision')
    print(json.dumps({'status':'VERIFIED','trial_receipts':count,'source_files':len(manifest)}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('prepare')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source-root',type=Path,default=ROOT)
    p.add_argument('--candidates',type=Path)
    p.add_argument('--profile',choices=['pilot','standard'],default='pilot')
    p.add_argument('--seed',type=int,default=20260915)
    p.add_argument('--cpu',type=int)
    p.add_argument('--timeout',type=float,default=60)
    p.add_argument('--max-processes',type=int,default=1800)
    p=commands.add_parser('propose')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--from-round',type=Path)
    for name in ('run','verify','status'):
        p=commands.add_parser(name);p.add_argument('--round',type=Path,required=True)
        if name=='run':p.add_argument('--stage',choices=['all']+STAGES,default='all')
    args=parser.parse_args()
    try:
        if args.command=='prepare':
            require(args.timeout>0 and args.max_processes>0,'positive limits required')
            prepare(args)
        elif args.command=='propose':
            if args.from_round:
                proposed,source=proposals_from_previous(args.from_round.resolve())
                write(args.out,proposed,exclusive=True)
                print(json.dumps({'candidates':str(args.out),'baseline_source_root':str(source),
                      'next_step':'Use this --source-root and a new campaign seed when preparing the next round.'}))
            else:
                write(args.out,candidates(),exclusive=True)
                print(args.out)
        elif args.command=='run':run_campaign(args)
        elif args.command=='verify':audit(args)
        else:
            root=args.round.resolve()
            print(json.dumps(read(root/'decision.json') if (root/'decision.json').exists() else
                read(root/'state.json') if (root/'state.json').exists() else {'status':'prepared'},indent=2))
    except (InvalidEvidence,ValueError,OSError,subprocess.SubprocessError) as exc:
        print(f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 1
    return 0


if __name__=='__main__':
    raise SystemExit(main())
