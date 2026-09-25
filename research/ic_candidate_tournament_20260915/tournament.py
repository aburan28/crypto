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
import threading

from oracle import Curve, InvalidEvidence, require, verify
from portfolio import retain
from driver_admission import (EVALUATOR, check_admission, freeze_admission, producer_metadata, run_record, online_table, distinct_candidates)

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


SPAWN = 'posix_spawn child without a fork of the evaluator; caps by prlimit before job delivery'


def execute(command, job, directory, timeout, memory, cpu):
    # Through round 0009 the child was created with a preexec_fn, which makes
    # CPython fork() the evaluator: the wall then charged a copy of the evaluator's
    # page tables to every job, a cost that grows with the evaluator's heap and
    # dominated the smallest cells' native times with a noise unrelated to either
    # arm. Without preexec_fn CPython uses vfork/posix_spawn. The child inherits
    # the calling thread's affinity, so it is pinned before the spawn; the memory
    # and core caps are applied with prlimit while the child is still blocked on
    # stdin, i.e. before it can allocate anything for the job; the watchdog thread
    # is started before the timing window opens.
    directory.mkdir(parents=True,exist_ok=True)
    payload = json.dumps(job,sort_keys=True).encode()
    lock = threading.Lock()
    holder = {'process':None,'fired':False}
    expired = threading.Event()
    def kill(process):
        try:
            os.killpg(process.pid,signal.SIGKILL)
            expired.set()
        except ProcessLookupError:
            pass
    def watchdog():
        with lock:
            holder['fired'] = True
            process = holder['process']
        if process is not None:
            kill(process)
    timer = threading.Timer(timeout,watchdog)
    timer.daemon = True
    status = 'EXITED'
    with (directory/'stdout.json').open('w') as stdout, (directory/'stderr.txt').open('w') as stderr:
        timer.start()
        inherited = os.sched_getaffinity(0)
        os.sched_setaffinity(0,{cpu})
        start_ns = time.monotonic_ns()
        try:
            process = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=stdout,stderr=stderr,
                env=child_env(),start_new_session=True)
        finally:
            os.sched_setaffinity(0,inherited)
        with lock:
            holder['process'] = process
            fired = holder['fired']
        if fired:
            kill(process)
        resource.prlimit(process.pid,resource.RLIMIT_AS,(memory,memory))
        resource.prlimit(process.pid,resource.RLIMIT_CORE,(0,0))
        try:
            # Blocking reap avoids communicate(timeout)'s exponential wait polling,
            # which quantized previous short native timings. Charge full process wall.
            # stdout/stderr already go to files, so no pipe-draining loop is
            # needed. Reap this exact child with wait4 to retain its RSS peak;
            # RUSAGE_CHILDREN would mix peaks from earlier jobs.
            try:
                process.stdin.write(payload)
            except BrokenPipeError:
                pass
            finally:
                try:
                    process.stdin.close()
                except BrokenPipeError:
                    pass
            _, wait_status, usage = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(wait_status)
        finally:
            wall_ns = time.monotonic_ns()-start_ns
            timer.cancel()
            timer.join()
        if expired.is_set(): status = 'TIMEOUT'
    return {'exit_code':process.returncode,'process_wall_seconds':wall_ns/1_000_000_000,
            'process_wall_ns':wall_ns, 'peak_rss_bytes':usage.ru_maxrss*1024,
            'process_status':status,'command':command,'memory_cap_bytes':memory,'cpu':cpu,'spawn':SPAWN}


def parse_profiles(directory, *, compressed=False, phase_schema=1):
    require(phase_schema in (1, 2, 3), 'unknown phase schema')
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
            phase = 'reporting_and_cleanup' if phase_schema == 1 else 'setup'
            terminated += 1
        else:
            require(trigger.startswith('Client Request: '),'unknown profiling boundary')
            phase = trigger.removeprefix('Client Request: ')
            if phase_schema in (2, 3):
                require(phase.startswith('ic_') or phase == 'reference_solve',
                        'legacy interval in scientific profile')
                phase = phase.removeprefix('ic_')
        if phase_schema == 1:
            require(phase in PHASES,'unknown cost phase')
        else:
            allowed = {'setup', 'factor_base', 'precompute', 'queries', 'pdp',
                'relation_check', 'matrix_build', 'relation_la', 'target_descent',
                'recovery_check', 'reference_solve'}
            if phase_schema == 3:
                allowed |= {'target_query', 'target_pdp', 'target_relation_check'}
            require(phase in allowed, 'unknown scientific cost phase')
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
    for name, expected in c['evaluator_sha256'].items():
        require(digest(HERE/name)==expected,'changed evaluator or checker: '+name)
    return c, read(round_dir/'fixtures.json'), read(round_dir/'candidates.json')


def built_worker(build_dir):
    # `[build] target` in a snapshot's cargo config moves the artifact under the
    # target triple; either way exactly one worker must have been produced.
    found = [p for p in [build_dir/'release/examples/ic_tournament_worker']
             +sorted(build_dir.glob('*/release/examples/ic_tournament_worker')) if p.is_file()]
    require(len(found)==1,'expected exactly one built worker, found %d'%len(found))
    return found[0]


def snapshot_build(source,destination, *, scientific=False):
    snap = destination/'source'
    snap.mkdir()
    # A cargo config inside the source root is part of the snapshot and its seal:
    # it is how a source tree selects link mode or target features, which rustc
    # cannot take from Cargo.toml.
    for name in ('Cargo.toml','Cargo.lock','build.rs','.cargo/config.toml'):
        if (source/name).exists():
            (snap/name).parent.mkdir(parents=True,exist_ok=True)
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
    if scientific:
        metadata = producer_metadata(source, manifest, subprocess.check_output(['rustc', '--version'], text=True).strip())
        write(destination/'producer.json', metadata, exclusive=True)
        write(destination/'preparation.json', metadata['preparation'], exclusive=True)
    with (destination/'build.log').open('w') as log:
        subprocess.run(['cargo','build','--release','--offline','--locked','--jobs','2',
            '--example','ic_tournament_worker','--target-dir',str(destination/'build')],
            cwd=snap,env=dict(child_env(), IC_SOURCE_MANIFEST_SHA256=objhash(manifest)),stdout=log,stderr=subprocess.STDOUT,check=True)
    binary = destination/'worker'
    shutil.copy2(built_worker(destination/'build'),binary)
    return binary, manifest


def prepare(args):
    out = args.out.resolve()
    require(args.targets == 1, 'scientific admission requires a single public target')
    require(args.candidates is not None, 'supply an explicit registry of admitted optimized candidates')
    require(platform.machine()=='x86_64','instruction protocol currently supports amd64 only')
    require(shutil.which('valgrind') is not None,'Valgrind required')
    version = subprocess.check_output(['valgrind','--version'],text=True).strip()
    require(version=='valgrind-3.22.0','version requires a new calibrated protocol')
    out.mkdir(parents=True,exist_ok=False)
    evaluator = out/'evaluator'
    evaluator.mkdir()
    for name in EVALUATOR:
        (evaluator/name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(HERE/name,evaluator/name)
    source = args.source_root.resolve()
    binary, manifest = snapshot_build(source,out,scientific=True)
    arms = read(args.candidates) if args.candidates else candidates()
    require(isinstance(arms,list) and 2 <= len(arms) <= 16,'candidate count must be 2..16')
    require(arms[0]['id']=='incumbent','first candidate must be incumbent')
    ids, hashes = set(),set()
    built_sources = {source:(out,binary,manifest)}
    for arm in arms:
        name = arm['id']
        require(re.fullmatch(r'[a-z][a-z0-9_-]{0,31}',name) is not None and name not in ids,'invalid candidate id')
        require(name not in {'rho','aa_control'},'reserved candidate id')
        if arm.get('source_root'):
            require(name!='incumbent','use --source-root for the incumbent')
            candidate_source=Path(arm['source_root']).resolve()
            if candidate_source in built_sources:
                destination,arm_binary,arm_manifest=built_sources[candidate_source]
            else:
                destination=out/'source_candidates'/name
                destination.mkdir(parents=True)
                arm_binary,arm_manifest=snapshot_build(candidate_source,destination,scientific=True)
                built_sources[candidate_source]=(destination,arm_binary,arm_manifest)
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
    # A stronger reference may live outside the IC incumbent's source tree.
    # Freeze it explicitly rather than silently weakening rho when changing IC.
    rho_reference = synthetic_arm('rho',arms[0])
    if args.rho_source_root:
        destination=out/'rho_reference'
        destination.mkdir()
        rho_binary,rho_manifest=snapshot_build(args.rho_source_root.resolve(),destination,scientific=True)
        rho_reference.update(binary_relative=str(rho_binary.relative_to(out)),
            source_directory=str((destination/'source').relative_to(out)),
            source_manifest_relative=str((destination/'source-manifest.json').relative_to(out)),
            source_manifest_sha256=objhash(rho_manifest))
    if args.rho_config:
        rho_reference['config']=read(args.rho_config)
        require(isinstance(rho_reference['config'],dict) and 'summands' in rho_reference['config'],
                'rho config must be a complete worker configuration')
    rho_reference['configuration_sha256']=objhash(rho_reference['config'])
    cpus = sorted(os.sched_getaffinity(0))
    cpu = args.cpu if args.cpu is not None else cpus[-1]
    require(cpu in cpus,'CPU outside permitted affinity')
    # A pilot deliberately limits sample count. The panel is declared, not
    # implied: --cells fixes the development and selection curve cells and
    # --holdout-cells the extra cells that appear only in confirmation and
    # replay. The defaults are the panel rounds 0006-0015 used, so a command
    # that does not name them reproduces those rounds' panel exactly. Widening
    # the panel is the only lever that shrinks the between-cell term of the
    # nested bootstrap, which the round-0016 pre-registration measures.
    cells = parse_cells(args.cells)
    holdout = parse_cells(args.holdout_cells)
    require(not (set(cells) & set(holdout)),'a holdout cell is already a development cell')
    profile = {'development':3,'selection':3,'confirmation':12}
    if args.profile=='standard':
        profile = {'development':30,'selection':30,'confirmation':100}
    # Per-cell confirmation allocation (additive; the default is the flat
    # profile above and reproduces every round before 0019).
    #
    # The strict rho gate asks every cell's ratio to be below one, and that
    # per-cell test is a point estimate of a geometric mean over the cell's
    # fixtures.  Its resolution is therefore the cell's own spread over its own
    # sample count, and rounds 0017 and 0018b measured that spread to differ by
    # a factor of nine across this panel: the winner/rho instruction ratio has
    # a per-case log spread of 0.037 at n13a0 and 0.329 at n23a1, because rho's
    # solve phase is where the variance lives and it dominates rho's cost only
    # at the large-subgroup cells.  A flat twelve therefore buys a 1e-8 failure
    # probability at one end of the panel and a 9% one at the other, and round
    # 0018b spent that 9%.
    #
    # `--confirmation-cases n23a1=48` raises a named cell's count.  Raises
    # only: a count below the flat profile is refused, so no round can weaken
    # a cell's evidence, and the allocation moves an estimate toward the truth
    # in whichever direction the truth lies -- it cannot buy a pass, only
    # resolve one.  The estimator and the gate are untouched; the cross-cell
    # mean is unweighted and stays unbiased under unequal counts.
    extra = parse_confirmation_allocation(args.confirmation_cases,profile['confirmation'],
                                          [f'n{n}a{a}' for n,a in cells+holdout])
    limits = {'timeout_seconds':args.timeout,'memory_bytes':8*1024**3,'cpu':cpu,
              'worker_threads':1,'max_profiled_jobs':args.max_processes}
    fixtures = {}
    used_targets = {}
    rng = random.Random(args.seed)
    for stage in STAGES:
        stage_cells = cells+holdout if stage in ('confirmation','replay') else cells
        count = profile.get(stage,1)
        if stage=='replay':
            fixtures[stage] = copy.deepcopy(fixtures['confirmation'])
            continue
        cases = []
        for degree,a in stage_cells:
            n_cases = extra.get(f'n{degree}a{a}',count) if stage=='confirmation' else count
            for index in range(n_cases):
                case = {'id':f'n{degree}a{a}-{index:03d}', 'cell':f'n{degree}a{a}',
                    'job':{'mode':'fixture','degree':degree,'curve_a':a,
                           'target_seeds':[],'algorithm_seed':rng.getrandbits(64),
                           'factor_base':{'kind':'subgroup_orbits','seed':43,'points':6*degree},
                           'config':arms[0]['config']}}
                used = used_targets.setdefault((degree,a),set())
                for attempt in range(1000):
                    case['job']['target_seeds']=[rng.getrandbits(64) for _ in range(args.targets)]
                    directory=out/'fixture_generation'/stage/case['id']/f'attempt-{attempt}'
                    process=execute([str(binary)],case['job'],directory,args.timeout,limits['memory_bytes'],cpu)
                    write(directory/'job.json',case['job'],exclusive=True)
                    write(directory/'process.json',process,exclusive=True)
                    require(process['exit_code']==0 and process['process_status']=='EXITED',
                            'fixture generation failed; raw preparation evidence retained')
                    fixture = read(directory/'stdout.json')['fixture']
                    curve = Curve(fixture)
                    require(curve.r-1-len(used)>=args.targets,'too few unused public targets for independent confirmation')
                    if reserve_targets(fixture,used):
                        case['fixture']=fixture
                        break
                else:
                    raise InvalidEvidence('could not sample distinct public targets within preparation budget')
                case['job']['public_targets'] = case['fixture']['targets']
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
        'paired_aa_gate':'same executable/config; both complete, no >=20% promotion, each cell within 5%'},exclusive=True)
    allocation = {f'n{n}a{a}':extra.get(f'n{n}a{a}',profile['confirmation']) for n,a in cells+holdout}
    confirmation_cases = sum(allocation.values())
    resources = {k: (str(v) if k=='timeout_seconds' else v) for k,v in limits.items() if k!='max_profiled_jobs'}
    all_admission_arms = arms+[synthetic_arm('aa_control', arms[0]), rho_reference]
    admissions = {}
    for stage, stage_cases in fixtures.items():
        for case in stage_cases:
            admitted_arms=[]
            # Admission is independent of selection; predeclare every possible arm.
            for arm in all_admission_arms:
                if stage=='aa' and arm['id'] not in ('incumbent', 'aa_control'):
                    continue
                if stage!='aa' and arm['id']=='aa_control':
                    continue
                relative = str(Path('admissions')/stage/case['id']/arm['id'])
                job = dict(copy.deepcopy(case['job']), mode='rho' if arm['id']=='rho' else 'ic', config=arm['config'])
                destination = (out/arm['source_manifest_relative']).parent
                admitted = freeze_admission(out/relative, binary=out/arm['binary_relative'], job=job,
                    fixture=case['fixture'], manifest=read(out/arm['source_manifest_relative']),
                    metadata=read(destination/'producer.json'), resources=resources,
                    worker_sha256=digest(out/arm['binary_relative']), execute=lambda binary, job, directory:
                        execute([str(binary)],job,directory,limits['timeout_seconds'],limits['memory_bytes'],limits['cpu']))
                admissions[relative] = objhash(admitted)
                admitted_arms.append((arm['id'],admitted))
            distinct_candidates(admitted_arms)
    c = {'schema_version':2, 'scientific_admission':True, 'reference_qualification':None, 'admissions':admissions, 'resources':resources,
        'run_number_base':int.from_bytes(os.urandom(16),'big') << 16,
        'run_aliases':[a['id'] for a in all_admission_arms],'profile':args.profile,'seed':args.seed,'created_unix':time.time(),
        'cells':[f'n{n}a{a}' for n,a in cells],'holdout_cells':[f'n{n}a{a}' for n,a in holdout],
        'confirmation_cases':confirmation_cases,
        'confirmation_cases_per_cell':allocation,
        'confirmation_allocation':('flat' if not extra else
            'per-cell; raised at the cells whose measured spread a flat count cannot resolve, '
            'never lowered, from frozen prior rounds only; estimator and gate unchanged'),
        'unit':UNIT,'evidence_scope':'bounded public-hash ECDLP configuration tournament',
        'metric_class':'implementation_instruction_cost','family_wide_or_scaling_claim':False,
        'equivalent_suite_reason':'Point-base collector/rank/descent API, not WDSat ANF/conflict protocol. Fresh reference/candidates, independent point and scalar-field certificates; the final suite has '+str(confirmation_cases)+' inputs.',
        'curve_diversity_limit':'One Koblitz curve per development degree; additional holdout curve in confirmation. Does not meet three curves per size for a broad family claim.',
        'limits':limits,'repetitions':3,'confirmation_ratio':0.8,'max_cell_ratio':1.1,
        'selection_width':args.selection_width,'exploration_slots':args.exploration_slots,
        'rho_reference':rho_reference,
        'require_native_progress':args.require_native_progress,'parity_margin':1.10,
        'objective':args.objective,'no_regression_ratio':0.98,
        'comparison_kind':args.comparison_kind,
        'support_contract':('Stable support within each arm/case; base policy varies across arms; same public ECDLPs and per-arm floors.'
                            if args.comparison_kind=='factor-base-policy' else 'Identical factor-base support across implementation arms.'),
        'native_timing_protocol':'blocking process reap with independent watchdog; complete cold process wall; '+SPAWN,
        'ci_level':0.95,'bootstrap_draws':2000,'host':platform.uname()._asdict(),
        'profiler_version':version,'compiler':subprocess.check_output(['rustc','--version'],text=True).strip(),
        'evaluator_sha256':{n:digest(evaluator/n) for n in EVALUATOR},
        'source_manifest_sha256':objhash(manifest),
        'target_count':args.targets,'workload':'one supplied public target; primary native online interval after reusable preparation through scalar replay; supplementary complete cold process','native_timings':'paired native reruns retained; no runtime claim from profiled elapsed time',
        'target_uniqueness':'Distinct public points within each curve across A/A, smoke, development, selection and confirmation; replay intentionally repeats confirmation.',
        'pinned_files':{n:digest(out/n) for n in set(['worker','source-manifest.json','candidates.json','fixtures.json','calibration.json']+[a['binary_relative'] for a in arms+[rho_reference]]+[a['source_manifest_relative'] for a in arms+[rho_reference]])}}
    c['pinned_files'].update({str(p.relative_to(out)):digest(p) for name in ('admissions','fixture_generation')
        for p in (out/name).rglob('*') if p.is_file()})
    for arm in all_admission_arms:
        directory=(out/arm['source_manifest_relative']).parent
        for name in ('producer.json','preparation.json','build.log'):
            path=directory/name
            c['pinned_files'][str(path.relative_to(out))]=digest(path)
    write(out/'contract.json',c,exclusive=True)
    write(out/'seal.json',{'contract_sha256':objhash(c)},exclusive=True)
    print(json.dumps({'status':'prepared','round':str(out),'command':f'python3 {evaluator / "tournament.py"} run --round {out}'}))


def reserve_targets(fixture,used):
    """New random seeds alone do not guarantee new points in a small group."""
    points=[tuple(map(int,p)) for p in fixture['targets']]
    if len(set(points))!=len(points) or any(p in used for p in points):
        return False
    used.update(points)
    return True


def parse_confirmation_allocation(text,floor,declared):
    """'n23a1=48' -> {'n23a1': 48}, against the profile's flat `floor`.

    Raises only.  A count below the floor is refused rather than clamped: the
    allocation exists to resolve a cell a flat count cannot, and a round that
    could also *lower* a count could buy a pass by measuring less.  Unknown
    cell names are refused too, so a typo shows up as a failed prepare rather
    than as a silently flat panel."""
    out = {}
    for item in (text or '').split(','):
        if not item.strip():
            continue
        name,sep,value = item.partition('=')
        require(bool(sep),f'--confirmation-cases wants cell=count, got {item!r}')
        name = name.strip()
        require(name not in out,f'{name}: named twice in --confirmation-cases')
        n = int(value)
        require(n>=floor,f'{name}: {n} confirmation cases is below the profile floor {floor}')
        out[name] = n
    unknown = sorted(set(out)-set(declared))
    require(not unknown,f'--confirmation-cases names cells outside the panel: {unknown}')
    return out


def parse_cells(text):
    """'13a0,17a1' -> [(13,0),(17,1)]. Order is preserved and duplicates are
    refused: a repeated cell would enter the outer bootstrap twice and silently
    understate the between-cell variance it is there to measure."""
    out = []
    for token in [t.strip() for t in str(text).split(',') if t.strip()]:
        m = re.fullmatch(r'n?(\d{1,2})a([01])',token)
        require(m is not None,f'cell {token!r} is not <degree>a<curve_a>')
        pair = (int(m.group(1)),int(m.group(2)))
        require(pair not in out,f'cell {token!r} named twice')
        out.append(pair)
    require(bool(out),'empty cell list')
    return out


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
        costs = parse_profiles(directory/'profile',compressed=True,phase_schema=3 if c.get('scientific_admission') else 1)
        require(costs==r['phase_costs'] and sum(costs.values())==r['total_operations'],'changed cost summary')
        native = read(directory/'native/stdout.json')
        native_proof = verify(native,case['fixture'],expected_mode=r['mode'],summands=arm['config']['summands'])
        require(native_proof['solutions']==proof['solutions'],'native/profile mismatch')
    if c.get('scientific_admission'):
        require(r['status'] in ('VERIFIED','TIMEOUT','OOM','INVALID_OR_INCOMPLETE'), 'unknown trial outcome')
        if r['status']=='VERIFIED':
            for key in ('profile_process','native_process'):
                process=r[key]
                require(process['process_status']=='EXITED' and process['exit_code']==0, 'failed process marked verified')
                require(process['memory_cap_bytes']==c['limits']['memory_bytes'] and process['cpu']==c['limits']['cpu'], 'changed measured resources')
                require(type(process['process_wall_ns']) is int and process['process_wall_ns']>0 and
                    process['process_wall_seconds']==process['process_wall_ns']/1e9, 'changed process clock units')
        else:
            require(all(r[k] is None for k in ('certificate','phase_costs','total_operations')), 'failed trial has verified costs')
        root = directory.parents[4]
        require(read(directory/'profile/process.json') == r['profile_process'], 'changed profile process record')
        if 'native_process' in r:
            require(read(directory/'native/process.json') == r['native_process'], 'changed native process record')
        require(set(r['artifacts']) == {str(p.relative_to(directory)) for p in directory.rglob('*')
                    if p.is_file() and p.name!='receipt.json'}, 'missing trial artifacts')
        record = scientific_trial(root, c, r, case, arm, directory)
        require(record == read(directory/'run.json') == r['measurement'], 'changed scientific run record')
    return r


def scientific_trial(root, c, row, case, arm, directory):
    stage = row['stage']
    relative = str(Path('admissions')/stage/case['id']/arm['id'])
    admitted = read(root/relative/'admission.json')
    require(objhash(admitted) == c['admissions'][relative], 'changed admission receipt')
    job = dict(copy.deepcopy(case['job']), mode='rho' if arm['id']=='rho' else 'ic', config=arm['config'])
    require(read(directory/'job.json') == job, 'changed executed job')
    manifest_path=root/arm['source_manifest_relative']
    check_admission(admitted, job=job, fixture=case['fixture'], manifest=read(manifest_path),
        metadata=read(manifest_path.parent/'producer.json'), resources=c['resources'],
        worker_sha256=digest(root/arm['binary_relative']))
    complete=row['status']=='VERIFIED'
    return run_record(admitted, number=c['run_number_base']+(STAGES.index(stage)*len(c['run_aliases'])+
        c['run_aliases'].index(arm['id']))*c['repetitions']+row['repetition'],
        host_id=objhash(c['host']), status='complete' if complete else
            {'TIMEOUT':'timeout','OOM':'oom'}.get(row['status'],'error'),
        native=read(directory/'native/stdout.json') if complete else None,
        profile=read(directory/'profile/stdout.json') if complete else None,
        costs=row['phase_costs'] if complete else None,
        process_wall_ns=row.get('native_process',{}).get('process_wall_ns'))


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
    if c.get('scientific_admission'):
        write(directory/'profile/process.json',run,exclusive=True)
    r = {'schema_version':2 if c.get('scientific_admission') else 1,'stage':stage,'arm':arm['id'],'case':case['id'],'cell':case['cell'],
         'case_sha256':objhash(case),'arm_sha256':objhash(arm),'repetition':repetition,
         'unit':c['unit'],'mode':job['mode'],'status':'ERROR','profile_process':run,
         'total_operations':None,'phase_costs':None,'certificate':None}
    try:
        require(run['process_status']!='TIMEOUT','TIMEOUT')
        require(run['exit_code']==0,'worker failed or incomplete')
        report = read(directory/'profile/stdout.json')
        proof = verify(report,case['fixture'],expected_mode=job['mode'],summands=job['config']['summands'])
        costs = parse_profiles(directory/'profile',phase_schema=3 if c.get('scientific_admission') else 1)
        expected = {'startup_and_input','curve_and_targets','final_verification','reporting_and_cleanup'}
        expected |= {'rho_solve'} if job['mode']=='rho' else PHASES-{'rho_solve'}
        if not c.get('scientific_admission'):
            require(set(costs)==expected,'missing exclusive phases')
        native = execute([str(binary)],job,directory/'native',lim['timeout_seconds'],lim['memory_bytes'],lim['cpu'])
        r['native_process'] = native
        if c.get('scientific_admission'):
            write(directory/'native/process.json',native,exclusive=True)
        require(native['exit_code']==0 and native['process_status']=='EXITED','native worker failed')
        nproof = verify(read(directory/'native/stdout.json'),case['fixture'],expected_mode=job['mode'],summands=job['config']['summands'])
        require(nproof['solutions']==proof['solutions'],'native/profile answers differ')
        total = sum(costs.values())
        r.update(status='VERIFIED',total_operations=total,phase_costs=costs,certificate=proof,
                 normalized_S=total/math.sqrt(int(case['fixture']['subgroup_order'])),
                 floor_operations=proof['rank'],
                 ratio_to_floor=total/proof['rank'] if proof['rank'] else None)
        if c.get('scientific_admission'):
            r['measurement']=scientific_trial(root,c,r,case,arm,directory)
            require(r['measurement']['total_operations']==total, 'scientific whole-cost mismatch')
    except (InvalidEvidence,ValueError,KeyError,TypeError) as exc:
        r.update(reason=str(exc), total_operations=None, phase_costs=None, certificate=None)
        for key in ('normalized_S','floor_operations','ratio_to_floor','measurement'):
            r.pop(key,None)
        if any(p.get('process_status')=='TIMEOUT' for p in (run,r.get('native_process',{}))):
            r['status']='TIMEOUT'
        elif any('memory allocation' in p.read_text() for p in directory.glob('*/stderr.txt')):
            r['status']='OOM'
        else:
            r['status']='INVALID_OR_INCOMPLETE'
    if c.get('scientific_admission'):
        if r['status']!='VERIFIED':
            r['measurement']=scientific_trial(root,c,r,case,arm,directory)
        write(directory/'run.json',r['measurement'],exclusive=True)
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


def comparison(rows, candidate_id, *, baseline='incumbent', draws=2000, match_support=True):
    selected=[r for r in rows if r['arm'] in (candidate_id,baseline)]
    cases={r['case'] for r in rows}
    grouped={}
    for r in selected:
        grouped.setdefault((r['case'],r['arm']),[]).append(r)
    logs={}
    native_logs={}
    online_logs={}
    scientific=any('measurement' in r for r in selected)
    for case in sorted(cases):
        a=grouped.get((case,baseline),[]); b=grouped.get((case,candidate_id),[])
        if not a or not b or len(a)!=len(b) or len({x['repetition'] for x in a})!=len(a) or len({x['repetition'] for x in b})!=len(b) or {x['repetition'] for x in a}!={x['repetition'] for x in b}:
            return {'candidate':candidate_id,'eligible':False,'reason':'missing paired trials'}
        if any(x['status']!='VERIFIED' or x['total_operations'] is None for x in a+b):
            return {'candidate':candidate_id,'eligible':False,'reason':'unverified or unpriced workload'}
        values=[v for x in a+b for v in (x['total_operations'],
                    x.get('native_process',{}).get('process_wall_seconds'))]
        if any(type(v) not in (int,float) or not math.isfinite(v) or v<=0 for v in values):
            return {'candidate':candidate_id,'eligible':False,'reason':'invalid or missing full cost'}
        require(len({x['case_sha256'] for x in a+b})==1,'changed paired fixtures')
        if candidate_id!='rho' and baseline!='rho':
            if match_support:
                require(len({x['certificate']['factor_base_sha256'] for x in a+b})==1,'changed paired factor-base support')
            else:
                require(all(len({x['certificate']['factor_base_sha256'] for x in arm})==1 for arm in (a,b)),
                        'changed factor-base support within a policy arm')
        cell=a[0]['cell']
        ratio=statistics.median(x['total_operations'] for x in b)/statistics.median(x['total_operations'] for x in a)
        logs.setdefault(cell,[]).append(math.log(ratio))
        wall_ratio=statistics.median(x['native_process']['process_wall_seconds'] for x in b)/statistics.median(x['native_process']['process_wall_seconds'] for x in a)
        native_logs.setdefault(cell,[]).append(math.log(wall_ratio))
        if scientific:
            values=[x.get('measurement',{}).get('native_timing') for x in a+b]
            if any(v is None or type(v['online']['wall_ns']) is not int or v['online']['wall_ns']<=0 for v in values):
                return {'candidate':candidate_id,'eligible':False,'reason':'missing verified online interval'}
            ratio=statistics.median(x['measurement']['native_timing']['online']['wall_ns'] for x in b)/statistics.median(x['measurement']['native_timing']['online']['wall_ns'] for x in a)
            online_logs.setdefault(cell,[]).append(math.log(ratio))
    require(bool(logs),'empty comparison')
    cell_values={cell:statistics.mean(xs) for cell,xs in logs.items()}
    estimate=math.exp(statistics.mean(cell_values.values()))
    rng=random.Random(751203)
    cells=sorted(logs)
    boot=[]
    native_boot=[]
    online_boot=[]
    for _ in range(draws):
        means=[]
        wall_means=[]
        online_means=[]
        for _ in cells:
            cell=rng.choice(cells)
            indices=rng.choices(range(len(logs[cell])),k=len(logs[cell]))
            means.append(statistics.mean(logs[cell][i] for i in indices))
            wall_means.append(statistics.mean(native_logs[cell][i] for i in indices))
            if scientific:
                online_means.append(statistics.mean(online_logs[cell][i] for i in indices))
        boot.append(math.exp(statistics.mean(means)))
        native_boot.append(math.exp(statistics.mean(wall_means)))
        if scientific:
            online_boot.append(math.exp(statistics.mean(online_means)))
    boot.sort()
    native_boot.sort()
    online_boot.sort()
    return {'candidate':candidate_id,'eligible':True,'candidate_over_baseline':estimate,
        'speedup':1/estimate,'ci95':[boot[int(.025*draws)],boot[min(draws-1,int(.975*draws))]],
        'per_cell':{cell:math.exp(v) for cell,v in cell_values.items()},
        'native_wall_candidate_over_baseline':math.exp(statistics.mean(statistics.mean(v) for v in native_logs.values())),
        'native_wall_ci95':[native_boot[int(.025*draws)],native_boot[min(draws-1,int(.975*draws))]],
        'native_wall_per_cell':{cell:math.exp(statistics.mean(v)) for cell,v in native_logs.items()},
        'native_wall_status':'paired cold-process timing; nested curve/target bootstrap',
        'paired_cases':len(cases),'independent_curve_blocks':len(cells),
        'online':dict(candidate_over_baseline=math.exp(statistics.mean(statistics.mean(v) for v in online_logs.values())),
            ci95=[online_boot[int(.025*draws)],online_boot[min(draws-1,int(.975*draws))]],
            per_cell={cell:math.exp(statistics.mean(v)) for cell,v in online_logs.items()},
            boundary='one supplied public point after reusable preparation through scalar replay') if scientific else None}


def gate(result,c):
    """Promotion over the incumbent.

    Objective `incumbent` (default): at least `1 - confirmation_ratio` lower
    cost with the paired interval below one and no cell worse than
    `max_cell_ratio`, in instructions and, when required, native wall.

    Objective `rho`: the incumbent gate only guards against regression — the
    candidate must be measurably cheaper than the incumbent (instruction
    ratio at most `no_regression_ratio` with the interval below one, native
    interval below one, no cell worse than `max_cell_ratio` in either
    metric) — and the decision additionally requires `rho_gate` on both final
    stages.  The margin keeps an A/A control from passing on noise.
    """
    if c.get('scientific_admission') and not c.get('reference_qualification'):
        return False  # Admission controls alone do not qualify a comparative reference.
    if not result.get('eligible'):
        return False
    cells_ok = max(result['per_cell'].values())<=c['max_cell_ratio']
    if c.get('objective','incumbent')=='rho':
        return bool(result['candidate_over_baseline']<=c.get('no_regression_ratio',0.98)
            and result['ci95'][1]<1 and cells_ok
            and result['native_wall_ci95'][1]<1
            and max(result['native_wall_per_cell'].values())<=c['max_cell_ratio'])
    return bool(result['candidate_over_baseline']<=c['confirmation_ratio']
        and result['ci95'][1]<1 and cells_ok
        and (not c.get('require_native_progress') or (result['native_wall_ci95'][1]<1
             and result['native_wall_candidate_over_baseline']<=c['confirmation_ratio']
             and max(result['native_wall_per_cell'].values())<=c['max_cell_ratio'])))


def rho_gate(paired):
    """Strictly below matched rho: both metrics' paired upper 95% limits and
    every curve cell below one.  A point estimate below one is not enough."""
    return bool(paired.get('eligible') and paired['ci95'][1]<1 and paired['native_wall_ci95'][1]<1
        and max(paired['per_cell'].values())<1 and max(paired['native_wall_per_cell'].values())<1)


def stage_arms(root, stage, arms):
    base=arms[0]
    contract=read(root/'contract.json') if (root/'contract.json').exists() else {}
    rho=contract.get('rho_reference',synthetic_arm('rho',base))
    if stage=='aa':
        return [base,synthetic_arm('aa_control',base)]
    if stage=='smoke':
        return arms+[rho]
    if stage=='development':
        smoke=read(root/'summaries/smoke.json')
        failed={r['arm'] for r in smoke['failures']}
        require('incumbent' not in failed,'incumbent failed smoke')
        return [a for a in arms if a['id'] not in failed]+[rho]
    if stage=='selection':
        d=read(root/'summaries/development.json')
        eligible=d.get('retained_portfolio')
        if eligible is None:
            eligible=sorted([r for r in d['comparisons'] if r.get('eligible')],key=lambda r:r['candidate_over_baseline'])[:2]
        return [base]+[next(a for a in arms if a['id']==r['candidate']) for r in eligible]+[rho]
    d=read(root/'summaries/selection.json')
    provisional=d.get('provisional_challenger')
    return [base]+([next(a for a in arms if a['id']==provisional)] if provisional else [])+[rho]


def summarize(root,c,stage,fixtures,arms,*,save=True):
    rows=load_stage(root,stage,fixtures,arms,c['repetitions'])
    comps=[comparison(rows,a['id'],draws=c['bootstrap_draws'],
                      match_support=c.get('comparison_kind','fixed-support')!='factor-base-policy') for a in arms if a['id'] not in ('incumbent','rho')]
    rho = comparison(rows,'rho',draws=c['bootstrap_draws']) if any(a['id']=='rho' for a in arms) else None
    result={'stage':stage,'runs':len(rows),'verified_runs':sum(r['status']=='VERIFIED' for r in rows),
            'comparisons':comps,'rho_over_incumbent':rho,
            'process_wall_seconds_including_profiling':sum(r['profile_process']['process_wall_seconds']+r.get('native_process',{}).get('process_wall_seconds',0) for r in rows),
            'failures':[{k:r.get(k) for k in ('case','arm','repetition','status','reason')} for r in rows if r['status']!='VERIFIED']}
    if c.get('scientific_admission'):
        result['single_target_online']=online_table(rows,fixtures,arms,c['repetitions'],
            ['rho'] if any(a['id']=='rho' for a in arms) else [])
        result['primary_metric']='single-target native online wall time'
        result['promotion_prerequisite']='Qualified IC and rho references plus the frozen familywise confirmation protocol; admission alone cannot promote.'
    if stage=='aa':
        aa=comps[0]
        result['passed']=bool(aa.get('eligible') and all(.95<=r<=1.05 for r in aa['per_cell'].values()) and not gate(aa,c))
    if stage=='development' and 'selection_width' in c:
        result['retained_portfolio']=retain(comps,arms,width=c['selection_width'],
            exploration=c['exploration_slots'],seed=c['seed'])
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
    challenger_over_rho={}
    if challenger and c.get('objective','incumbent')=='rho':
        for stage in ('confirmation','replay'):
            active=stage_arms(root,stage,all_arms)
            rows=load_stage(root,stage,fixtures[stage],active,c['repetitions'])
            challenger_over_rho[stage]=comparison(rows,challenger,baseline='rho',draws=c['bootstrap_draws'])
        passed=passed and all(rho_gate(p) for p in challenger_over_rho.values())
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
        if challenger_over_rho and not all(rho_gate(p) for p in challenger_over_rho.values()):
            reasons.append('The provisional challenger did not beat matched rho on both metrics with every upper 95% limit and every cell below one.')
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
    if choice:
        parity=[]
        for stage in ('confirmation','replay'):
            active=stage_arms(root,stage,all_arms)
            rows=load_stage(root,stage,fixtures[stage],active,c['repetitions'])
            paired=comparison(rows,choice,baseline='rho',draws=c['bootstrap_draws'])
            parity.append(paired)
        result['winner_over_rho']={'confirmation':parity[0],'replay':parity[1]}
        result['rho_parity']=all(p.get('eligible') and p['ci95'][1]<=c['parity_margin']
            and p['native_wall_ci95'][1]<=c['parity_margin']
            and max(p['per_cell'].values())<=c['parity_margin']
            and max(p['native_wall_per_cell'].values())<=c['parity_margin'] for p in parity)
        result['parity_definition']='Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.'
        result['beats_rho_strict']=all(rho_gate(p) for p in parity)
        result['beats_rho_strict_definition']='Winner/rho upper paired 95% limits and every cell ratio < 1 in both instructions and native process wall, on confirmation and replay. beats_rho alone is the confirmation point estimate of rho/winner exceeding one.'
    if challenger_over_rho:
        result['challenger_over_rho']=challenger_over_rho
    result['objective']=c.get('objective','incumbent')
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
    sources=arms+([c['rho_reference']] if 'rho_reference' in c else [])
    source_pairs={(a.get('source_manifest_relative','source-manifest.json'),a.get('source_directory','source')) for a in sources}
    source_files=0
    for manifest_path,source_path in source_pairs:
        source_manifest=read(root/manifest_path)
        source_files+=len(source_manifest)
        for name,h in source_manifest.items():
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
    print(json.dumps({'status':'VERIFIED','trial_receipts':count,'source_files':source_files}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('prepare')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--source-root',type=Path,default=ROOT)
    p.add_argument('--rho-source-root',type=Path,help='Freeze a separately qualified rho worker source.')
    p.add_argument('--rho-config',type=Path,help='Freeze the rho configuration selected on development data.')
    p.add_argument('--comparison-kind',choices=['fixed-support','factor-base-policy'],default='fixed-support')
    p.add_argument('--candidates',type=Path)
    p.add_argument('--profile',choices=['pilot','standard'],default='pilot')
    p.add_argument('--confirmation-cases',default='',
        help='cell=count,... raising named cells above the profile floor. Lowering is refused.')
    p.add_argument('--cells',default='13a0,17a1,19a0,23a0',
        help='development and selection curve cells, <degree>a<curve_a> comma separated')
    p.add_argument('--holdout-cells',default='19a1',
        help='cells added in confirmation and replay only; must not repeat a --cells entry')
    p.add_argument('--seed',type=int,default=20260915)
    p.add_argument('--cpu',type=int)
    p.add_argument('--timeout',type=float,default=60)
    p.add_argument('--max-processes',type=int,default=1800)
    p.add_argument('--require-native-progress',action='store_true')
    p.add_argument('--selection-width',type=int,default=6)
    p.add_argument('--exploration-slots',type=int,default=1)
    p.add_argument('--targets',type=int,default=1)
    p.add_argument('--objective',choices=['incumbent','rho'],default='incumbent',
                   help='rho: promote only a challenger that is measurably cheaper than the incumbent AND strictly below matched rho in both metrics (see gate/rho_gate)')
    p=commands.add_parser('propose')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--from-round',type=Path)
    for name in ('run','verify','status'):
        p=commands.add_parser(name);p.add_argument('--round',type=Path,required=True)
        if name=='run':p.add_argument('--stage',choices=['all']+STAGES,default='all')
    args=parser.parse_args()
    try:
        if args.command=='prepare':
            require(2<=args.selection_width<=15 and 0<=args.exploration_slots<args.selection_width,
                    'invalid selection/exploration budget')
            require(args.timeout>0 and args.max_processes>0,'positive limits required')
            require(1<=args.targets<=100,'target count must be 1..100')
            prepare(args)
        elif args.command=='propose':
            if args.from_round:
                proposed,source=proposals_from_previous(args.from_round.resolve())
                write(args.out,proposed,exclusive=True)
                parent_contract=read(args.from_round.resolve()/'contract.json')
                rho_reference=parent_contract.get('rho_reference',{})
                print(json.dumps({'candidates':str(args.out),'baseline_source_root':str(source),
                      'rho_source_root':str(args.from_round.resolve()/rho_reference.get('source_directory','source')),
                      'rho_config':rho_reference.get('config'),
                      'target_count':parent_contract.get('target_count',1),
                      'comparison_kind':parent_contract.get('comparison_kind','fixed-support'),
                      'require_native_progress':parent_contract.get('require_native_progress',False),
                      'next_step':'Use this source, preserve target count, comparison kind and metric gates, and choose a new campaign seed. A different target count is a separate workload panel.'}))
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
