"""Freeze targets, then run exact coverage, matched solvers and complete DLP sweeps."""
import argparse
import importlib.util
import json
import hashlib
from pathlib import Path
import platform
import random
import subprocess
import sys
import time
import traceback

import core

HERE=Path(__file__).resolve().parent
ROOT=core.ROOT


def readRows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def freeze(contract):
    excluded={}
    files=list(HERE.parent.glob('solver_*/raw.jsonl'))+[HERE.parent/'coefficient_pullback/results/raw.jsonl',HERE.parent/'pullback_incremental/results_v2/raw.jsonl']
    for path in files:
        for row in readRows(path):
            if row.get('kind') in ('trial','stage'):
                excluded.setdefault(row['n'],set()).add(tuple(row['target']))
    out=[]
    for profile in contract['profiles']:
        truth=core.Truth(profile);f=truth.f
        points=sorted(([f.toCoords(x) for x in p] for p in truth.points))
        if profile['curve_b_coords']==(1<<profile['bits'])-1:
            points=[p for p in points if tuple(p) not in excluded.get(profile['bits'],set())]
        rng=random.Random(contract['target_seed']+profile['bits']+profile['curve_b_coords'])
        selected=rng.sample(points,contract['stage_targets_per_profile'])
        out.append({'profile':profile['name'],'curve_order':truth.full,'subgroup_order':truth.order,
                    'targets':[{'target':p,'stratum':'development' if i<contract['stage_development_targets'] else 'holdout'}
                               for i,p in enumerate(selected)]})
    return {'contract_sha256':hashlib.sha256((HERE/'contract.json').read_bytes()).hexdigest(),'profiles':out}


def checkWitness(truth,ell,k,targetCoords,result):
    f=truth.f;curve=truth.curve;target=tuple(f.fromCoords(v) for v in targetCoords)
    expected=truth.census(ell,k)[1][target]
    witness=result.get('witness')
    if witness is None:
        if result['status']=='empty' and expected:raise ArithmeticError('false empty decomposition')
        return
    if not expected:raise ArithmeticError('unexpected decomposition')
    points=[tuple(f.fromCoords(v) for v in p) for p in witness['points']]
    allowed=set(core.support(f,truth.profile,ell));xs=[p[0] for p in points]
    if len(points)!=k or len(set(xs))!=k or any(not x or x==target[0] or x not in allowed for x in xs):
        raise ArithmeticError('invalid support certificate')
    total=None
    for p in points:
        if not curve.onCurve(p):raise ArithmeticError('off-curve certificate')
        total=curve.add(total,p)
    if total!=target:raise ArithmeticError('wrong certificate sum')


def validate(contract):
    profile=contract['profiles'][0];truth=core.Truth(profile);checks=0
    for ell in profile['dimensions']:
        for k in contract['summands']:
            for variant in contract['stage_variants']:
                f,curve=core.fieldCurve(profile)
                search=core.Search(f,curve,profile,ell,k,variant,time.perf_counter()+120)
                for point in truth.points:
                    coords=[truth.f.toCoords(x) for x in point]
                    found=search.find(tuple(f.fromCoords(v) for v in coords))
                    checkWitness(truth,ell,k,coords,{'status':'found' if found else 'empty','witness':core.certificate(f,found)})
                    checks+=1
    subfieldChecks=0
    for profile in contract['profiles'][2:]:
        f,curve=core.fieldCurve(profile);scalars=[f.fromCoords(x) for x in range(512) if f.frob(f.fromCoords(x),3)==f.fromCoords(x)]
        if f.frob(curve.b,3)!=curve.b:raise ArithmeticError('curve coefficient outside GF8')
        for ell in profile['dimensions']:
            values=set(core.support(f,profile,ell))
            if len(values)!=8**ell:raise ArithmeticError('wrong subfield-space size')
            for a in scalars:
                for x in values:
                    if f.mul(a,x) not in values:raise ArithmeticError('space not GF8-linear')
            subfieldChecks+=1
    return {'exhaustive_five_bit_solver_checks':checks,'gf8_support_closure_checks':subfieldChecks,
            'all_passed':True}


def algebraicCell(profile,ell,target,variant,budget,truth):
    f=truth.f;point=tuple(f.fromCoords(v) for v in target)
    oracle=core.prior.previous.PairOracle(profile['bits'],ell)
    expected=oracle.expected(point)
    if variant=='cached-pullback':
        spec=importlib.util.spec_from_file_location('joint_pullback_control',HERE.parent/'pullback_incremental/optimized.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        row=module.cell(profile['bits'],ell,target,'first',budget)
    else:
        # Expected sets belong only to the post-solve check in the inherited adapter.
        row=core.prior.previous.s4.runCell(oracle.f,oracle.curve,ell,point,expected,variant,'first',budget)
    if not {tuple(xs) for xs in row['solutions']}<=expected:
        raise ArithmeticError('algebraic control returned an invalid exact point decomposition')
    counts=truth.census(ell,3)[1]
    if row['solutions'] and not counts[point]:raise ArithmeticError('algebraic control oracle mismatch')
    if row['status']=='complete' and not row['solutions'] and counts[point]:raise ArithmeticError('algebraic false empty')
    return {'profile':profile['name'],'ell':ell,'k':3,'variant':variant,'target':target,
            'status':'found' if row['solutions'] else ('empty' if row['status']=='complete' else 'timeout'),
            'cold_seconds':row['all_phase_seconds'],'within_budget':row['all_phase_seconds']<=budget,'adapter_result':row,
            'field_api_counts':row.get('field_api_counts'),'full_dlp_S':None,'rho_ratio':None}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',type=Path);p.add_argument('--targets',type=Path)
    p.add_argument('--output',type=Path);p.add_argument('--validate-only',action='store_true');args=p.parse_args()
    contract=json.loads((HERE/'contract.json').read_text())
    if args.freeze:
        frozen=freeze(contract)
        with args.freeze.open('x') as out:json.dump(frozen,out,indent=2);out.write('\n')
        print('Frozen profiles:',[(r['profile'],r['curve_order'],r['subgroup_order']) for r in frozen['profiles']]);return
    if args.targets is None or args.output is None:p.error('--targets and --output required')
    targets=json.loads(args.targets.read_text())
    if targets['contract_sha256']!=hashlib.sha256((HERE/'contract.json').read_bytes()).hexdigest():raise ValueError('contract changed after freeze')
    args.output.mkdir(parents=True,exist_ok=False)
    paths=list(HERE.glob('*.py'))+[HERE/'contract.json',args.targets]+list(core.prior.CODE.glob('*.py'))
    for folder in ['coefficient_pullback','pullback_incremental','solver_02','solver_04','solver_05','solver_06','solver_07','solver_08']:
        paths+=list((HERE.parent/folder).glob('*.py'))
    provenance={'contract':contract,'targets':targets,'commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                'python':platform.python_version(),'platform':platform.platform(),'pycryptosat':core.prior.importlib.metadata.version('pycryptosat'),
                'source_sha256':{str(path.resolve().relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                'command':[sys.executable]+sys.argv}
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    errors=[];totals={'stage':0,'algebraic':0,'dlp':0,'verified_dlp':0,'incomplete_dlp':0};rng=random.Random(contract['target_seed'])
    with (args.output/'raw.jsonl').open('x') as out:
        def record(row):out.write(json.dumps(row)+'\n');out.flush()
        try:validation=validate(contract);record({'kind':'validation',**validation});print(validation,flush=True)
        except Exception:record({'kind':'error','phase':'validation','traceback':traceback.format_exc()});raise
        if not args.validate_only:
            for profile in contract['profiles']:
                started=time.perf_counter();truth=core.Truth(profile)
                for ell in profile['dimensions']:
                    for k in contract['summands']:record({'kind':'census',**truth.census(ell,k)[0]})
                record({'kind':'external_validation_cost','profile':profile['name'],'seconds':time.perf_counter()-started,'field_api_counts':truth.f.report()})
                selected=next(p for p in targets['profiles'] if p['profile']==profile['name'])['targets']
                for ell in profile['dimensions']:
                    for k in contract['summands']:
                        for entry in selected:
                            variants=list(contract['stage_variants'])
                            if k==3 and profile['q_bits']==1:variants+=['cached-pullback','s4-symmetric','chained-s3']
                            rng.shuffle(variants)
                            for variant in variants:
                                try:
                                    if variant in contract['stage_variants']:
                                        row=core.stageCell(profile,ell,k,entry['target'],variant,contract['stage_seconds'])
                                        checkWitness(truth,ell,k,entry['target'],row);kind='stage'
                                    else:
                                        row=algebraicCell(profile,ell,entry['target'],variant,contract['stage_seconds'],truth);kind='algebraic'
                                    record({'kind':kind,'stratum':entry['stratum'],**row});totals[kind]+=1
                                except Exception:
                                    error={'kind':'error','profile':profile['name'],'ell':ell,'k':k,'variant':variant,**entry,'traceback':traceback.format_exc()}
                                    errors.append(error);record(error)
                        print('stage',profile['name'],ell,k,flush=True)
                        for seed in contract['e2e']['seeds']:
                            variants=list(contract['e2e']['variants']);rng.shuffle(variants)
                            for variant in variants:
                                try:
                                    row=core.dlp(profile,ell,k,variant,seed,contract['e2e'])
                                    for attempt in row['attempts']:checkWitness(truth,ell,k,attempt['target'],attempt)
                                    if row['status']=='verified':
                                        G=tuple(truth.f.fromCoords(v) for v in row['generator']);Q=tuple(truth.f.fromCoords(v) for v in row['target'])
                                        if truth.curve.mul(G,row['recovered_scalar'])!=Q:raise ArithmeticError('independent DLP check failed')
                                        totals['verified_dlp']+=1
                                    else:totals['incomplete_dlp']+=1
                                    totals['dlp']+=1;record({'kind':'dlp',**row})
                                    print('dlp',profile['name'],ell,k,seed,variant,row['status'],row.get('error'),flush=True)
                                except Exception:
                                    error={'kind':'error','profile':profile['name'],'ell':ell,'k':k,'seed':seed,'variant':variant,'traceback':traceback.format_exc()}
                                    errors.append(error);record(error)
                for seed in contract['e2e']['seeds']:
                    row=core.dlp(profile,None,None,'rho',seed,contract['e2e'])
                    if row['status']!='verified':raise ArithmeticError('rho reference failed')
                    record({'kind':'rho',**row})
    (args.output/'summary.json').write_text(json.dumps({'counts':totals,'errors':errors,'validation':validation},indent=2)+'\n')
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
