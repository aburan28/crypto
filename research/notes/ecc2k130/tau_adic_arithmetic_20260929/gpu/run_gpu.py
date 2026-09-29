"""One-shot CuPy runner for fixed synthetic known-scalar fixtures only."""
import argparse
import json
from pathlib import Path
import statistics
import sys
import time
import traceback
import numpy as np
import fixtures as fx
from receipt_io import reserve, save

VARIANTS=(('binary_naf',0),('reduced_tau_naf',0),('binary_naf',1),('reduced_tau_naf',1))


def label(v):return f'{v[0]}:linear_square={v[1]}'


def load_kernel(cp,modules,m,fast):
    start=time.perf_counter();key=(m,fast);was_cached=key in modules
    options=('--std=c++17',f'-DFIELD_M={m}',f'-DFAST_SQUARE={fast}')
    if not was_cached:
        modules[key]=cp.RawModule(code=(fx.HERE/'arithmetic.cuh').read_text(),options=options)
    kernel=modules[key].get_function('evaluate')
    return kernel,{'seconds':time.perf_counter()-start,'reused_module':was_cached,
                   'explicit_options':options,'attributes':kernel.attributes}


def check_device(cp,panels,modules,receipt,output):
    """Run the 96 unique archived cases in all four arms before any timing."""
    receipt['stage']='device_smoke';save(output,receipt)
    for panel in panels:
        m=panel['m'];prepared={method:fx.prepare(panel,method,repeat=1)
                              for method in ('binary_naf','reduced_tau_naf')}
        row={k:panel[k] for k in ('m','holdout','seed','identity','input_sha256','sage_output_sha256')}
        row.update(invocations={},unique_cases=len(panel['cases']))
        receipt['smoke']['panels'].append(row)
        dc=cp.asarray(fx.columns(m))
        for method,fast in VARIANTS:
            key=label((method,fast));hp,hd,hl,expected,meta=prepared[method]
            kernel,compiled=load_kernel(cp,modules,m,fast)
            dp,dd,dl=cp.asarray(hp),cp.asarray(hd),cp.asarray(hl)
            out=cp.empty_like(dp);n=len(hl)
            row['invocations'][key]={'status':'launch_pending','compile':compiled,
                                     'prepare':meta,'block_size':128,'evaluations':n}
            save(output,receipt)
            kernel(((n+127)//128,),(128,),
                   (dp,dd,dl,np.int32(n),dc,out,np.int32(method=='reduced_tau_naf')))
            receipt['gpu_launch_submitted']=True;save(output,receipt)
            cp.cuda.Stream.null.synchronize();receipt['gpu_executed']=True
            row['invocations'][key]['status']='verification_pending';save(output,receipt)
            got=out.get()
            if not np.array_equal(got,expected):
                row['invocations'][key]['status']='mismatch'
                raise AssertionError((m,key,'device smoke mismatch'))
            row['invocations'][key].update(status='passed',output_sha256=fx.sha(got.tobytes()))
            receipt['smoke']['verified_outputs']+=n;save(output,receipt)
    receipt['smoke']['status']='passed';save(output,receipt)


def hardware():
    return fx.ref.command(['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total,clocks.sm,clocks.mem,power.draw,temperature.gpu','--format=csv,noheader'])


def run(output,mode='benchmark'):
    if mode not in ('smoke','benchmark'):raise ValueError(mode)
    receipt={'status':'started','stage':'initializing','mode':mode,
             'gpu_executed':False,'gpu_launch_submitted':False,
             'classification':'known-scalar GPU arithmetic stage diagnostic',
             'smoke':{'status':'pending','verified_outputs':0,'panels':[]},
             'sources':{p.name:fx.sha(p.read_bytes()) for p in fx.HERE.iterdir() if p.suffix in ('.py','.cuh','.md')},
             'dependency_sources':{
                 'parent_benchmark.py':fx.sha((fx.PARENT/'benchmark.py').read_bytes()),
                 'curve_identity.py':fx.sha((fx.HERE.parents[4]/'tools/curve_identity.py').read_bytes()),
                 'parent_run-02.json':fx.sha((fx.PARENT/'results/run-02.json').read_bytes())},
             'rho_speedup':None,'walk_iterations_per_second':None,'panels':[]}
    reserve(output,receipt)
    try:
        import cupy as cp
        receipt['host']=fx.ref.host()
        receipt['before']=hardware()
        receipt['device_properties']=cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
        receipt['cuda_runtime']=cp.cuda.runtime.runtimeGetVersion()
        receipt['cuda_driver']=cp.cuda.runtime.driverGetVersion()
        receipt['cupy_version']=cp.__version__
        modules={};start_all=time.perf_counter()
        panels=list(fx.panels())
        receipt['oracle_seconds']=sum(p['oracle_seconds'] for p in panels)
        check_device(cp,panels,modules,receipt,output)
        receipt['stage']='benchmark' if mode=='benchmark' else 'device_smoke_complete'
        save(output,receipt)
        for panel in (panels if mode=='benchmark' else ()):
            m=panel['m'];prepared={};prep_meta={}
            for method in ('binary_naf','reduced_tau_naf'):
                arrays=fx.prepare(panel,method,repeat=1024)
                prepared[method]=arrays[:4];prep_meta[method]=arrays[4]
            start=time.perf_counter();cols=fx.columns(m);table_seconds=time.perf_counter()-start
            functions={};compile_meta={}
            for fast in (0,1):
                functions[fast],compile_meta[str(fast)]=load_kernel(cp,modules,m,fast)
            row={k:panel[k] for k in ('m','holdout','seed','identity','field','curve','source_receipt_sha256','input_sha256','sage_output_sha256','oracle_seconds')}
            row.update(prepare=prep_meta,table_seconds=table_seconds,compile=compile_meta,
                       block_size=128,unique_cases=24,evaluations_per_launch=24576,aa=[],rounds=[],invocations={})
            receipt['panels'].append(row);save(output,receipt)
            resident={}
            for v in VARIANTS:
                method,fast=v;hp,hd,hl,expected=prepared[method]
                start=time.perf_counter()
                dp,dd,dl=cp.asarray(hp),cp.asarray(hd),cp.asarray(hl)
                dc=cp.asarray(cols);out=cp.empty_like(dp);cp.cuda.Stream.null.synchronize()
                uploaded=time.perf_counter()
                args=(dp,dd,dl,np.int32(len(hl)),dc,out,np.int32(method=='reduced_tau_naf'))
                kernel=functions[fast];grid=((len(hl)+127)//128,)
                kernel(grid,(128,),args);cp.cuda.Stream.null.synchronize()
                launched=time.perf_counter();got=out.get();downloaded=time.perf_counter()
                if not np.array_equal(got,expected):raise AssertionError((m,label(v),'warm invocation'))
                verified=time.perf_counter();receipt['gpu_executed']=True
                row['invocations'][label(v)]={'upload_allocate_seconds':uploaded-start,
                    'kernel_host_seconds':launched-uploaded,'download_seconds':downloaded-launched,
                    'verification_seconds':verified-downloaded,'warm_gpu_invocation_seconds':downloaded-start,
                    'output_sha256':fx.sha(got.tobytes()),
                    'accounted_stage_seconds':prep_meta[method]['prepare_seconds']+table_seconds+compile_meta[str(fast)]['seconds']+verified-start}
                resident[v]=(kernel,grid,args,out,expected)
                save(output,receipt)

            def measure(v,repeats):
                kernel,grid,args,out,expected=resident[v]
                begin=cp.cuda.Event();end=cp.cuda.Event();begin.record()
                for _ in range(repeats):kernel(grid,(128,),args)
                end.record();end.synchronize();ms=cp.cuda.get_elapsed_time(begin,end)
                got=out.get()
                if not np.array_equal(got,expected):raise AssertionError((m,label(v),'timed result'))
                return {'event_ms':ms,'launches':repeats,'ms_per_launch':ms/repeats,
                        'known_scalar_evaluations_per_second':24576*repeats/(ms/1000),
                        'output_sha256':fx.sha(got.tobytes())}

            pilot=measure(VARIANTS[0],1)
            repeats=min(16,max(1,int(np.ceil(50/max(pilot['event_ms'],.001)))))
            row.update(pilot=pilot,launches_per_sample=repeats)
            for _ in range(5):
                row['aa'].append([measure(VARIANTS[0],repeats),measure(VARIANTS[0],repeats)])
                save(output,receipt)
            for rep in range(7):
                order=VARIANTS if rep%2==0 else VARIANTS[::-1]
                row['rounds'].append({'order':[label(v) for v in order],
                                      'samples':{label(v):measure(v,repeats) for v in order}})
                save(output,receipt)
            baseline=label(VARIANTS[0]);row['summary']={}
            noise=max(abs(a['event_ms']/b['event_ms']-1) for a,b in row['aa'])
            row['aa_max_deviation']=noise
            for v in VARIANTS:
                k=label(v);samples=[r['samples'][k]['ms_per_launch'] for r in row['rounds']]
                ratios=[r['samples'][k]['event_ms']/r['samples'][baseline]['event_ms'] for r in row['rounds']]
                ratio=statistics.median(ratios)
                row['summary'][k]={'median_ms_per_launch':statistics.median(samples),'min_ms_per_launch':min(samples),
                    'median_paired_kernel_cost_ratio':ratio,'screen_passed':ratio<=.90 and 1-ratio>noise,
                    'short_samples':any(r['samples'][k]['event_ms']<50 for r in row['rounds'])}
            row['after']=hardware();save(output,receipt)
            del resident
            cp.get_default_memory_pool().free_all_blocks()
        receipt['total_benchmark_seconds']=time.perf_counter()-start_all
        receipt['status']='passed';receipt['stage']='complete';receipt['after']=hardware()
    except Exception:
        receipt['status']='failed';receipt['error']=traceback.format_exc()
        if receipt['stage']=='device_smoke':receipt['smoke']['status']='failed'
    save(output,receipt)
    print(json.dumps({'status':receipt['status'],'gpu_executed':receipt['gpu_executed'],'output':str(output)}))
    return 0 if receipt['status']=='passed' else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--mode',choices=('smoke','benchmark'),default='benchmark')
    args=ap.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    raise SystemExit(run(args.output,args.mode))
