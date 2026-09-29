"""Compare matched compiler-resource receipts; never infer device speed."""
import argparse
import json
from pathlib import Path
import fixtures as fx
from receipt_io import reserve


def compare(control,candidate):
    inputs=[json.loads(p.read_text()) for p in (control,candidate)]
    a,b=inputs
    if any(r['status']!='passed' or r['gpu_executed'] for r in inputs):
        raise ValueError('Expected two complete offline compiler receipts')
    for key in ('source_sha256','checker_sha256','ptxas_sha256','nvrtc_version'):
        if a[key]!=b[key]:raise ValueError(f'Unmatched {key}')
    if a['linear_square_no_unroll'] or not b['linear_square_no_unroll']:
        raise ValueError('Expected default control and no-unroll candidate')
    key=lambda c:(c['arch'],c['m'],c['fast_square'])
    left={key(c):c for c in a['compilations']};right={key(c):c for c in b['compilations']}
    expected={(arch,m,fast) for arch in ('compute_89','compute_120') for m in (83,131) for fast in (0,1)}
    if set(left)!=expected or set(right)!=expected or any(len(r['compilations'])!=8 for r in inputs):
        raise ValueError('Missing or duplicate compiler configurations')
    rows=[]
    for config in sorted(expected):
        baseline=left[config];trial=right[config]
        if trial['options']!=baseline['options']+['-DLINEAR_SQUARE_NOUNROLL=1']:
            raise ValueError('Compiler options differ by more than the frozen candidate')
        row={'arch':config[0],'m':config[1],'fast_square':config[2]}
        for label,value in (('control',baseline),('candidate',trial)):
            assembly=value['assembly']
            if value['returncode']!=0 or assembly['status']!='passed':
                raise ValueError('Compilation did not pass')
            keys=('registers_per_thread','stack_bytes','spill_store_bytes','spill_load_bytes')
            if any(assembly.get(k) is None for k in keys):raise ValueError('Missing resource statistics')
            row[label]={k:assembly[k] for k in keys}
        row['register_ratio']=row['candidate']['registers_per_thread']/row['control']['registers_per_thread']
        rows.append(row)
    return {'classification':'compiler-resource diagnostic; not a runtime measurement',
            'gpu_executed':False,'gpu_speedup':None,'rho_speedup':None,
            'inputs':{p.name:fx.sha(p.read_bytes()) for p in (control,candidate)},
            'source_sha256':a['source_sha256'],'rows':rows,
            'static_screen_passed':all(r['candidate']['spill_store_bytes']==0 and r['candidate']['spill_load_bytes']==0 for r in rows),
            'device_correctness':'pending','runtime_comparison':'pending'}


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--control',type=Path,required=True)
    ap.add_argument('--candidate',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();result=compare(args.control,args.candidate)
    reserve(args.output,result)
    print(json.dumps(result,indent=2))
