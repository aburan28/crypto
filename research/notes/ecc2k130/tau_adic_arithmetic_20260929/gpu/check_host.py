"""CPU compilation + independent arithmetic validation of the CUDA header."""
import argparse
import ctypes as C
import json
from pathlib import Path
import random
import subprocess
import tempfile
import time
import traceback
import numpy as np
import fixtures as fx


def main(output):
    with output.open('x') as f:f.write('{}\n')
    result={'status':'started','gpu_executed':False,'source_sha256':fx.sha((fx.HERE/'arithmetic.cuh').read_bytes()),
            'field_checks':0,'point_checks':0,'panels':[],'compiler':subprocess.check_output(['g++','--version'],text=True)}
    try:
        panels=list(fx.panels())
        with tempfile.TemporaryDirectory(prefix='tau-host-') as folder:
            for m in (83,131):
                field=fx.Field(m);w=(m+63)//64;columns=fx.columns(m)
                for fast in (0,1):
                    so=Path(folder)/f'field{m}-{fast}.so'
                    cmd=['g++','-std=c++17','-O2','-shared','-fPIC','-x','c++',f'-DFIELD_M={m}',f'-DFAST_SQUARE={fast}',str(fx.HERE/'arithmetic.cuh'),'-o',str(so)]
                    subprocess.run(cmd,check=True,capture_output=True,text=True)
                    lib=C.CDLL(str(so))
                    lib.host_evaluate.argtypes=[C.c_void_p,C.c_void_p,C.c_void_p,C.c_int,C.c_void_p,C.c_void_p,C.c_int]
                    lib.host_field.argtypes=[C.c_void_p]*4
                    rng=random.Random(20260929+m)
                    for a in [0,1,(1<<m)-1,*[rng.getrandbits(m) for _ in range(64)]]:
                        words=np.array([(a>>(64*j))&((1<<64)-1) for j in range(w)],dtype=np.uint64)
                        sq=np.zeros(w,dtype=np.uint64);inv=sq.copy()
                        lib.host_field(words.ctypes.data,columns.ctypes.data,sq.ctypes.data,inv.ctypes.data)
                        decode=lambda x:sum(int(v)<<(64*j) for j,v in enumerate(x))
                        assert decode(sq)==field.mul(a,a)
                        assert decode(inv)==(field.inv(a) if a else 0)
                        result['field_checks']+=2
                    chosen=[p for p in panels if p['m']==m]
                    point=chosen[0]['cases'][0][0]
                    edges=[(p,k) for p in (None,(0,1),point) for k in (0,1,2,3,-1)]
                    edge={'m':m,'seed':20260929,'cases':edges,'expected':[field.multiply(p,k) for p,k in edges]}
                    for panel in chosen+[edge]:
                        for method in ('binary_naf','reduced_tau_naf'):
                            points,ds,lengths,expected,meta=fx.prepare(panel,method)
                            out=np.empty_like(points)
                            lib.host_evaluate(points.ctypes.data,ds.ctypes.data,lengths.ctypes.data,len(lengths),columns.ctypes.data,out.ctypes.data,int(method=='reduced_tau_naf'))
                            assert np.array_equal(out,expected),(m,fast,method)
                            result['point_checks']+=len(lengths)
                    result['panels'].append({'m':m,'fast_square':fast,'status':'passed'})
        result['status']='passed'
    except Exception:
        result['status']='failed';result['error']=traceback.format_exc()
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('compiler',)}))
    return 0 if result['status']=='passed' else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    raise SystemExit(main(args.output))
