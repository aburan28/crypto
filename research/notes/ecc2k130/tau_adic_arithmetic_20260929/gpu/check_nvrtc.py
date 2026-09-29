"""Compile the actual CUDA kernels to PTX without a GPU; not device execution."""
import argparse
import ctypes as C
import json
from pathlib import Path
import sys
import traceback
import fixtures as fx


def main(output):
    with output.open('x') as f:f.write('{}\n')
    result={'status':'started','gpu_executed':False,'compilations':[],
            'source_sha256':fx.sha((fx.HERE/'arithmetic.cuh').read_bytes())}
    try:
        paths=list(Path(sys.prefix).glob('lib/python*/site-packages/nvidia/cuda_nvrtc/lib/libnvrtc.so*'))
        if not paths:raise RuntimeError('Install nvidia-cuda-nvrtc-cu12==12.8.93')
        nv=C.CDLL(str(paths[0]));program=C.c_void_p
        nv.nvrtcCreateProgram.argtypes=[C.POINTER(program),C.c_char_p,C.c_char_p,C.c_int,C.c_void_p,C.c_void_p]
        nv.nvrtcCompileProgram.argtypes=[program,C.c_int,C.POINTER(C.c_char_p)]
        nv.nvrtcGetProgramLogSize.argtypes=[program,C.POINTER(C.c_size_t)]
        nv.nvrtcGetProgramLog.argtypes=[program,C.c_void_p]
        nv.nvrtcGetPTXSize.argtypes=[program,C.POINTER(C.c_size_t)]
        nv.nvrtcGetPTX.argtypes=[program,C.c_void_p]
        nv.nvrtcDestroyProgram.argtypes=[C.POINTER(program)]
        major=C.c_int();minor=C.c_int();nv.nvrtcVersion(C.byref(major),C.byref(minor))
        result['nvrtc_version']=[major.value,minor.value]
        for arch in ('compute_89','compute_120'):
            for m in (83,131):
                for fast in (0,1):
                    p=program();src=(fx.HERE/'arithmetic.cuh').read_bytes()
                    assert nv.nvrtcCreateProgram(C.byref(p),src,b'arithmetic.cuh',0,None,None)==0
                    try:
                        opts=[b'--std=c++17',f'--gpu-architecture={arch}'.encode(),f'-DFIELD_M={m}'.encode(),f'-DFAST_SQUARE={fast}'.encode()]
                        rc=nv.nvrtcCompileProgram(p,len(opts),(C.c_char_p*len(opts))(*opts))
                        size=C.c_size_t();nv.nvrtcGetProgramLogSize(p,C.byref(size))
                        log=C.create_string_buffer(size.value);nv.nvrtcGetProgramLog(p,log)
                        row={'arch':arch,'m':m,'fast_square':fast,'returncode':rc,'log':log.value.decode()}
                        result['compilations'].append(row)
                        assert rc==0,row
                        assert nv.nvrtcGetPTXSize(p,C.byref(size))==0
                        ptx=C.create_string_buffer(size.value);assert nv.nvrtcGetPTX(p,ptx)==0
                        row.update(ptx_bytes=size.value,ptx_sha256=fx.sha(ptx.value))
                    finally:nv.nvrtcDestroyProgram(C.byref(p))
        result['status']='passed'
    except Exception:
        result['status']='failed';result['error']=traceback.format_exc()
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return 0 if result['status']=='passed' else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    raise SystemExit(main(args.output))
