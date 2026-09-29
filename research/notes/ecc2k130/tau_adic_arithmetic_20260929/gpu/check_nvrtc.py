"""Compile the actual CUDA kernels to PTX without a GPU; not device execution."""
import argparse
import ctypes as C
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import traceback
import fixtures as fx
from receipt_io import reserve, save


def assemble_ptx(ptx, arch, ptxas):
    with tempfile.TemporaryDirectory(prefix='tau-ptxas-') as folder:
        source=Path(folder)/'kernel.ptx';binary=Path(folder)/'kernel.cubin'
        source.write_bytes(ptx)
        cmd=[str(ptxas),'--verbose',f'--gpu-name={arch.replace("compute_", "sm_")}',str(source),'-o',str(binary)]
        try:
            proc=subprocess.run(cmd,capture_output=True,text=True,timeout=60)
        except subprocess.TimeoutExpired as error:
            decode=lambda s:s.decode(errors='replace') if isinstance(s,bytes) else (s or '')
            return {'status':'timeout','timeout_seconds':60,
                    'stdout':decode(error.stdout),'stderr':decode(error.stderr)}
        log=proc.stdout+proc.stderr
        registers=re.search(r'Used (\d+) registers',log)
        stack=re.search(r'(\d+) bytes stack frame, (\d+) bytes spill stores, (\d+) bytes spill loads',log)
        row={'status':'passed' if proc.returncode==0 else 'failed','returncode':proc.returncode,
             'options':cmd[1:3],'stdout':proc.stdout,'stderr':proc.stderr,
             'registers_per_thread':int(registers[1]) if registers else None,
             'stack_bytes':int(stack[1]) if stack else None,
             'spill_store_bytes':int(stack[2]) if stack else None,
             'spill_load_bytes':int(stack[3]) if stack else None}
        if proc.returncode==0:
            data=binary.read_bytes();row.update(cubin_bytes=len(data),cubin_sha256=fx.sha(data))
        return row


def main(output, assemble=False):
    result={'status':'started','gpu_executed':False,'compilations':[],
            'source_sha256':fx.sha((fx.HERE/'arithmetic.cuh').read_bytes()),
            'checker_sha256':fx.sha(Path(__file__).read_bytes()),'assemble':assemble}
    reserve(output,result)
    try:
        ptxas=None
        if assemble:
            paths=list(Path(sys.prefix).glob('lib/python*/site-packages/nvidia/cuda_nvcc/bin/ptxas'))
            if not paths:raise RuntimeError('Install nvidia-cuda-nvcc-cu12==12.8.93')
            ptxas=paths[0]
            result['ptxas_version']=subprocess.check_output([str(ptxas),'--version'],text=True)
            result['ptxas_sha256']=fx.sha(ptxas.read_bytes())
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
                        if assemble:
                            row['assembly']=assemble_ptx(ptx.value,arch,ptxas)
                            if row['assembly']['status']!='passed':
                                raise RuntimeError(f'ptxas failed for {arch}, m={m}, fast_square={fast}')
                        save(output,result)
                    finally:nv.nvrtcDestroyProgram(C.byref(p))
        result['status']='passed'
    except Exception:
        result['status']='failed';result['error']=traceback.format_exc()
    save(output,result);print(json.dumps(result))
    return 0 if result['status']=='passed' else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--assemble',action='store_true',help='Also assemble PTX and retain register/spill reports; no GPU needed')
    args=ap.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    raise SystemExit(main(args.output,args.assemble))
