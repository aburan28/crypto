from client_process import capture
"""Compare complete CUDA-client state, DP output and mixed-binary resume."""
import argparse
import array
import struct
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

def verify(reference,candidate,label,candidate_batch=16,run_ids=(0,139),resume_run_id=142,logical_walks=8192):
    assert logical_walks>0 and logical_walks%16==0 and logical_walks%candidate_batch==0
    reference_workers=logical_walks//16
    paths={'reference':str(Path(reference).resolve()),'candidate':str(Path(candidate).resolve())}
    results=[]
    def state_image(checkpoint):
        data=checkpoint.read_bytes()
        magic,version,m,threads,batch,lanes,runid,iters=struct.unpack_from('<8s6IQ',data)
        assert (magic,version,m,lanes)==(b'ECC2K130',2,131,1)
        count=threads*batch
        assert count==logical_walks and len(data)==40+count*60
        if batch==16: return data
        # Normalize to the reference's 16-slot checkpoint layout.
        result=bytearray(struct.pack('<8s6IQ',magic,version,m,reference_workers,16,lanes,runid,iters))
        offset=40
        for _ in range(2):
            source=array.array('I');source.frombytes(data[offset:offset+count*20])
            dest=array.array('I',[0])*(count*5)
            for i in range(count):
                for w in range(5):
                    dest[((i//reference_workers)*5+w)*reference_workers+i%reference_workers]=source[((i//threads)*5+w)*threads+i%threads]
            result.extend(dest.tobytes());offset+=count*20
        result.extend(data[offset:])
        return bytes(result)
    def run(name,tag,checkpoint,launches=4,runid=0,bench=False):
        corpus=ROOT/'build'/f'{label}-{tag}.dp'
        corpus.unlink(missing_ok=True)
        command=[paths[name],'--packed','--threads',str(logical_walks//(candidate_batch if name=='candidate' else 16)),'--steps','64','--launches',str(launches),
            '--run-id',str(runid),'--checkpoint',str(checkpoint)]
        command+=['--bench','--verify','0'] if bench else ['--dp-weight','52','--verify','16','--dp-file',str(corpus)]
        p=capture(command,OUT/f'{label}-{tag}.log')
        assert p.returncode==0 and 'MISMATCH' not in p.stdout,p.stdout
        assert '0 dropped' in p.stdout
        if not bench: assert '16 verified against the reference' in p.stdout
        data=corpus.read_bytes() if not bench else b''
        assert len(data)%32==0
        return hashlib.sha256(state_image(checkpoint)).hexdigest(),sorted(data[i:i+32] for i in range(0,len(data),32))
    for runid in run_ids:
        checks=[]
        for name in paths:
            checkpoint=ROOT/'build'/f'{label}-{name}-{runid}.ckpt'
            checkpoint.unlink(missing_ok=True)
            checks.append(run(name,f'{name}-{runid}',checkpoint,runid=runid))
        assert checks[0]==checks[1]
        results.append(dict(run_id=runid,checkpoint_sha256=checks[0][0],matching_dp_records=len(checks[0][1]),reference_replays_per_binary=16))
    # Mixed-binary continuation in both directions, compared with uninterrupted state.
    expected=ROOT/'build'/f'{label}-uninterrupted.ckpt';expected.unlink(missing_ok=True)
    end,_=run('reference','uninterrupted',expected,runid=resume_run_id,bench=True)
    for first,second in ((('reference','candidate'),('candidate','reference')) if candidate_batch==16 else ()):
        checkpoint=ROOT/'build'/f'{label}-resume-{first}.ckpt';checkpoint.unlink(missing_ok=True)
        run(first,f'resume-{first}-start',checkpoint,launches=2,runid=resume_run_id,bench=True)
        digest,_=run(second,f'resume-{first}-end',checkpoint,launches=2,runid=resume_run_id,bench=True)
        assert digest==end
    result=dict(logical_walks=logical_walks,validation=results,bidirectional_resume=(candidate_batch==16),candidate_batch=candidate_batch,binary_sha256={k:hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in paths.items()})
    (OUT/f'{label}-validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('reference');p.add_argument('candidate');p.add_argument('label');p.add_argument('--candidate-batch',type=int,default=16)
    p.add_argument('--logical-walks',type=int,default=8192)
    verify(**vars(p.parse_args()))
