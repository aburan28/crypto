#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess
from client_process import capture
from experiment import ROOT,OUT,sha

def run(binary,batch,tag,checkpoint,launches,run_id):
    command=[str(binary),'--packed','--threads',str(8064//batch),'--steps','64',
             '--launches',str(launches),'--run-id',str(run_id),'--checkpoint',str(checkpoint),
             '--bench','--verify','0']
    record=capture(command,OUT/f'{tag}.log')
    assert record.returncode==0 and 'MISMATCH' not in record.stdout and '0 dropped' in record.stdout
    return hashlib.sha256(checkpoint.read_bytes()).hexdigest(),record.stdout

rows=[]
for batch in (24,32):
    binary=ROOT/'build'/f'g7-split-batch{batch}';run_id=58231+batch
    uninterrupted=ROOT/'build'/f'batch{batch}-uninterrupted.ckpt';uninterrupted.unlink(missing_ok=True)
    expected,_=run(binary,batch,f'batch{batch}-uninterrupted',uninterrupted,4,run_id)
    resumed=ROOT/'build'/f'batch{batch}-resumed.ckpt';resumed.unlink(missing_ok=True)
    run(binary,batch,f'batch{batch}-resume-start',resumed,2,run_id)
    actual,text=run(binary,batch,f'batch{batch}-resume-end',resumed,2,run_id)
    assert expected==actual and 'resumed from' in text and 'at iteration 128' in text
    rows.append({'batch':batch,'passed':True,'final_checkpoint_sha256':actual,'final_iteration':256})

source=ROOT/'build/batch16-rejection-source.ckpt';source.unlink(missing_ok=True)
run(ROOT/'build/g7-split-batch16',16,'batch16-rejection-source',source,1,58301)
for batch in (24,32):
    log=OUT/f'batch{batch}-reject-batch16.log'
    command=[str(ROOT/'build'/f'g7-split-batch{batch}'),'--packed','--threads',str(8064//batch),
             '--steps','1','--launches','1','--run-id','58301','--checkpoint',str(source),
             '--bench','--verify','0']
    record=capture(command,log)
    assert record.returncode==6 and 'incompatible or incomplete' in record.stdout
    rows.append({'batch':batch,'rejects_batch16_checkpoint':True,'returncode':record.returncode,'log_sha256':sha(log)})
(OUT/'persistence.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
