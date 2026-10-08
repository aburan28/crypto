"""Observe GPU compute processes while retaining ownership of the benchmark PID."""
from pathlib import Path
import datetime, json, os, signal, subprocess, time
from experiment import ROOT, OUT, sha


def finite(command, log, timeout=900, cwd=Path('/tmp'), env=None):
    start=time.monotonic();timed_out=False;errors=[];external=[];observations=0
    receipt=log.with_name(log.stem+'-pid.json')
    samples=log.with_name(log.stem+'-gpu-processes.jsonl')
    assert not log.exists() and not receipt.exists() and not samples.exists()
    with log.open('x') as stream,samples.open('x') as observed:
        child=subprocess.Popen(command,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        identity=dict(pid=child.pid,command=command,binary_sha256=sha(Path(command[0])),
                      launched_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                      source_sha256=sha(Path(__file__)))
        receipt.write_text(json.dumps(identity,indent=2)+'\n')
        while child.poll() is None:
            # No wait/poll/reap of this child occurs during the query. Its PID
            # cannot be recycled, even if it exits while nvidia-smi is running.
            row=dict(at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                     elapsed_seconds=time.monotonic()-start,owned_child_pid=child.pid,
                     child_alive_before_query=True,child_not_reaped_during_query=True)
            try:
                query=subprocess.run(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader'],capture_output=True,text=True,timeout=10)
                assert query.returncode==0,query.stderr
                processes=[dict(pid=int(line.split(',',1)[0]),path=line.split(',',1)[1].strip()) for line in query.stdout.splitlines() if line.strip()]
                row['processes']=processes
                outside=[p for p in processes if p['pid']!=child.pid]
                if outside:external.append(dict(at_utc=row['at_utc'],processes=outside))
                for process in processes:
                    if process['pid']==child.pid:
                        assert process['path'] in {command[0],'[No data]'},process
            except Exception as exc:
                row['observation_error']=repr(exc);errors.append(row.copy())
            observed.write(json.dumps(row)+'\n');observed.flush();observations+=1
            remaining=timeout-(time.monotonic()-start)
            if remaining<=0:
                timed_out=True;os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                break
            time.sleep(min(0.5,remaining))
    record=dict(passed=child.returncode==0 and not timed_out and not errors and not external and observations>0,
                pid=child.pid,observations=observations,nominal_sampling_interval_seconds=0.5,
                external_gpu_processes=external,observation_errors=errors,
                receipt=str(receipt.relative_to(OUT)),receipt_sha256=sha(receipt),
                samples=str(samples.relative_to(OUT)),samples_sha256=sha(samples),
                scope='Each query retains the launched child PID until completion. No process observation occurs after that child is reaped. Sampling cannot exclude activity between queries.',
                source_sha256=sha(Path(__file__)))
    return dict(command=command,returncode=child.returncode,timed_out=timed_out,timeout_seconds=timeout,
                elapsed_seconds=time.monotonic()-start,log=str(log.relative_to(OUT)),log_sha256=sha(log),gpu_observation=record)


def audit_observation(record):
    assert record['passed'] and record['observations']>0 and not record['external_gpu_processes'] and not record['observation_errors']
    assert record['source_sha256']==sha(Path(__file__))
    assert sha(OUT/record['receipt'])==record['receipt_sha256'] and sha(OUT/record['samples'])==record['samples_sha256']
    receipt=json.loads((OUT/record['receipt']).read_text())
    assert receipt['pid']==record['pid'] and receipt['source_sha256']==sha(Path(__file__))
    assert sha(Path(receipt['command'][0]))==receipt['binary_sha256']
    rows=[json.loads(line) for line in (OUT/record['samples']).read_text().splitlines()]
    assert len(rows)==record['observations']
    for row in rows:
        assert row['owned_child_pid']==record['pid'] and row['child_alive_before_query'] and row['child_not_reaped_during_query'] and 'observation_error' not in row
        assert all(p['pid']==record['pid'] and p['path'] in {receipt['command'][0],'[No data]'} for p in row['processes'])
    return True
