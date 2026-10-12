"""Retain and transport-replay a terminal SAT v3 development control once."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import tarfile

from audit_static_sat_full_v3 import audit
from identity import sha256, write_immutable
from oracle import require
from sat_runtime_execution_v3 import read


def digest(data):
    return hashlib.sha256(data).hexdigest()


def collect(root,prefix):
    root=Path(root).resolve()
    files={}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(),'symlinked SAT control publication input')
        if path.is_file():
            require(path.resolve().is_relative_to(root),'SAT publication input escapes root')
            role=prefix+'/'+path.relative_to(root).as_posix()
            files[role]=(path.read_bytes(),path.stat().st_mode&0o777)
    require(files,'empty SAT control publication input')
    return files


def inventory(files):
    return {role:dict(bytes=len(data),sha256=digest(data),mode=mode)
            for role,(data,mode) in sorted(files.items())}


def publish(registration,execution,audit_file,output):
    registration,execution,output=map(Path,(registration,execution,output))
    require(not output.exists(),'SAT control publication exists; never overwrite evidence')
    spec=read(registration/'execution.json')
    result=audit(execution,spec)
    require(result==read(audit_file),'SAT control independent receipt differs from fresh replay')
    require(spec['arguments']['panel']['question']=='development-source-control',
            'this publisher does not qualify a paired tournament result')
    files=collect(registration,'registration')
    files.update(collect(execution,'execution'))
    files['independent-audit.json']=(Path(audit_file).read_bytes(),0o444)
    files['publisher.py']=(Path(__file__).read_bytes(),0o444)
    output.mkdir(parents=True)
    archive=output/'evidence.tar.gz'
    with archive.open('xb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as compressed:
        with tarfile.open(fileobj=compressed,mode='w|') as tar:
            for role,(data,mode) in sorted(files.items()):
                item=tarfile.TarInfo(role)
                item.size,item.mtime,item.mode=len(data),0,mode
                tar.addfile(item,io.BytesIO(data))
    receipt=dict(schema_version=3,execution_sha256=sha256(spec),
                 candidate_id=result['candidate_id'],workload_id=result['workload_id'],
                 run_id=result['run_id'],archive_sha256=digest(archive.read_bytes()),
                 archive_bytes=archive.stat().st_size,inventory=inventory(files),
                 result=result,promotion_eligible=False,online_speedup=None)
    write_immutable(output/'receipt.json',receipt)
    write_immutable(output/'AUDIT.json',result)
    return receipt


def replay(bundle,output,expected_execution_sha256):
    bundle,output=map(Path,(bundle,output))
    receipt=read(bundle/'receipt.json')
    require(receipt['execution_sha256']==expected_execution_sha256,
            'SAT control differs from externally frozen invocation')
    data=(bundle/'evidence.tar.gz').read_bytes()
    require(digest(data)==receipt['archive_sha256'] and len(data)==receipt['archive_bytes'],
            'SAT control published archive changed')
    require(not output.exists(),'SAT control replay output exists')
    files={}
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as tar:
        for item in tar:
            role=PurePosixPath(item.name)
            require(item.isfile() and not role.is_absolute() and '..' not in role.parts
                    and role.as_posix()==item.name and item.name not in files,
                    'unsafe or duplicate SAT control archive member')
            require(0 <= item.mode <= 0o777 and not item.mode&0o022,
                    'unsafe SAT control archive permissions')
            files[item.name]=(tar.extractfile(item).read(),item.mode)
    require(inventory(files)==receipt['inventory'],'SAT control published inventory changed')
    output.mkdir(parents=True)
    for role,(data,mode) in files.items():
        path=output/role
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(data)
        path.chmod(mode)
    spec=read(output/'registration/execution.json')
    require(sha256(spec)==receipt['execution_sha256'], 'SAT control invocation hash changed')
    result=audit(output/'execution',spec)
    require(result==receipt['result']==read(output/'independent-audit.json')==read(bundle/'AUDIT.json')
            and result['candidate_id']==receipt['candidate_id']
            and result['workload_id']==receipt['workload_id']
            and result['run_id']==receipt['run_id']
            and receipt['promotion_eligible'] is False and receipt['online_speedup'] is None,
            'SAT control transported mathematical/source replay differs')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    packing=sub.add_parser('publish')
    packing.add_argument('--registration',type=Path,required=True)
    packing.add_argument('--execution',type=Path,required=True)
    packing.add_argument('--audit',type=Path,required=True)
    packing.add_argument('--out',type=Path,required=True)
    transport=sub.add_parser('replay')
    transport.add_argument('--bundle',type=Path,required=True)
    transport.add_argument('--out',type=Path,required=True)
    transport.add_argument('--expected-execution-sha256',required=True)
    args=parser.parse_args()
    result=(publish(args.registration,args.execution,args.audit,args.out)
            if args.command=='publish' else replay(args.bundle,args.out,args.expected_execution_sha256))
    print(json.dumps({k:v for k,v in result.items() if k!='inventory'},sort_keys=True))
