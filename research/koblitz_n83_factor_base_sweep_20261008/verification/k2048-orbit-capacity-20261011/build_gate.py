import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

root = Path('/Volumes/SSD990/crypto/worktrees/codex-n83-upstream-integration-20261010')
out = Path('/private/tmp/n83-k2048-orbit-capacity-20261011')
target = Path('/private/tmp/n83-k2048-capacity-linux-target')
registry = Path('/Users/adamburan/.cargo/registry')
image = 'sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31'
out.mkdir(exist_ok=False)
target.mkdir(exist_ok=False)
command = ['docker','run','--rm','--name','n83-k2048-build-20261011','--pull','never','--platform','linux/arm64','--network','none','--cpus','4','--memory','8192m','--memory-swap','8192m','--pids-limit','256','--cap-drop','ALL','--security-opt','no-new-privileges','--user',f'{os.getuid()}:{os.getgid()}','--tmpfs','/cargo:rw,mode=1777','--mount',f'type=bind,src={registry},dst=/cargo/registry,readonly','--mount',f'type=bind,src={root},dst=/src,readonly','--mount',f'type=bind,src={target},dst=/target','--env','CARGO_HOME=/cargo','--env','CARGO_TARGET_DIR=/target','--workdir','/src',image,'cargo','build','--release','--offline','--locked','--no-default-features','--example','koblitz_n83_factor_base_export']
start = time.monotonic()
status = 'UNKNOWN'
code = None
with (out/'build.stdout').open('wb') as stdout, (out/'build.stderr').open('wb') as stderr:
    try:
        p = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=600, check=False)
        code = p.returncode
        status = 'PASS_build' if code == 0 else 'PRODUCER_FAILURE_build'
    except subprocess.TimeoutExpired:
        subprocess.run(['docker','kill','n83-k2048-build-20261011'],capture_output=True,timeout=15,check=False)
        status = 'UNKNOWN_build_wall_cap'
wall = time.monotonic()-start
remaining = subprocess.run(['docker','ps','-aq','--filter','name=^n83-k2048-build-20261011$'],capture_output=True,text=True,timeout=15,check=False)
cleanup_verified = remaining.returncode == 0 and not remaining.stdout.strip()
if not cleanup_verified:
    status = 'PRODUCER_FAILURE_build_cleanup'
binary = target/'release/examples/koblitz_n83_factor_base_export'
h = hashlib.sha256(binary.read_bytes()).hexdigest() if status == 'PASS_build' and binary.exists() else None
receipt = {'schema':'n83.k2048-orbit-capacity-build/v1','status':status,'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'command':command,'image_id':image,'wall_elapsed_seconds_L0':wall,'wall_cap_seconds':600,'exit_code':code,'cleanup_verified':cleanup_verified,'binary':str(binary),'binary_sha256':h,'binary_bytes':binary.stat().st_size if h else None}
(out/'build.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:receipt[k] for k in ['status','wall_elapsed_seconds_L0','exit_code','binary_sha256']}))
