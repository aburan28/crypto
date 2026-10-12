"""Rebuild the pinned public WDSat source; requires GCC and libm."""
from pathlib import Path
import hashlib, json, subprocess

root=Path(__file__).resolve().parent
for item in json.loads((root/'source_manifest.json').read_text()):
    p=root/item['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'], f'source hash mismatch: {p}'
cmd=['gcc','-O3','-Wall']+[str(p) for p in sorted((root/'vendor/WDSat/src').glob('*.c'))]+['-lm','-o',str(root/'vendor/WDSat/wdsat_solver')]
out=root/'results';out.mkdir(exist_ok=True)
p=subprocess.run(cmd,capture_output=True,text=True)
(out/'build.stdout').write_text(p.stdout)
(out/'build.stderr').write_text(p.stderr)
print(json.dumps({'exit_code':p.returncode,'compiler_command':cmd,'diagnostics':'pilot/results/build.stderr'}))
raise SystemExit(p.returncode)
