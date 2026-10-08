"""Read-only historical registration replay using validated frozen sources."""
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import tempfile

from identity import sha256
from oracle import require

HERE = Path(__file__).resolve().parent
SNAPSHOT = HERE/'goal_20260924/static-sat-runtime'
SOLVER_BUILD_BUNDLE = Path(
    'research/sat_factor_base_review_20260908/continuation-05-sota-gates/'
    'stage-20-phase-b-terminal-evidence-successor-04-20260910')


def verified_sources(snapshot=SNAPSHOT):
    snapshot = Path(snapshot)
    receipt = json.loads((snapshot/'receipt.json').read_text())
    data = (snapshot/'runtime.tar.gz').read_bytes()
    require(hashlib.sha256(data).hexdigest() == receipt['archive_sha256'],
            'frozen SAT runtime archive changed')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            name = PurePosixPath(item.name)
            require(item.isfile() and not name.is_absolute()
                    and '..' not in name.parts and item.name not in files,
                    'unsafe or duplicate frozen SAT source')
            files[item.name] = tar.extractfile(item).read()
    require({name: hashlib.sha256(value).hexdigest() for name, value in files.items()}
            == receipt['files'], 'frozen SAT source inventory changed')
    for version, registration in receipt['registrations'].items():
        source = json.loads((HERE/'goal_20260924'/registration['directory']
                             /'source-manifest.json').read_text())
        expected = {item['role']: item['sha256'] for item in source['components']
                    if item['role'].endswith('.py')}
        require(registration['source_manifest_sha256'] == sha256(source)
                and registration['files'] == expected
                and all(receipt['files'].get(name) == value
                        for name, value in expected.items()),
                'frozen '+version+' runtime differs from immutable registration')
        require(registration['complete_preexecution_python_manifest'] is False
                and registration['supplemental_files']
                and all(name not in expected and receipt['files'].get(name) == value
                        for name, value in registration['supplemental_files'].items()),
                'historical SAT supplemental source limitation changed')
    return files, receipt


def replay(version):
    require(version in ('v1', 'v2'), 'unknown frozen SAT registration')
    files, receipt = verified_sources()
    suffix = '' if version == 'v1' else '_v2'
    script = '''import importlib,json,sys
from identity import sha256
from oracle import require
from tournament import read
suffix=sys.argv[1]
r=importlib.import_module('static_sat_registration'+suffix)
m=importlib.import_module('run_static_sat_full'+suffix)
panel=read(m.PANEL)
source,method,candidate,workload=r.identities(panel)
require(source==read(r.REGISTRATION/'source-manifest.json')
        and method==read(r.REGISTRATION/'method.json')
        and candidate==read(r.REGISTRATION/'candidate.json')
        and workload==read(r.REGISTRATION/'workload.json'),
        'frozen SAT registration failed source/method/identity replay')
_,curve,_,target,_,matrix,seal=m.admit(panel,require_local_binary=False)
counter,replayed=(m.target_from_seed(curve,2026092938) if not suffix else
                  m.target_from_seed(curve,2026092948,'ic-paired-target-v1'))
require(target==replayed and seal['candidate_id']==candidate['candidate_id'],
        'frozen SAT target or candidate changed')
print(json.dumps(dict(source_manifest_sha256=sha256(source),
    candidate_id=candidate['candidate_id'],workload_id=workload['workload_id'],
    target=list(target),target_counter=counter,columns=len(matrix.columns),
    actual_usable_points=candidate['record']['factor_base']['inventory']['usable_point_count']),
    sort_keys=True))
'''
    with tempfile.TemporaryDirectory(prefix='frozen-sat-runtime-') as temporary:
        root = Path(temporary)
        for name, data in files.items():
            target = root/name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        frozen_here = root/'research/ic_candidate_tournament_20260915'
        # Shared inputs are immutable registered data; executed Python comes
        # entirely from the hash-checked snapshot above.
        for name in ('goal_20260924', 'evidence'):
            (frozen_here/name).symlink_to(HERE/name, target_is_directory=True)
        bundle = root/SOLVER_BUILD_BUNDLE
        bundle.parent.mkdir(parents=True, exist_ok=True)
        bundle.symlink_to(HERE.parents[1]/SOLVER_BUILD_BUNDLE,
                          target_is_directory=True)
        environment = {key: value for key, value in os.environ.items()
                       if key not in ('PYTHONPATH', 'PYTHONSTARTUP')}
        command = [sys.executable, '-c',
                   'import sys;sys.path.insert(0,sys.argv.pop(1));'+script,
                   str(frozen_here), suffix]
        process = subprocess.run(command, cwd=root, env=environment,
                                 capture_output=True, text=True, timeout=60)
        require(process.returncode == 0,
                'frozen SAT registration replay failed: '+process.stderr[-6000:])
        result = json.loads(process.stdout)
    require(result['source_manifest_sha256']
            == receipt['registrations'][version]['source_manifest_sha256'],
            'frozen SAT replay used a different source manifest')
    result['complete_preexecution_python_manifest'] = False
    result['promotion_eligible'] = False
    return result
