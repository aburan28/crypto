"""Immutable non-Python inputs for the complete source-bound SAT pipeline."""
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import tarfile

from identity import sha256, write_immutable
from oracle import require


def manifest_for(files, executable_roles):
    require(set(executable_roles) <= set(files), 'missing executable SAT asset')
    return dict(schema_version=3, components=[
        dict(role=role, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
             executable=role in executable_roles)
        for role, data in sorted(files.items())])


def safe_role(role):
    require(type(role) is str, 'non-string SAT asset role')
    name = PurePosixPath(role)
    require(role and role != '.' and name.as_posix() == role
            and not name.is_absolute() and '..' not in name.parts
            and name.suffix not in ('.py', '.pyc'),
            'unsafe or Python-source SAT input role')


def freeze_assets(files, executable_roles, output):
    """Retain exact bytes once; source snapshots remain a separate contract."""
    output = Path(output)
    require(not output.exists(), 'SAT asset snapshot already exists')
    for role, data in files.items():
        safe_role(role)
        require(type(data) is bytes, 'SAT asset input must be exact bytes')
    manifest = manifest_for(files, executable_roles)
    output.mkdir(parents=True)
    archive = output/'assets.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(
            fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for role, data in sorted(files.items()):
                item = tarfile.TarInfo(role)
                item.size, item.mtime = len(data), 0
                item.mode = 0o555 if role in executable_roles else 0o444
                tar.addfile(item, io.BytesIO(data))
    seal = dict(schema_version=3, manifest_sha256=sha256(manifest),
                archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                registration_stage='before-execution')
    write_immutable(output/'manifest.json', manifest)
    write_immutable(output/'seal.json', seal)
    return manifest, seal


def verified_assets(snapshot, expected_manifest, expected_seal):
    snapshot = Path(snapshot)
    require(json.loads((snapshot/'manifest.json').read_text()) == expected_manifest
            and json.loads((snapshot/'seal.json').read_text()) == expected_seal
            and expected_manifest['schema_version'] == expected_seal['schema_version'] == 3
            and sha256(expected_manifest) == expected_seal['manifest_sha256']
            and expected_seal['registration_stage'] == 'before-execution',
            'SAT inputs differ from registered asset seal')
    archive = (snapshot/'assets.tar.gz').read_bytes()
    require(hashlib.sha256(archive).hexdigest() == expected_seal['archive_sha256'],
            'SAT asset archive changed')
    files, executable = {}, set()
    with tarfile.open(fileobj=io.BytesIO(archive), mode='r:gz') as tar:
        for item in tar:
            safe_role(item.name)
            require(item.isfile() and item.name not in files
                    and item.mode in (0o444, 0o555),
                    'unsafe or duplicate SAT asset member')
            files[item.name] = tar.extractfile(item).read()
            if item.mode == 0o555:
                executable.add(item.name)
    require(manifest_for(files, executable) == expected_manifest,
            'SAT asset inventory changed')
    return files


def extract_assets(snapshot, output, expected_manifest, expected_seal):
    files = verified_assets(snapshot, expected_manifest, expected_seal)
    output = Path(output)
    require(not output.exists(), 'SAT asset extraction already exists')
    output.mkdir(parents=True)
    executable = {item['role'] for item in expected_manifest['components']
                  if item['executable']}
    for role, data in files.items():
        path = output/role
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(0o555 if role in executable else 0o444)
    return output


def check_extracted_assets(root, expected_manifest):
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), 'SAT asset tree missing or symlinked')
    actual, executable = {}, set()
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'symlinked extracted SAT asset')
        if path.is_file():
            role = path.relative_to(root).as_posix()
            safe_role(role)
            actual[role] = path.read_bytes()
            if path.stat().st_mode & 0o111:
                executable.add(role)
    require(manifest_for(actual, executable) == expected_manifest,
            'extracted SAT assets differ from registered inputs')
    return actual
