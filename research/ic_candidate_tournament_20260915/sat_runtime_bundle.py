"""Preexecution Python source snapshots for future complete SAT registrations.

The bounded surface includes every top-level tournament module, the entire
producer package, and the three Python child scripts used by the SAT runner.
An execution-time module gate rejects repository imports outside that surface.
Historical v1/v2 registrations must not be rewritten to use this manifest.
"""
import gzip
import hashlib
import importlib.util
import io
from pathlib import Path, PurePosixPath
import sys
import sysconfig
import tarfile

from identity import sha256, write_immutable
from oracle import require

DIRECTORY = Path('research/ic_candidate_tournament_20260915')
PACKAGES = ('producer',)
CHILD_SCRIPTS = (
    'scripts/process_meter.py',
    'scripts/run_koblitz_pdp_matrix.py',
    'scripts/verify_stage20_phase_b_terminal_evidence.py',
)


def source_files(repository):
    """Read a deliberately complete local surface instead of guessing imports."""
    repository = Path(repository).resolve()
    directory = repository/DIRECTORY
    require(directory.is_dir(), 'SAT runtime module directory missing')
    paths = set(directory.glob('*.py'))
    for name in PACKAGES:
        package = directory/name
        require(package.is_dir() and not package.is_symlink(),
                'SAT runtime package directory missing or symlinked')
        require(not any(path.is_symlink() for path in package.rglob('*')),
                'symlinked SAT runtime package member')
        paths.update(package.rglob('*.py'))
    paths.update(repository/role for role in CHILD_SCRIPTS)
    files = {}
    for path in sorted(paths):
        require(path.is_file() and not path.is_symlink()
                and path.resolve().is_relative_to(repository),
                'SAT runtime source missing, symlinked or outside repository')
        files[path.relative_to(repository).as_posix()] = path.read_bytes()
    require(files, 'SAT runtime surface is empty')
    return files


def manifest_for(files):
    return dict(schema_version=2,
                scope='complete bounded Python surface with loaded-module gate',
                module_directory=DIRECTORY.as_posix(),
                package_directories=list(PACKAGES),
                child_script_roles=list(CHILD_SCRIPTS),
                components=[dict(role=name, bytes=len(data),
                                 sha256=hashlib.sha256(data).hexdigest())
                            for name, data in sorted(files.items())])


def source_manifest(repository):
    return manifest_for(source_files(repository))


def freeze(repository, output):
    """Create a new source snapshot before execution; never replace a receipt."""
    output = Path(output)
    require(not output.exists(), 'SAT runtime snapshot output already exists')
    files = source_files(repository)
    manifest = manifest_for(files)
    output.mkdir(parents=True)
    archive = output/'runtime.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(
            fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for name, data in sorted(files.items()):
                item = tarfile.TarInfo(name)
                item.size, item.mode, item.mtime = len(data), 0o644, 0
                tar.addfile(item, io.BytesIO(data))
    write_immutable(output/'manifest.json', manifest)
    require(source_manifest(repository) == manifest,
            'SAT runtime source changed during snapshot creation')
    seal = dict(schema_version=1, manifest_sha256=sha256(manifest),
                archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                registration_stage='before-execution')
    write_immutable(output/'seal.json', seal)
    return manifest, seal


def verified_files(snapshot, expected_manifest, expected_seal):
    """Verify against the candidate's sealed values, independently of live code."""
    import json
    snapshot = Path(snapshot)
    require(json.loads((snapshot/'manifest.json').read_text()) == expected_manifest
            and json.loads((snapshot/'seal.json').read_text()) == expected_seal
            and sha256(expected_manifest) == expected_seal['manifest_sha256']
            and expected_manifest['schema_version'] == 2
            and expected_manifest['module_directory'] == DIRECTORY.as_posix()
            and expected_manifest['package_directories'] == list(PACKAGES)
            and expected_manifest['child_script_roles'] == list(CHILD_SCRIPTS)
            and expected_seal['registration_stage'] == 'before-execution',
            'SAT runtime snapshot differs from candidate registration')
    data = (snapshot/'runtime.tar.gz').read_bytes()
    require(hashlib.sha256(data).hexdigest() == expected_seal['archive_sha256'],
            'SAT runtime snapshot archive changed')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            name = PurePosixPath(item.name)
            require(item.isfile() and not name.is_absolute()
                    and '..' not in name.parts and item.name not in files,
                    'unsafe or duplicate SAT runtime archive entry')
            files[item.name] = tar.extractfile(item).read()
    require(manifest_for(files) == expected_manifest,
            'SAT runtime archive source inventory changed')
    return files


def extract(snapshot, output, expected_manifest, expected_seal):
    files = verified_files(snapshot, expected_manifest, expected_seal)
    output = Path(output)
    require(not output.exists(), 'SAT runtime extraction output already exists')
    output.mkdir(parents=True)
    for name, data in files.items():
        path = output/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(0o444)
    return output


def check_loaded_modules(repository, manifest, modules=None):
    """Reject an unbound local import before a result can be admitted.

    Call before measurement and again at termination to include lazy imports.
    External Python imports are limited to this interpreter's standard library;
    live-checkout fallbacks and site-packages require separate future admission.
    The snapshot launcher must also keep the source tree read-only during the
    run and set PYTHONDONTWRITEBYTECODE=1 in its fresh extraction directory.
    """
    repository = Path(repository).resolve()
    stdlib = Path(sysconfig.get_path('stdlib')).resolve()
    site_packages = {Path(sysconfig.get_path(key)).resolve()
                     for key in ('purelib', 'platlib')}
    expected = {item['role']: item for item in manifest['components']}
    require(len(expected) == len(manifest['components']),
            'duplicate SAT runtime source roles')
    loaded = {}
    for name, module in (sys.modules if modules is None else modules).items():
        origin = getattr(module, '__file__', None)
        if not origin:
            continue
        path = Path(origin).resolve()
        if not path.is_relative_to(repository):
            require(path.is_relative_to(stdlib)
                    and not any(path.is_relative_to(site) for site in site_packages),
                    'unregistered external SAT Python import: '+name)
            continue
        if path.suffix == '.pyc':
            try:
                path = Path(importlib.util.source_from_cache(str(path)))
            except ValueError:
                path = path.with_suffix('.py')
        role = path.relative_to(repository).as_posix()
        require(role in expected, 'unregistered local SAT runtime import: '+role)
        data = path.read_bytes()
        require(len(data) == expected[role]['bytes']
                and hashlib.sha256(data).hexdigest() == expected[role]['sha256'],
                'loaded SAT runtime source changed: '+role)
        loaded[name] = role
    return loaded
