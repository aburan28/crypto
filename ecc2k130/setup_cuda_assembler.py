"""Prepare the measured CUDA 13.3 front end with the pinned 13.4 assembler.

This optional local setup never changes a system CUDA installation. Existing
prefixes must already match; an incompatible directory is left untouched.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
import urllib.request

VERSION = '13.4.59'
URL = 'https://developer.download.nvidia.com/compute/cuda/redist/cuda_nvcc/linux-x86_64/cuda_nvcc-linux-x86_64-13.4.59-archive.tar.xz'
ARCHIVE_SHA256 = '0c08d1df80b5d0bd081778d392446ab50a6a54b047c0136b39c8dbfdf2cfed3f'
PTXAS_SHA256 = '9c2c084df7eb9be48f0502b6938ce2571679b2a78cad1c564b583bd7d9ea3c50'
BASE_PTXAS_SHA256 = 'afd8d1e1fa6e310f7faee44f6621e4c1315fb7fd6da7d4d87414358e12a651dc'
MEMBER = 'cuda_nvcc-linux-x86_64-13.4.59-archive/bin/ptxas'

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def manifest(prefix):
    return {str(p.relative_to(prefix)): digest(p) for p in prefix.rglob('*') if p.is_file()}

def verify(base_files, prefix):
    files = manifest(prefix)
    expected = dict(base_files, **{'bin/ptxas': PTXAS_SHA256})
    if files != expected:
        raise RuntimeError(f'{prefix} does not match the pinned assembler-only prefix; left unchanged')
    return files

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('build/cuda133/nvidia/cu13'))
    parser.add_argument('--prefix', type=Path, default=Path('build/assembler134-toolchain'))
    parser.add_argument('--archive', type=Path, help='Use an already downloaded, checksum-verified archive')
    parser.add_argument('--check', action='store_true', help='Verify an existing prefix without downloading or writing')
    args = parser.parse_args()
    base = args.base.resolve()
    prefix = args.prefix.absolute()
    base_files = manifest(base)
    if base_files.get('bin/ptxas') != BASE_PTXAS_SHA256 or not (base/'bin/nvcc').is_file():
        parser.error('The pinned CUDA 13.3.73 base is missing; run make setup-local-cuda133 first')
    if prefix.exists():
        files = verify(base_files, prefix)
        created = False
    elif args.check:
        parser.error('The assembler prefix does not exist')
    else:
        prefix.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='assembler134-', dir=prefix.parent) as temporary:
            staging = Path(temporary)
            archive = args.archive
            if archive is None:
                archive = staging/'cuda_nvcc.tar.xz'
                with urllib.request.urlopen(URL, timeout=60) as source, archive.open('wb') as target:
                    shutil.copyfileobj(source, target)
            if digest(archive) != ARCHIVE_SHA256:
                raise RuntimeError('CUDA assembler archive checksum mismatch')
            candidate = staging/'prefix'
            shutil.copytree(base, candidate, symlinks=True)
            with tarfile.open(archive) as package:
                member = package.getmember(MEMBER)
                if not member.isfile():
                    raise RuntimeError('Expected a regular assembler executable')
                with package.extractfile(member) as source, (candidate/'bin/ptxas').open('wb') as target:
                    shutil.copyfileobj(source, target)
            (candidate/'bin/ptxas').chmod(0o755)
            files = verify(base_files, candidate)
            os.rename(candidate, prefix)
        created = True
    print(json.dumps(dict(prefix=str(prefix),nvcc=str(prefix/'bin/nvcc'),assembler_version=VERSION,
                          assembler_sha256=files['bin/ptxas'],changed_files=['bin/ptxas'],created=created),indent=2))

if __name__ == '__main__':
    main()
