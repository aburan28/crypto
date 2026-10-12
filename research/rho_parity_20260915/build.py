#!/usr/bin/env python3
"""Isolate package outputs by revision; only dependency artifacts may be reused."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

ap = argparse.ArgumentParser()
ap.add_argument('--manifest', type=Path, required=True)
ap.add_argument('--target', type=Path, required=True)
ap.add_argument('--dependency-seed', type=Path)
ap.add_argument('--output', type=Path, required=True)
ap.add_argument('--cargo', type=Path, required=True)
ap.add_argument('--log', type=Path, required=True)
ap.add_argument('--cost-only', action='store_true')
args = ap.parse_args()
assert not args.target.exists(), 'new isolated target required'
if args.dependency_seed:
    shutil.copytree(args.dependency_seed, args.target)
else:
    args.target.mkdir(parents=True)
# Package-relative dep-info can survive a copied target. Remove all crypto
# fingerprints explicitly, before Cargo gets a chance to call them fresh.
fingerprints = args.target/'release/.fingerprint'
for p in fingerprints.glob('crypto-*'):
    shutil.rmtree(p)
for pattern in ['release/deps/libcrypto_lib-*', 'release/deps/crypto_lib-*',
                'release/libcrypto_lib.*', 'release/examples/rho_parity_cost*',
                'release/examples/weil_factor_composition*', 'release/examples/f4_linear_algebra_bench*']:
    for p in args.target.glob(pattern):
        if p.is_file():
            p.unlink()
env = os.environ.copy()
env['CARGO_TARGET_DIR'] = str(args.target.resolve())
common = [str(args.cargo.absolute()), '--manifest-path', str(args.manifest.resolve())]
names = ['rho_parity_cost'] if args.cost_only else ['weil_factor_composition','f4_linear_algebra_bench','rho_parity_cost']
args.log.parent.mkdir(parents=True, exist_ok=True)
with args.log.open('w') as log:
    # This is mandatory even if a dependency-seed directory was copied.
    subprocess.run([common[0],'clean',*common[1:],'--release','-p','crypto'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run([common[0],'build',*common[1:],'--locked','--release','--features','redis-cache',
                    *[x for name in names for x in ['--example',name]]],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
args.output.mkdir(parents=True, exist_ok=True)
for name in names:
    shutil.copy2(args.target/'release/examples'/name, args.output/name)
print('Saved isolated binaries:', ', '.join(names))
