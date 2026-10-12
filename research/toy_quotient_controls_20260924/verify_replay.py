"""Audit frozen evidence and optionally compare a fresh complete replay."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from run import summarize, digest

HERE=Path(__file__).resolve().parent


def read(folder):
    data=json.loads(gzip.decompress((folder/'raw.json.gz').read_bytes()))
    summary=json.loads((folder/'summary.json').read_text())
    if summarize(data)!=summary:
        raise AssertionError('frozen summary is inconsistent')
    if data['contract_sha256']!=digest(json.loads((HERE/'contract.json').read_text())):
        raise AssertionError('contract hash changed')
    for path,expected in data['source_sha256'].items():
        if hashlib.sha256((HERE.parent/path).read_bytes()).hexdigest()!=expected:
            raise AssertionError(f'source hash changed: {path}')
    if len(data['cases'])!=1920 or summary['verified']!=1920:
        raise AssertionError('missing or failed cases')
    if summary['independent_audits']!=240 or summary['independent_verified']!=240:
        raise AssertionError('independent Buchberger audit incomplete')
    for c in data['cases']:
        if len(set(c['repetition_sha256']))!=1:
            raise AssertionError('repetition hashes differ')
        if c['verified_solution_sha256']!=c['oracle_solution_sha256']:
            raise AssertionError('point oracle mismatch')
        if set(c['phase_opcodes'])!={'setup','encoding','solve_extract','reconstruct','verify'} or any(v<=0 for v in c['phase_opcodes'].values()):
            raise AssertionError('missing phase counts')
        if sum(c['phase_opcodes'].values())!=c['total_opcodes']:
            raise AssertionError('nonexclusive phase sum')
        del c['uninstrumented_phase_seconds']
    # Counts are a CPython 3.12 metric; patch metadata is not an algebraic result.
    if not data['python'].startswith('3.12.') or data['python_implementation']!='CPython':
        raise AssertionError('frozen opcode comparison requires CPython 3.12')
    del data['python']
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('fresh',type=Path,nargs='?')
    args=parser.parse_args()
    frozen=read(HERE/'results/run-002')
    if args.fresh and read(args.fresh)!=frozen:
        raise SystemExit('deterministic evidence differs')
    print('Evidence verified'+('; exact full replay matches' if args.fresh else ''))
