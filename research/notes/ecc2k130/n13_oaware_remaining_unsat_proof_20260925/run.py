#!/usr/bin/env python3
"""Fail-preserving remaining 22 fixed n13 CaDiCaL-to-DRAT-trim proof pipeline."""
from __future__ import annotations

import argparse
import contextlib
import gzip
import json
import platform
import resource
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from common import HERE, file_hashes, freeze, query_bytes, sha, sha_bytes, checker_selftest


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def peak_child_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return value if sys.platform == 'darwin' else value * 1024


def archive_size(root: Path) -> int:
    return sum(path.stat().st_size for path in root.rglob('*') if path.is_file())


@contextlib.contextmanager
def retained_scratch(root: Path):
    """Keep a failed child's exact query and available proof for audit."""
    path = root / 'scratch'
    path.mkdir()
    try:
        yield path
    except Exception:
        raise
    else:
        shutil.rmtree(path)


def child(name: str, argv: list[str], cap: int, output: Path, receipt: dict) -> dict:
    entry = {'name': name, 'argv': argv, 'started_utc': utc(),
             'timeout_seconds': cap + 15}
    receipt['commands'].append(entry)
    start = time.perf_counter()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    stdout = output / f'{name}.stdout.txt'
    stderr = output / f'{name}.stderr.txt'
    with stdout.open('wb') as out, stderr.open('wb') as err:
        try:
            process = subprocess.run(argv, stdout=out, stderr=err,
                                     timeout=cap + 15, check=False)
            entry['exit_code'] = process.returncode
        except subprocess.TimeoutExpired:
            entry['exit_code'] = 'TIMEOUT'
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    entry.update({'finished_utc': utc(), 'wall_seconds': time.perf_counter() - start,
                  'child_user_cpu_seconds': usage_after.ru_utime - usage_before.ru_utime,
                  'child_system_cpu_seconds': usage_after.ru_stime - usage_before.ru_stime,
                  'peak_child_rss_bytes_upper': peak_child_rss(),
                  'stdout_sha256': sha(stdout), 'stderr_sha256': sha(stderr)})
    return entry


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = {'status': 'started', 'started_utc': utc(), 'commands': [], 'cases': [],
               'python': sys.version, 'platform': platform.platform()}
    try:
        data, frozen, base, schema = freeze()
        assert data['parent_proof_merge_commit'] is not None, 'dependent #790 merge gate is not released'
        subprocess.run(['git', 'merge-base', '--is-ancestor', data['parent_proof_head'],
                        data['parent_proof_merge_commit']], cwd=HERE.parents[3], check=True, timeout=10)
        subprocess.run(['git', 'merge-base', '--is-ancestor', data['parent_proof_merge_commit'],
                        'origin/main'], cwd=HERE.parents[3], check=True, timeout=10)
        receipt['freeze_sha256'] = sha(HERE / 'FROZEN.json')
        receipt['source_sha256'] = frozen
        receipt['parent_main_merge'] = data['parent_main_merge']
        receipt['parent_proof_head'] = data['parent_proof_head']
        receipt['parent_proof_merge_commit'] = data['parent_proof_merge_commit']
        receipt['parent_proof_receipt_sha256'] = data['parent_proof_receipt_sha256']
        receipt['checker_upstream_commit'] = data['checker_upstream_commit']
        solver, checker = Path(data['solver_path']), Path(data['checker_build'][-1])
        assert solver.is_file() and sha(solver) == data['solver_binary_sha256']
        assert subprocess.check_output([str(solver), '--version'], text=True).strip() == data['solver_version']
        assert checker.is_file() and sha(checker) == data['checker_mac_binary_sha256']
        receipt['solver_binary_sha256'] = sha(solver)
        receipt['checker_binary_sha256'] = sha(checker)
        receipt['checker_build_recipe'] = data['checker_build']
        receipt['checker_compiler'] = subprocess.check_output(['cc', '--version'], text=True).splitlines()[0]
        control = args.out / 'control'
        receipt['mutated_proof_control'] = checker_selftest(checker, control)
        assert archive_size(args.out) <= data['caps']['archive_bytes']
        proofs = args.out / 'proofs'
        proofs.mkdir()
        with retained_scratch(args.out) as temp:
            first_proof = None
            for case in data['cases']:
                ident = case['id']
                result = {'id': ident, 'literal': case['literal'],
                          'point_oracle_positive': False,
                          'prior_query_sha256': case['prior_query_sha256'],
                          'status': 'started'}
                receipt['cases'].append(result)
                query = temp / f'{ident}.cnf'
                proof = temp / f'{ident}.drat'
                raw_query = query_bytes(base, case['literal'], schema['variables'], schema['clauses'])
                assert sha_bytes(raw_query) == case['prior_query_sha256']
                query.write_bytes(raw_query)
                result.update({'query_sha256': sha(query), 'query_bytes': query.stat().st_size})
                solver_entry = child(f'{ident}-solver', [str(solver), *data['solver_options'],
                                                          str(query), str(proof)],
                                     data['caps']['solver_wall_seconds'], args.out, receipt)
                result['solver_command'] = solver_entry['name']
                assert solver_entry['exit_code'] == 20
                assert solver_entry['wall_seconds'] <= data['caps']['solver_wall_seconds']
                assert solver_entry['peak_child_rss_bytes_upper'] <= data['caps']['rss_bytes']
                solver_text = (args.out / f'{ident}-solver.stdout.txt').read_text(errors='replace')
                assert 's UNSATISFIABLE' in solver_text
                assert proof.is_file() and 0 < proof.stat().st_size <= data['caps']['proof_bytes']
                raw_proof = proof.read_bytes()
                compressed = proofs / f'{ident}.drat.gz'
                compressed.write_bytes(gzip.compress(raw_proof, mtime=0))
                result.update({'raw_proof_bytes': len(raw_proof),
                               'raw_proof_sha256': sha_bytes(raw_proof),
                               'gzip_proof_bytes': compressed.stat().st_size,
                               'gzip_proof_sha256': sha(compressed)})
                assert archive_size(args.out) <= data['caps']['archive_bytes']
                assert archive_size(temp) <= data['caps']['temporary_bytes']
                check_entry = child(f'{ident}-checker', [str(checker), str(query), str(proof)],
                                    data['caps']['checker_wall_seconds'], args.out, receipt)
                result['checker_command'] = check_entry['name']
                assert check_entry['exit_code'] == 0
                assert check_entry['wall_seconds'] <= data['caps']['checker_wall_seconds']
                assert check_entry['peak_child_rss_bytes_upper'] <= data['caps']['rss_bytes']
                check_text = (args.out / f'{ident}-checker.stdout.txt').read_text(errors='replace')
                assert 's VERIFIED' in check_text and 's NOT VERIFIED' not in check_text
                result['status'] = 'PASS'
                if ident == 'Q0T1':
                    first_proof = proof
                else:
                    proof.unlink()
                query.unlink()
            assert first_proof is not None
            positive = data['positive_control']
            sat_query = temp / 'Q0T3.cnf'
            raw_query = query_bytes(base, positive['literal'], schema['variables'], schema['clauses'])
            assert sha_bytes(raw_query) == positive['prior_query_sha256']
            sat_query.write_bytes(raw_query)
            assert archive_size(temp) <= data['caps']['temporary_bytes']
            cross = child('Q0T1-proof-on-Q0T3', [str(checker), str(sat_query), str(first_proof)],
                          data['caps']['checker_wall_seconds'], args.out, receipt)
            receipt['wrong_input_control'] = {'proof_id': 'Q0T1', 'sat_query_id': 'Q0T3',
                                              'sat_query_sha256': sha(sat_query),
                                              'command': cross['name']}
            cross_text = (args.out / 'Q0T1-proof-on-Q0T3.stdout.txt').read_text(errors='replace')
            assert cross['exit_code'] != 0 and 's NOT VERIFIED' in cross_text
            assert cross['wall_seconds'] <= data['caps']['checker_wall_seconds']
            assert cross['peak_child_rss_bytes_upper'] <= data['caps']['rss_bytes']
        assert archive_size(args.out) <= data['caps']['archive_bytes']
        receipt['status'] = 'success'
    except Exception as error:
        receipt['status'] = 'failed'
        receipt['failure'] = repr(error)
    finally:
        receipt['finished_utc'] = utc()
        receipt['archive_bytes'] = archive_size(args.out)
        receipt['files_sha256'] = file_hashes(args.out)
        (args.out / 'receipt.json').write_text(
            json.dumps(receipt, sort_keys=True, separators=(',', ':')) + '\n')
    return int(receipt['status'] != 'success')


if __name__ == '__main__':
    raise SystemExit(main())
