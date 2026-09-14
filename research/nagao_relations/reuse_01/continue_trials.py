"""Preserve an interrupted run and finish only its missing frozen cells."""
import faulthandler
import hashlib
import json
import multiprocessing
import platform
import queue
import random
import subprocess
import time
import traceback
import trial

HERE, ROOT, base = trial.HERE, trial.ROOT, trial.base


def worker(output, item, variant, mode, budget, watchdog):
    faulthandler.dump_traceback_later(watchdog, exit=True)
    try:
        result = trial.cell(item, variant, mode, budget)
        output.put({'row': result})
    except BaseException:
        output.put({'error': traceback.format_exc()})
    finally:
        faulthandler.cancel_dump_traceback_later()


def main():
    original = (HERE / 'raw.jsonl').read_bytes()
    rows = list(map(json.loads, original.splitlines()))
    original_trials = [r for r in rows if r['kind'] == 'trial']
    contract = json.loads((HERE / 'contract.json').read_text())
    addendum = json.loads((HERE / 'runner_addendum.json').read_text())
    targets = json.loads((HERE / 'targets.json').read_text())
    jobs = []
    for item in targets:
        for mode in contract['modes']:
            variants = contract['variants'][:]
            random.Random(item['instance_sha256'] + mode).shuffle(variants)
            jobs += [(item, variant, mode) for variant in variants]
    expected_keys = [(i['instance_sha256'], v, m) for i, v, m in jobs]
    actual_keys = [(r['instance_sha256'], r['variant'], r['mode']) for r in original_trials]
    if len(original_trials) != 198 or actual_keys != expected_keys[:198]:
        raise ArithmeticError('interruption point differs from recorded incident')
    if expected_keys[198] != (addendum['affected_instance_sha256'], addendum['affected_variant'], addendum['affected_mode']):
        raise ArithmeticError('next cell differs from incident')
    paths = [HERE / p for p in ('continue_trials.py', 'runner_addendum.json', 'trial.py', 'reuse.py', 'contract.json', 'targets.json')]
    provenance = {'kind': 'continuation_provenance',
                  'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'original_raw_sha256': hashlib.sha256(original).hexdigest(),
                  'sha256': {str(p.relative_to(ROOT)): base.digest(p) for p in paths},
                  'python': platform.python_version(), 'hardware': platform.uname()._asdict(),
                  'start_job': 198, 'pending_jobs': 90, 'worker_threads': 1}
    instances = {r['instance_sha256'] for r in rows if r['kind'] == 'instance'}
    oracles, checkers = {}, {}
    ctx = multiprocessing.get_context('spawn')
    with (HERE / 'continuation_raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record(provenance)
        record({'kind': 'runner_incident', **addendum})
        for item, variant, mode in jobs[198:]:
            digest = item['instance_sha256']
            space_key = item['n'], tuple(item['basis']), tuple(item['coefficients'])
            if space_key not in oracles:
                f, c, v = base.old.context(item, False)
                oracles[space_key] = base.old.Oracle(f, c, v)
            oracle = oracles[space_key]
            if digest not in checkers:
                checkers[digest] = base.Checker(item)
            if digest not in instances:
                expected = oracle.expected(tuple(oracle.f.fromCoords(x) for x in item['target']))
                record({'kind': 'instance', **item, 'expected': [list(x) for x in sorted(expected)]})
                instances.add(digest)
            channel = ctx.Queue()
            process = ctx.Process(target=worker, args=(channel, item, variant, mode,
                                  contract['cold_budget_seconds'], addendum['hard_watchdog_seconds']))
            start = time.perf_counter()
            process.start()
            try:
                result = channel.get(timeout=addendum['hard_watchdog_seconds'] + 1)
                if 'error' in result:
                    raise ArithmeticError(result['error'])
                row = result['row']
            except queue.Empty:
                row = {'variant': variant, 'mode': mode, 'status': 'watchdog-timeout',
                       'within_budget': False, 'budget_seconds': contract['cold_budget_seconds'],
                       'all_phase_seconds': time.perf_counter() - start, 'first_verified_seconds': None,
                       'phase_seconds': None, 'solutions': [], 'certificates': [], 'rejection_certificates': [],
                       'verified_unique_relations': 0, 'counters': None, 'common_operation_speedup': None}
            finally:
                if process.is_alive():
                    process.join(timeout=.5)
                if process.is_alive():
                    process.terminate()
                process.join()
                channel.close()
            controller_seconds = time.perf_counter() - start
            t = time.perf_counter()
            faulthandler.dump_traceback_later(20)
            count = base.check(row, item, oracle, checkers[digest])
            faulthandler.cancel_dump_traceback_later()
            affected = (digest, variant, mode) == expected_keys[198]
            record({'kind': 'trial', **item, **row, 'expected_count': count,
                    'runner': 'isolated-continuation', 'controller_seconds': controller_seconds,
                    'worker_start_and_transfer_seconds': max(0., controller_seconds - row['all_phase_seconds']),
                    'harness_audit_seconds': time.perf_counter() - t,
                    'prior_unresolved_runner_incident': affected})
            print(item['n'], item['d'], item['cohort'], mode, variant, row['status'],
                  len(row['solutions']), round(row['all_phase_seconds'], 3), flush=True)


if __name__ == '__main__':
    main()
