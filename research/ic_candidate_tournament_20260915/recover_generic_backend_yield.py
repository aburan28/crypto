#!/usr/bin/env python3
"""Repair one frozen natural-yield auditor's float-only report-hash defect.

The 2026-09-29 campaign keeps its original evaluator and measured receipts.
Its natural-yield auditor mistakenly applies the candidate-identity hash to a
worker report containing floating-point timing diagnostics. This read-only
wrapper loads the sealed evaluator, replaces only that hash function with the
evaluator's ordinary JSON-object hash, and writes a separately labelled audit.
It is not a measurement retry or a replacement for `tournament.py verify`.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


FROZEN_AUDITOR_SHA256 = '627b5b84c2bbb9b9a9abb922a0ae2b6ed39ea70a9f212f087de7064e95637dc5'


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def recover(bundle):
    bundle = Path(bundle).resolve()
    round_dir = bundle / 'tournament'
    evaluator = round_dir / 'evaluator'
    auditor_path = evaluator / 'generic_backend_yield.py'
    contract_path = round_dir / 'contract.json'
    contract = json.loads(contract_path.read_text())
    if (contract.get('qualification_reference_schema') != 1
            or contract.get('seed') != 2026092901
            or contract.get('evaluator_sha256', {}).get('generic_backend_yield.py')
            != FROZEN_AUDITOR_SHA256
            or file_sha256(auditor_path) != FROZEN_AUDITOR_SHA256):
        raise ValueError('not the exact frozen 2026-09-29 auditor and contract')

    # Import all checker dependencies from the immutable archived evaluator.
    # `frozen_inputs` then checks every source and evaluator digest in the
    # contract; the raw receipts are checked by its separate `verify` command.
    if any(name in sys.modules for name in ('tournament', 'identity', 'oracle',
                                           'generic_stages', 'generic_phases')):
        raise RuntimeError('run this recovery in a fresh Python process')
    sys.path.insert(0, str(evaluator))
    spec = importlib.util.spec_from_file_location('frozen_generic_backend_yield', auditor_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from tournament import objhash
    if Path(sys.modules['tournament'].__file__).resolve() != (evaluator / 'tournament.py').resolve():
        raise RuntimeError('archived evaluator was not imported')

    if module.sha256.__module__ != 'identity':
        raise ValueError('frozen auditor has an unexpected report-hash helper')
    module.sha256 = objhash
    result = module.audit_campaign(round_dir)
    result['posthoc_audit_repair'] = dict(
        classification='post-hoc evaluator repair; no measured input or outcome changed',
        frozen_auditor_sha256=FROZEN_AUDITOR_SHA256,
        repair_script_sha256=file_sha256(__file__),
        contract_file_sha256=file_sha256(contract_path),
        change='report_sha256 uses sealed tournament.objhash, which accepts finite JSON timing floats; '
               'all natural-query, build, group, matrix and phase checks remain in the frozen auditor',
        prerequisite='run the unmodified archived tournament.py verify on all retained receipts')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = recover(args.bundle)
    with args.out.open('x') as target:
        json.dump(result, target, sort_keys=True, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(status=result['status'], observations=len(result['observations']),
                          rows=len(result['rows']), repaired=True)))


if __name__ == '__main__':
    main()
