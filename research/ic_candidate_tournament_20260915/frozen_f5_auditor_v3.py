"""Replay mathematical admission from the retained postexecution auditor source.

This source context is independent analysis. It is not a retroactive change
to the preregistered method, controller, Python runtime or native invocation.
"""
import argparse
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--execution', type=Path, required=True)
    parser.add_argument('--spec', type=Path, required=True)
    parser.add_argument('--expected-execution-sha256', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--gate-out', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    sys.path.insert(0, str(root/'research/ic_candidate_tournament_20260915'))
    from audit_f5_runtime_v3 import audit
    from identity import sha256, write_immutable
    from oracle import require
    from sat_runtime_bundle import check_loaded_modules, source_manifest
    manifest = json.loads(args.manifest.read_text())
    require(source_manifest(root) == manifest, 'retained F5 auditor source differs')
    before = check_loaded_modules(root, manifest)
    spec = json.loads(args.spec.read_text())
    require(sha256(spec) == args.expected_execution_sha256,
            'retained F5 auditor invocation hash differs')
    result = audit(args.execution, spec)
    after = check_loaded_modules(root, manifest)
    require(source_manifest(root) == manifest, 'F5 auditor source changed during replay')
    write_immutable(args.out, result)
    write_immutable(args.gate_out, dict(schema_version=1,
        source_manifest_sha256=sha256(manifest), before=before, after=after,
        python_version=list(sys.version_info[:3]), isolated=True, site_disabled=True,
        stage='postexecution-independent-auditor-context', native_execution=False))


if __name__ == '__main__':
    main()
