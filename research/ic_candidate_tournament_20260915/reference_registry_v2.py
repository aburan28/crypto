"""Write exact v2 reference declarations from the accepted retained archive.

This transports declarations only. It never launches a worker, generates a
target, registers a candidate panel, or consumes an improvement attempt.
"""
import argparse
import copy
import json
from pathlib import Path

import campaign_rules_v2 as rules
from generic_build import digest, verify_build_record
from generic_reference_qualification import prepared_sources
from identity import sha256
from oracle import require


def declarations(bundle):
    bundle = Path(bundle).resolve()
    prepared = prepared_sources(bundle / 'controls')
    build_dir = bundle.parent / 'ic-generic-qualification-build'
    source = json.loads((build_dir / 'source-manifest.json').read_text())
    build = json.loads((build_dir / 'build-record.json').read_text())
    verify_build_record(build, source)
    require(sha256(source) == rules.GENERIC_SOURCE and build['build_sha256'] == rules.GENERIC_BUILD
            and digest(build_dir / 'worker') == build['worker_sha256'],
            'changed accepted generic reference build')
    report = json.loads((bundle / 'tournament/qualification.json').read_text())
    observer = json.loads((bundle / 'observer/summary.json').read_text())
    expected = rules.expected_binding()['bindings']
    bound = []
    registry = []
    for name, spec in expected.items():
        arm = dict(id=name, config=copy.deepcopy(spec['configuration']),
                   source_manifest_sha256=spec['source_manifest_sha256'])
        declaration = dict(id=name, config=copy.deepcopy(spec['configuration']))
        if name != 'incumbent':
            arm['kind'] = declaration['kind'] = 'ic-reference' if name == 'ic_online' else 'rho-reference'
        if spec['adapter']:
            arm.update(adapter=spec['adapter'], build_sha256=build['build_sha256'])
            declaration.update(adapter=spec['adapter'], generic_build=str(build_dir))
        else:
            label = 'both' if name == 'ic_online' else 'pairinv'
            actual = json.loads((prepared[label] / 'source-manifest.json').read_text())
            require(sha256(actual) == spec['source_manifest_sha256'], 'changed prepared reference source')
            declaration['source_root'] = str(prepared[label] / 'source')
        bound.append(arm)
        if name != 'incumbent':
            registry.append(declaration)
    rules.qualified_binding(report, observer, bound[0], bound[1:])
    rules.validate_reference_declarations(registry)
    return registry, dict(incumbent_source=str(prepared['pairinv'] / 'source'),
        qualified_report=str(bundle / 'tournament/qualification.json'),
        qualified_observer=str(bundle / 'observer/summary.json'),
        workers_executed=0, targets_generated=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True,
                        help='Restored ic-generic-reference-qualification directory')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    registry, inputs = declarations(args.bundle)
    with args.out.open('x') as stream:
        json.dump(registry, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(status='DECLARED', reference_registry=str(args.out.resolve()), **inputs)))


if __name__ == '__main__':
    main()
