"""Replay transported closed Job schema evidence without native execution."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

import f5_job_schema_control
import generic_stages
import identity
import oracle
from identity import sha256
from oracle import require
from replay_paired_n17_evidence import retained_files

PREFIX = 'analysis-src/research/ic_candidate_tournament_20260915/'
REGISTRATION = PREFIX+'goal_20260924/f5-source-bound-runtime-v2/interface-control/'


def source_archive(data, expected):
    observed = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            require(item.isfile() and item.name not in observed, 'non-file or duplicate source member')
            observed[item.name] = hashlib.sha256(tar.extractfile(item).read()).hexdigest()
    require(observed == expected, 'retained source/dependency archive differs from compiled manifest')


def replay_files(files):
    record = json.loads(files['preexecution.json'])
    registration = json.loads(files['native-registration.json'])
    result = json.loads(files['RESULT.json'])
    protocol = json.loads(files[REGISTRATION+'protocol.json'])
    source = json.loads(files['source-manifest.json'])
    require(record['protocol_sha256'] == hashlib.sha256(files[REGISTRATION+'protocol.json']).hexdigest()
            and record['protocol_document_sha256'] == hashlib.sha256(files[REGISTRATION+'PROTOCOL.md']).hexdigest(),
            'registered protocol/document differs')
    require(record['limits'] == dict(build_seconds=protocol['build_seconds'], native_seconds=protocol['native_seconds'])
            and registration['native_seconds'] == protocol['native_seconds'], 'registered limit differs')
    require(hashlib.sha256(files['preexecution.json']).hexdigest() == result['preexecution_sha256']
            == registration['preexecution_sha256'], 'preexecution seal differs')
    require(hashlib.sha256(files['native-registration.json']).hexdigest() == result['native_registration_sha256'],
            'native registration seal differs')
    require(hashlib.sha256(files['diagnostic']).hexdigest() == registration['binary_sha256'],
            'native executable differs from preexecution registration')
    require(record['analysis_sources'] == result['analysis_sources'], 'executed analysis inventory differs')
    for name, expected in record['analysis_sources'].items():
        require('/research/' in name, 'unsupported retained analysis path')
        role = 'analysis-src/research/'+name.split('/research/', 1)[1]
        require(hashlib.sha256(files[role]).hexdigest() == expected, 'executed analysis source differs')
    for module in (f5_job_schema_control, generic_stages, identity, oracle):
        require(hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
                == hashlib.sha256(files[PREFIX+Path(module.__file__).name]).hexdigest(),
                'replay mathematics source differs from executed checker')
    context = json.loads(files['independent-auditor-context.json'])
    require(context['stage'] == 'postexecution-source-retention'
            and context['preregistered_commit'] == '8bf6181bd0f39ae1296378730267de8d771ec18a'
            and context['generic_stages_sha256'] == hashlib.sha256(files[PREFIX+'generic_stages.py']).hexdigest(),
            'postexecution auditor context differs')
    base = json.loads(files['base-source-manifest.json'])
    kernel_name = 'examples/ic_tournament_worker.rs'
    original = files['kernel-original.rs']
    append = files[REGISTRATION+'kernel_append.rs']
    require(hashlib.sha256(original).hexdigest() == base['root_files'][kernel_name]
            and hashlib.sha256(append).hexdigest() == protocol['kernel_append_sha256'],
            'original kernel or diagnostic append differs')
    compiled_hash = hashlib.sha256(original+append).hexdigest()
    expected = json.loads(json.dumps(base))
    expected['root_files'][kernel_name] = compiled_hash
    expected['diagnostic_overlay'] = {
        'examples/f5_job_schema_control.rs': hashlib.sha256(files[REGISTRATION+'export.rs']).hexdigest(),
        kernel_name: dict(base_sha256=base['root_files'][kernel_name],
                         append_sha256=protocol['kernel_append_sha256'], compiled_sha256=compiled_hash)}
    require(source == expected, 'compiled source differs from exact append-only diagnostic overlay')
    require(sha256(base) == protocol['historical_source_manifest_sha256']
            == record['historical_source_manifest_sha256'] == result['historical_source_manifest_sha256'],
            'historical Rust source identity differs')
    require(sha256(source) == record['compiled_source_manifest_sha256']
            == result['compiled_source_manifest_sha256'], 'compiled source identity differs')
    source_archive(files['root-source.tar.gz'], dict(source['root_files'], **{
        'examples/f5_job_schema_control.rs': source['diagnostic_overlay']['examples/f5_job_schema_control.rs']}))
    dependencies = {item['package']+'-'+item['version']+'/'+name: expected
                    for item in base['dependencies'] for name, expected in item['files'].items()}
    source_archive(files['dependency-source.tar.gz'], dependencies)
    inputs = files[REGISTRATION+'inputs.json']
    require(hashlib.sha256(inputs).hexdigest() == record['input_sha256']
            == registration['inputs_sha256'] == protocol['input_sha256'], 'input identity differs')
    for phase in ('build', 'native'):
        process = json.loads(files[phase+'-process.json'])
        require(process['exit_code'] == 0 and process['timeout'] is False
                and process['owned_group_cleanup'] is True, 'missing successful terminal '+phase+' process')
        require(process['argv'] == (record['build_argv'] if phase == 'build' else registration['argv'])
                and process['watchdog_seconds'] == protocol[phase+'_seconds'], 'process argv/limit differs')
        for stream in ('stdout', 'stderr'):
            require(hashlib.sha256(files[phase+'.'+stream]).hexdigest() == process[stream+'_sha256'],
                    'raw process output differs')
    require(hashlib.sha256(files['native.stdout']).hexdigest() == result['raw_export_sha256'],
            'native export identity differs')
    derived = f5_job_schema_control.audit(json.loads(inputs), json.loads(files['native.stdout']))
    require(all(result.get(key) == value for key, value in derived.items()), 'retained mathematical receipt differs')
    return derived


def replay(bundle, expected_archive_sha256):
    receipt = json.loads((Path(bundle)/'receipt.json').read_text())
    require(receipt['archive_sha256'] == expected_archive_sha256,
            'Job schema evidence differs from externally retained archive seal')
    return replay_files(retained_files(bundle))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--expected-archive-sha256', required=True)
    args = parser.parse_args()
    print(json.dumps(replay(args.bundle, args.expected_archive_sha256), sort_keys=True))
