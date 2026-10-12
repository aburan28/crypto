"""Certify reusable preparation from two accepted, closed n17 IC controls.

This reads retained bytes and reconstructs ordinary relations and column logs.
It launches no native executable, generates no target, and awards no speedup.
ICP1 identifies mathematical preparation, not a complete IC candidate or run.
"""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile

from identity import (curve_record, factor_base_inventory, fields,
                      natural, sha256, write_immutable)
from oracle import Curve, require
from static_sat_matrix import RelationMatrix


HERE = Path(__file__).resolve().parent
POLICY = 'ordinary-relations-only-cofactor-sign-frobenius-v1'
STATE_SHA256 = 'edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107'
ACCEPTED = {
    'f5': dict(
        bundle='f5-source-bound-runtime-v2/results-20260930',
        archive_sha256='d62ff4b0e5848b7e36474d1a6af2be7fef1ee060d6cc05ae1e02c86c2957c1b2',
        archive_bytes=38314958,
        execution_sha256='9409cb6c40548816153f5b9978e8fa0637b5aaa2801472575e7cc55a9f61d56d',
        candidate_id='IC1N17Ckb1fb62PDP3f5RCsampleLAgaussTDpdpISO0h84b504d19844',
        candidate_record_sha256='84b504d19844b94096c94b4cf60b3fb279936ec357df5bb3bb10767823a50891',
        inputs_sha256='6678815b78b9414d81f94b595d3e10ab38f388511fa3a6831891f685e96f833d',
        merge_commit='07d7c63636965a47035c66e37f665d5611b6a7a0'),
    'sat': dict(
        bundle='static-sat-runtime-v3/full-development-20260929/results-20260930',
        archive_sha256='86e84fa1a31d152ef93a88837103ca578862269f333c89f44df26c472030ed84',
        archive_bytes=28871942,
        execution_sha256='4e6547c03cfe1fa54ae9a994c160e5d4afa987659dc72880c651885a0e9e8fa9',
        candidate_id='IC1N17Ckb1fb62PDP3satRCsampleLAgaussTDpdpISO0h0ef1e080d97e',
        candidate_record_sha256='0ef1e080d97ebec2cb3f079ece7ea7262fbb36c91c70fccb3fc25c6b3e71a9cf',
        inputs_sha256='af043f2d1a4a77a5840a73b081075e19ea7fe058355d58ff2e5ee701abc02ce0',
        merge_commit='9f1980d20fd98318fca549d2881f1d75b512f3be'),
}


def file_digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def accepted_files(family, bundle=None):
    """Check the externally pinned archive and every member without extracting it."""
    require(family in ACCEPTED, 'unaccepted preparation family')
    expected = ACCEPTED[family]
    bundle = Path(bundle) if bundle is not None else HERE/'goal_20260924'/expected['bundle']
    receipt = json.loads((bundle/'receipt.json').read_text())
    archive = bundle/'evidence.tar.gz'
    require(archive.stat().st_size == expected['archive_bytes']
            and file_digest(archive) == expected['archive_sha256']
            and receipt['archive_bytes'] == expected['archive_bytes']
            and receipt['archive_sha256'] == expected['archive_sha256']
            and receipt['execution_sha256'] == expected['execution_sha256'],
            'preparation archive differs from accepted external binding')
    # Only these inputs survive the read. Every archive member still has its
    # bytes, mode and identity checked, including unused target evidence.
    selected = {'registration/execution.json', 'independent-audit.json',
                'execution/entry-output/candidate.json',
                'execution/entry-output/pipeline.stdout'} if family == 'f5' else {
                'registration/execution.json', 'independent-audit.json',
                'execution/entry-output/candidate.json',
                'execution/entry-output/summary.json',
                'execution/asset-files/fixture.json', 'execution/asset-files/base.json'}
    files, observed = {}, {}
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            role = PurePosixPath(member.name)
            require(member.isfile() and not role.is_absolute()
                    and '..' not in role.parts and role.as_posix() == member.name
                    and member.name not in observed
                    and 0 <= member.mode <= 0o777 and not member.mode & 0o022,
                    'unsafe or duplicate preparation archive member')
            data = tar.extractfile(member).read()
            observed[member.name] = dict(bytes=len(data), mode=member.mode,
                                        sha256=hashlib.sha256(data).hexdigest())
            if member.name in selected:
                files[member.name] = json.loads(data)
    require(observed == receipt['inventory'] and set(files) == selected,
            'preparation archive inventory differs')
    spec, audit = files['registration/execution.json'], files['independent-audit.json']
    require(sha256(spec) == expected['execution_sha256']
            and audit == receipt['result'] == json.loads((bundle/'AUDIT.json').read_text())
            and audit['candidate_id'] == expected['candidate_id']
            and audit['final_rank'] == 29 and audit['promotion_eligible'] is False
            and audit['online_speedup'] is None,
            'accepted preparation invocation or admission differs')
    admitted = (audit['complete_ic_admitted'] is True if family == 'f5' else
                audit['complete_source_bound'] is True and audit['solved_targets'] == 1)
    require(admitted, 'source control was not completely admitted')
    return files


def ordinary_inputs(family, files):
    """Whitelist preparation fields; discard all previous target answers and costs."""
    if family == 'f5':
        report = files['execution/entry-output/pipeline.stdout']
        fixture = report['fixture']
        base, logs = report['factor_base'], report['column_logs']
        attempts = []
        for batch in report['collection_reports']:
            for row in batch['attempts']:
                require(type(row['b']) is int and row['b'] == 0,
                        'ordinary preparation query depends on a target')
                attempts.append(dict(trial=row['trial'], scalar=row['a'],
                    outcome=row['pdp']['outcome'], indices=row['pdp']['points']))
    else:
        require(family == 'sat', 'unknown preparation family')
        report = files['execution/entry-output/summary.json']
        fixture = files['execution/asset-files/fixture.json']
        base = files['execution/asset-files/base.json']['factor_base']
        logs = report['column_logs']
        attempts = [dict(trial=row['trial'], scalar=row['probe_scalar'],
            outcome=row['status'], indices=(row['point_witness']['point_indices']
                if row['status'] == 'VALID_POINT_WITNESS' else None))
            for row in report['collection']]
    # Both accepted preparations were sampled independently of their supplied
    # target. Their previous workload and recovered scalar never enter ICP1.
    fixture = {key: copy.deepcopy(fixture[key]) for key in (
        'degree', 'curve_a', 'irreducible', 'subgroup_order', 'group_order',
        'cofactor', 'generator', 'lambda')}
    fixture.update(targets=[], target_seeds=[], target_scalar_constructed=False)
    return dict(fixture=fixture, base=base, logs=logs, attempts=attempts)


def reconstruct(inputs):
    fields(inputs, 'fixture base logs attempts', 'preparation inputs')
    fixture, raw_base = inputs['fixture'], inputs['base']
    fields(fixture, 'degree curve_a irreducible subgroup_order group_order cofactor '
           'generator lambda targets target_seeds target_scalar_constructed', 'preparation fixture')
    require(fixture['targets'] == [] and fixture['target_seeds'] == []
            and fixture['target_scalar_constructed'] is False,
            'preparation contains target input')
    curve = Curve(fixture)
    require(curve.n == 17 and curve.a == 1 and curve.r == 65587,
            'preparation certificate is limited to the accepted n17 control')
    base = [curve.decode(point) for point in raw_base]
    require(len(base) == 63 and None not in base and len(set(base)) == 63,
            'preparation geometry is missing, repeated or identity')
    matrix = RelationMatrix(curve, base)
    require(len(matrix.columns) == 29, 'preparation column count differs')
    inventory = factor_base_inventory(dict(factor_base=raw_base, columns=29), fixture)
    outcomes, trajectory = Counter(), []
    for trial, attempt in enumerate(inputs['attempts']):
        fields(attempt, 'trial scalar outcome indices', 'ordinary preparation attempt')
        require(type(attempt['trial']) is int and attempt['trial'] == trial,
                'ordinary preparation chronology differs')
        scalar = natural(attempt['scalar'], 'ordinary scalar')
        require(scalar < curve.r, 'ordinary preparation scalar outside subgroup')
        outcome, indices = attempt['outcome'], attempt['indices']
        require(outcome in ('witness', 'proved_unsat', 'VALID_POINT_WITNESS',
                           'SOURCE_UNSAT', 'CONFLICT_BUDGET_INCONCLUSIVE'),
                'unaccepted ordinary preparation outcome')
        outcomes[outcome] += 1
        if outcome in ('witness', 'VALID_POINT_WITNESS'):
            require(type(indices) is list, 'ordinary witness lacks indices')
            matrix.push(scalar, indices)
        else:
            require(indices is None, 'failed ordinary query contributes a relation')
        trajectory.append(matrix.rank)
    require(matrix.rank == len(matrix.columns), 'preparation is rank deficient')
    logs = matrix.solve()  # Solve from the ordinary rows, then replay every log.
    require(len(inputs['logs']) == len(logs), 'preparation log count differs')
    for item, point, log in zip(inputs['logs'], matrix.columns, logs):
        fields(item, 'point log', 'retained column log')
        require(curve.decode(item['point']) == point and int(item['log']) == log,
                'retained preparation column or logarithm differs')
    record = dict(policy=POLICY, curve=curve_record(fixture),
        factor_base=dict(points=[list(point) for point in base], inventory=inventory),
        projection=[None if row is None else list(row) for row in matrix.projected],
        column_logs=[dict(point=list(point), log=log)
                     for point, log in zip(matrix.columns, logs)])
    proof = dict(ordinary_query_count=len(inputs['attempts']),
        ordinary_outcome_mix=dict(sorted(outcomes.items())), rank_trajectory=trajectory,
        matrix=matrix.snapshot(), all_column_logs_independently_verified=True,
        preparation_target_count=0, previous_target_evidence_retained_in_state=False)
    return record, proof


def create(family, files):
    require(family in ACCEPTED, 'unaccepted preparation family')
    candidate = files['execution/entry-output/candidate.json']
    expected = ACCEPTED[family]
    require(candidate['candidate_id'] == expected['candidate_id']
            and sha256(candidate['record']) == candidate['record_sha256'] == expected['candidate_record_sha256']
            and candidate['candidate_id'].endswith('h'+candidate['record_sha256'][:12]),
            'preparation parent candidate identity differs')
    inputs = ordinary_inputs(family, files)
    require(sha256(inputs) == expected['inputs_sha256'],
            'ordinary preparation inputs differ from accepted source evidence')
    record, proof = reconstruct(inputs)
    require(record['curve'] == {key:candidate['record'][key] for key in ('field', 'curve')}
            and record['factor_base']['inventory'] == candidate['record']['factor_base']['inventory'],
            'preparation differs from parent curve or factor-base identity')
    h = sha256(record)
    return dict(schema_version=1, state_id='ICP1h'+h[:12], record_sha256=h,
        record=record, certificate=dict(inputs=inputs, proof=proof),
        provenance=dict(family=family, parent_candidate_id=candidate['candidate_id'],
            parent_candidate_record_sha256=candidate['record_sha256'],
            archive_sha256=expected['archive_sha256'], execution_sha256=expected['execution_sha256'],
            accepted_merge_commit=expected['merge_commit']),
        scope='mathematical preparation certificate; not a measured IC candidate or target result',
        native_execution=False, fresh_targets_generated=0, promotion_eligible=False,
        online_wall_ns=None, online_speedup=None)


def verify(document, expected_sha256):
    """Require an external whole-certificate seal as well as exact group replay."""
    require(sha256(document) == expected_sha256, 'preparation differs from external certificate seal')
    fields(document, 'schema_version state_id record_sha256 record certificate provenance scope '
           'native_execution fresh_targets_generated promotion_eligible online_wall_ns online_speedup',
           'preparation document')
    require(type(document['schema_version']) is int and document['schema_version'] == 1
            and document['native_execution'] is False
            and type(document['fresh_targets_generated']) is int and document['fresh_targets_generated'] == 0
            and document['promotion_eligible'] is False
            and document['online_wall_ns'] is None and document['online_speedup'] is None
            and document['scope'] == 'mathematical preparation certificate; not a measured IC candidate or target result',
            'preparation certificate claims a measurement or promotion')
    fields(document['certificate'], 'inputs proof', 'preparation certificate')
    record, proof = reconstruct(document['certificate']['inputs'])
    h = sha256(record)
    require(record == document['record'] and proof == document['certificate']['proof']
            and h == document['record_sha256'] == STATE_SHA256
            and document['state_id'] == 'ICP1h'+h[:12],
            'preparation state or independent relation proof differs')
    provenance = document['provenance']
    fields(provenance, 'family parent_candidate_id parent_candidate_record_sha256 archive_sha256 '
           'execution_sha256 accepted_merge_commit', 'preparation provenance')
    require(provenance['family'] in ACCEPTED, 'unaccepted preparation provenance family')
    expected = ACCEPTED[provenance['family']]
    require(all(provenance[key] == expected[key] for key in (
                'archive_sha256', 'execution_sha256'))
            and provenance['parent_candidate_id'] == expected['candidate_id']
            and provenance['accepted_merge_commit'] == expected['merge_commit']
            and provenance['parent_candidate_record_sha256'] == expected['candidate_record_sha256'],
            'preparation provenance differs from accepted parent')
    require(sha256(document['certificate']['inputs']) == expected['inputs_sha256'],
            'ordinary preparation inputs differ from accepted source evidence')
    return dict(status='VERIFIED_REUSABLE_PREPARATION', state_id=document['state_id'],
        certificate_sha256=expected_sha256, actual_usable_base_size=62,
        geometric_base_size=63, effective_columns=29, rank=proof['matrix']['rank'],
        ordinary_queries=proof['ordinary_query_count'],
        verified_relation_rows=proof['matrix']['accepted_rows'],
        target_input_present=False, native_execution=False, promotion_eligible=False,
        online_speedup=None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    freeze = sub.add_parser('freeze')
    freeze.add_argument('--family', choices=sorted(ACCEPTED), required=True)
    freeze.add_argument('--out', type=Path, required=True)
    replay = sub.add_parser('verify')
    replay.add_argument('--certificate', type=Path, required=True)
    replay.add_argument('--expected-sha256', required=True)
    args = parser.parse_args()
    if args.command == 'freeze':
        document = create(args.family, accepted_files(args.family))
        result = verify(document, sha256(document))
        write_immutable(args.out, document)
    else:
        result = verify(json.loads(args.certificate.read_text()), args.expected_sha256)
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
