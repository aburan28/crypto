"""Freeze source, method, curve, and one-target workload for static SAT IC.

The registration files are generated once, before any complete pipeline run.
Their source manifest deliberately hashes the whole tournament Python surface:
this binds transitive helpers without guessing which imports are material.
"""
import argparse
import ast
import copy
import hashlib
from pathlib import Path

from identity import (candidate_manifest, canonical, run_id, sha256,
                      workload_manifest, write_immutable)
from oracle import Curve, require
from run_generic_exact_yield_audit import PANEL as EXACT_PANEL, load_evidence, record
from tournament import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REGISTRATION = HERE/'goal_20260924/static-sat-full'
PANEL = REGISTRATION/'panel.json'


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def runtime_paths():
    """Hash the local transitive Python imports of the executed controller."""
    roots = [HERE/'run_static_sat_full.py',
             HERE/'static_sat_registration.py',
             ROOT/'scripts/process_meter.py',
             ROOT/'scripts/run_koblitz_pdp_matrix.py',
             ROOT/'scripts/verify_stage20_phase_b_terminal_evidence.py']
    queued = list(roots)
    seen = set()
    while queued:
        path = queued.pop()
        if path in seen:
            continue
        require(path.is_file(), 'executed Python component is missing')
        seen.add(path)
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                modules = [node.module]
            else:
                continue
            for module in modules:
                name = module.split('.')[0]
                for directory in (HERE, ROOT/'scripts'):
                    local = directory/(name+'.py')
                    if local.is_file() and local not in seen:
                        queued.append(local)
    return sorted(seen)


def source_manifest():
    panel = read(PANEL)
    parent_files = load_evidence(read(EXACT_PANEL))
    parent_source = record(parent_files, 'build/source-manifest.json')
    components = [dict(role=str(path.relative_to(ROOT)),
                       sha256=file_digest(path)) for path in runtime_paths()]
    components.extend((
        dict(role='source-bound-Rust-exporter',
             sha256=panel['exporter_source_sha256']),
        dict(role='source-receipted-static-CryptoMiniSat-executable',
             sha256=panel['cms_executable_sha256']),
        dict(role='static-CryptoMiniSat-build-receipt',
             sha256=panel['cms_build_receipt_sha256']),
        dict(role='static-CryptoMiniSat-build-bundle-seal',
             sha256=panel['cms_build_bundle_seal_sha256']),
    ))
    return dict(schema_version=1, components=components,
                parent_rust_source_commit=panel['source_commit'],
                parent_rust_source_manifest_sha256=sha256(parent_source),
                parent_exact_result_sha256=panel['parent_exact_result_sha256'],
                stage_encoding=panel['source_encoding'])


def method_record(panel, source_sha256, fixture):
    source = source_sha256
    return dict(
        isogeny='none',
        endomorphism=dict(order_conductor=None,
                          frobenius_order_conductor=None,
                          volcano_levels=[]),
        factor_base=dict(construction=dict(
            recipe='standard-subspace', dimension=6,
            basis_bitmasks=[1, 2, 4, 8, 16, 32],
            ordering='constructor-abscissa-order-then-half-trace-root-and-negative',
            subgroup_map='cofactor-multiplication-remove-identity-deduplicate'),
            nominal_bound=6),
        point_decomposition=dict(
            summands=3, solver='sat', summation_polynomial='symmetrised-S4',
            encoding='source-bound-wide-S4-circuit-XOR-DIMACS',
            equation_order='source-defined-wide-S4-circuit-order',
            monomial_order='none', internal_matrix_kernel='CryptoMiniSat-native-XOR-CDCL',
            limits=dict(conflict_budget=panel['cms_conflict_budget'],
                        watchdog_seconds=panel['cms_timeout_seconds'],
                        exporter_watchdog_seconds=panel['export_timeout_seconds'],
                        models_per_query=panel['cms_max_models_per_query']),
            cache_policy='source-exporter-process-per-query; SAT-process-per-query',
            source_sha256=source),
        relation_collection=dict(
            collector='sample',
            query_distribution='StdRng08-trial-keyed-nonzero-uniform-subgroup-scalar',
            query_rule='public-aG; one-ordinary-query-per-trial; fresh-run-seed',
            filtering='independent-source-clause-XOR-model-lift-and-group-readd',
            verification='complete-source-model-and-exact-group-readd',
            duplicates='scalar-and-sorted-base-indices',
            dependencies='retain-unique-dependent-rows',
            stop_rule=dict(max_queries=panel['max_relation_queries'],
                           success='full-29-column-rank-and-all-logs-scalar-replayed'),
            source_sha256=source),
        relation_linear_algebra=dict(
            solver='gauss', modulus=int(fixture['subgroup_order']),
            matrix_construction='cofactor-project; sorted-sign-Frobenius-columns; rhs=h*a; sum-signed-lambda-powers',
            orbit_quotient='sign-and-Frobenius',
            rank_criterion='exact-modular-full-column-rank-and-all-rows-checked',
            block_parameters='none', preconditioner='none',
            source_sha256=source),
        target_descent=dict(
            method='pdp', policy='StdRng08-independent-aG+bQ-pairs',
            recursive_solvers='same-source-SAT-PDP; no-recursion',
            success_rule='relation-derived-scalar-and-independent-group-replay',
            stop_rule=dict(max_queries=panel['max_descent_queries']),
            source_sha256=source),
        implementation=dict(
            source_manifest_sha256=source,
            components=[dict(role='complete-executed-Python-and-external-source-manifest',
                             sha256=source)],
            flags=dict(cms_binary_sha256=panel['cms_executable_sha256'],
                       cms_build_receipt_sha256=panel['cms_build_receipt_sha256'],
                       rust_exporter_source_sha256=panel['exporter_source_sha256'],
                       model_limit=panel['cms_max_models_per_query'],
                       source_encoding=panel['source_encoding'])))


def identities(panel):
    files = load_evidence(read(EXACT_PANEL))
    report = record(files, 'jobs/n17a1/f5/stdout.json')
    fixture = report['fixture']
    curve = Curve(fixture)
    target = curve.decode(panel['target_input']['point'])
    require(target is not None and curve.mul(target, curve.r) is None,
            'registered target is not a subgroup public point')
    source = source_manifest()
    method = method_record(panel, sha256(source), fixture)
    candidate = candidate_manifest(fixture, report, method)
    target_fixture = copy.deepcopy(fixture)
    target_fixture.update(targets=[list(target)],
                          target_seeds=[panel['target_input']['seed']],
                          target_scalar_constructed=False)
    workload = workload_manifest(
        target_fixture,
        input_law='one-supplied-public-point; SHA256-x-lift-cofactor-v1-fixture',
        algorithm_seed=panel['target_input']['seed'],
        resource_envelope=dict(host_class='physical-macos-arm64',
                               cpu_workers=1, target_count=1,
                               memory_limit_bytes=None,
                               total_wall_limit_seconds=21600),
        cache_policy='warm')
    return source, method, candidate, workload


def register():
    panel = read(PANEL)
    require(panel['candidate_id'] is None
            and all(not (REGISTRATION/name).exists() for name in (
                'source-manifest.json', 'method.json', 'candidate.json',
                'workload.json', 'seal.json')),
            'static SAT candidate is already registered')
    source, method, candidate, workload = identities(panel)
    panel['candidate_id'] = candidate['candidate_id']
    panel['workload_id'] = workload['workload_id']
    panel['run_id'] = run_id(candidate['candidate_id'], workload['workload_id'], 0)
    write_immutable(REGISTRATION/'source-manifest.json', source)
    write_immutable(REGISTRATION/'method.json', method)
    write_immutable(REGISTRATION/'candidate.json', candidate)
    write_immutable(REGISTRATION/'workload.json', workload)
    PANEL.write_bytes(canonical(panel)+b'\n')
    seal = dict(schema_version=1, panel_sha256=file_digest(PANEL),
                method_sha256=file_digest(REGISTRATION/'method.json'),
                candidate_sha256=file_digest(REGISTRATION/'candidate.json'),
                workload_sha256=file_digest(REGISTRATION/'workload.json'),
                source_manifest_file_sha256=file_digest(
                    REGISTRATION/'source-manifest.json'),
                candidate_id=candidate['candidate_id'],
                workload_id=workload['workload_id'], run_id=panel['run_id'])
    write_immutable(REGISTRATION/'seal.json', seal)
    return seal


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register', action='store_true', required=True)
    parser.parse_args()
    print(register())
