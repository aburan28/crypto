"""Scientific stage admission using the explicit bounded negative-proof adapter."""
from generic_stages import effective_config, verify_base, dispatch, matrix_audit
from generic_query_law_exact_v1 import verify_query_law
from oracle import require, verify


def verify_stages(report, fixture, job):
    require(job.get('exclusive_phases') is True and type(report.get('generic_admission_schema')) is int
            and report['generic_admission_schema'] == 1,
            'missing scientific stage evidence')
    require(report.get('generic_runtime_policy') == 'default-environment-one-rayon-v1',
            'missing runtime override policy')
    require(report['fixture'] == fixture and job.get('public_targets') == fixture['targets']
            and len(fixture['targets']) == 1, 'stage fixture mismatch')
    require(job['mode'] in {'ic', 'inventory'} and report['mode'] == 'ic', 'stage mode mismatch')
    cfg = effective_config(job, report)
    base = verify_base(report, fixture, job)
    inventory = job['mode'] == 'inventory'
    law = observed_dispatch = certificate = None
    if not inventory:
        law = verify_query_law(report, fixture, job)
        observed_dispatch = dispatch(report, cfg)
        if report['status'] == 'complete':
            certificate = verify(report, fixture, summands=cfg['summands'])
    matrix = matrix_audit(report, fixture, cfg, inventory=inventory)
    return dict(schema_version=1, status='PASS', base=base, query_law=law,
                dispatch=observed_dispatch, matrix=matrix, certificate=certificate,
                scope='base, query, dispatch, matrix and correctness; build/accounting required separately',
                promotion_eligible=False)
