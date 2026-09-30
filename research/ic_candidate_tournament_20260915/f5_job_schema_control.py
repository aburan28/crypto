"""Independent expected Job deserialization results; no IC query or oracle."""
import copy
from pathlib import Path

from generic_stages import DEFAULTS
from oracle import require

HERE = Path(__file__).resolve().parent
CONTROL = HERE/'goal_20260924/f5-source-bound-runtime-v2/interface-control'


def cases(job):
    """Freeze the corrected supplied-point job and adversarial schema controls."""
    require(job['target_seeds'] == [], 'v2 unseeded job must have an empty input seed list')
    rows = [dict(label='corrected-empty-seeds', input=copy.deepcopy(job))]
    omitted = copy.deepcopy(job)
    del omitted['target_seeds']
    rows.append(dict(label='omitted-default-seeds', input=omitted))
    for label, mutate in (
        ('historical-null-seed', lambda value: value.update(target_seeds=[None])),
        ('negative-seed', lambda value: value.update(target_seeds=[-1])),
        ('algorithm-seed-overflow', lambda value: value.update(algorithm_seed=2**64)),
        ('unknown-job-field', lambda value: value.update(unknown_flag=True)),
        ('unknown-config-field', lambda value: value['config'].update(unknown_flag=True)),
        ('wrong-coordinate-shape', lambda value: value.update(public_targets=[['52411']])),
    ):
        value = copy.deepcopy(job)
        mutate(value)
        rows.append(dict(label=label, input=value))
    return rows


def audit(inputs, exported):
    require(exported['schema_version'] == 1 and exported['status'] == 'NATIVE_JOB_SCHEMA_CONTROL'
            and len(exported['results']) == len(inputs) == 8,
            'native Job control result domain differs')
    corrected = inputs[0]['input']
    require(inputs == cases(corrected), 'native Job control input cases differ')
    expected = dict(copy.deepcopy(corrected),
                    config=dict(copy.deepcopy(DEFAULTS), **corrected['config']))
    for number, (case, result) in enumerate(zip(inputs, exported['results'])):
        require(result['label'] == case['label'], 'native Job case order differs')
        if number < 2:
            require(result == dict(label=case['label'], status='PARSED', observed=expected),
                    'native Job/Config deserialization differs from registered job/defaults')
        else:
            require(set(result) == {'label', 'status', 'error'} and result['status'] == 'REJECTED'
                    and type(result['error']) is str and result['error'],
                    'invalid Job case was accepted or lost its error')
    return dict(schema_version=1, status='PASS_NATIVE_JOB_SCHEMA_ONLY', parsed_cases=2,
        rejected_cases=6, scope='actual pinned Job/Config deserialization only; execute_job never called',
        complete_ic_admitted=False, candidate_id=None, measured_costs=None,
        online_speedup=None, promotion_eligible=False)
