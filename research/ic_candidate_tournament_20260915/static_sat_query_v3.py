"""Complete-source SAT query with independently replayable native invocations."""
import time

from oracle import require
from run_cms_s4_controls import (digest, lift_source_assignment, parse_cms_model,
                                 sha_bytes, validate_xor_dimacs)
from run_static_cms_s4_controls import validate_export
from static_sat_native_v3 import meter
from tournament import read


def export_arguments(panel, item):
    return ['17','6','standard',str(panel['export_nonce']),
            str(panel['cms_conflict_budget']),'instance','1','0',
            '--target-x',str(item['point'][0]),'--target-y',str(item['point'][1]),
            '--blind-instance-id',f"control-{item['trial']:02d}",'--export-only']


def cms_arguments(panel):
    return ['--verb','1','--threads','1','--random','1','--maxsol','1',
            '--maxconfl',str(panel['cms_conflict_budget']),'instance/instance.xor.cnf']


def one_query(panel, item, execution, curve, base, out):
    """Keep capped/failed attempts and measure verification separately.

    The caller measures the whole function interval. ``verification_wall_ns``
    is an exclusive subset of that interval; the remainder belongs to PDP
    construction/solver/dispatch, including fresh child-process overhead.
    """
    trial=item['trial']
    directory=out/f'trial-{trial:02d}'
    directory.mkdir()
    if item['point'] is None:
        return dict(trial=trial,probe_scalar=item['probe_scalar'], public_point=None,
                    exporter=None,cms=None,source_model_valid=None,point_witness=None,
                    status='IDENTITY_QUERY',verification_wall_ns=0)
    instance=directory/'instance'
    command=export_arguments(panel,item)
    exported=meter(execution,'bin/exporter',command,directory,'export',
                   panel['export_timeout_seconds'])
    row=dict(trial=trial,probe_scalar=item['probe_scalar'],
             public_point=item['point'],exporter=exported,cms=None,
             source_model_valid=None,point_witness=None,
             status='EXPORT_FAILURE',verification_wall_ns=0)
    if exported['timed_out'] or exported['returncode']!=0:
        return row
    check_start=time.monotonic_ns()
    try:
        manifest=read(instance/'manifest.json')
        row['exports']=validate_export(manifest,item,panel,base,instance)
        row['manifest_sha256']=digest(instance/'manifest.json')
    except (OSError,ValueError,KeyError,TypeError) as error:
        row.update(status='INVALID_EXPORT',
                   reason=f'{type(error).__name__}: {error}')
        row['verification_wall_ns']+=time.monotonic_ns()-check_start
        return row
    row['verification_wall_ns']+=time.monotonic_ns()-check_start
    command=cms_arguments(panel)
    measured=meter(execution,'bin/cms',command,directory,'cms',
                   panel['cms_timeout_seconds'])
    row['cms']=measured
    check_start=time.monotonic_ns()
    stdout=(directory/'cms.stdout').read_text()
    if measured['timed_out']:
        row['status']='TIMEOUT'
    elif measured['returncode']==10 and 's SATISFIABLE' in stdout:
        maximum=manifest['exports']['cryptominisat_xor_dimacs']['variables']
        model=parse_cms_model(stdout,maximum)
        if model is None or not validate_xor_dimacs(
                instance/'instance.xor.cnf',model):
            row['status']='INVALID_SOURCE_MODEL'
        else:
            row['source_model_valid']=True
            row['source_model_sha256']=sha_bytes(bytes(model))
            witness=lift_source_assignment(model,manifest,curve,base,
                                           tuple(item['point']))
            row['point_witness']=witness
            row['status']=('VALID_POINT_WITNESS' if witness['group_replay']
                           else 'SOURCE_MODEL_NONLIFTING')
    elif measured['returncode']==20 and 's UNSATISFIABLE' in stdout:
        row['status']='SOURCE_UNSAT'
    elif (measured['returncode']==15 and 's INDETERMINATE' in stdout
          and not measured['timed_out']):
        row['status']='CONFLICT_BUDGET_INCONCLUSIVE'
    elif measured['returncode']==0:
        row['status']='UNKNOWN_INCONCLUSIVE'
    else:
        row['status']='SOLVER_ERROR'
    row['verification_wall_ns']+=time.monotonic_ns()-check_start
    require(row['verification_wall_ns']>=0,'negative query verification clock')
    return row
