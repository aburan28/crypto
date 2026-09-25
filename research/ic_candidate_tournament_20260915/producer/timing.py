"""Check native worker intervals; never substitute them for instruction costs."""
from identity import natural
from measurement import PHASES, report_sha256
from oracle import require

TARGET_CHILDREN = ('target_query', 'target_pdp', 'target_relation_check')
RAW_PHASES = (set(PHASES)-{'isogeny'}) | set(TARGET_CHILDREN) | {'reference_solve'}


def native_intervals(report, process_wall_ns):
    require(report.get('phase_schema') == 3 and report.get('status') == 'complete',
            'native interval requires a complete schema-v3 report')
    require(report.get('mode') in ('ic', 'rho'), 'unknown timed algorithm')
    require(len(report['fixture']['targets']) == 1, 'online interval is single-target only')
    source = report.get('diagnostics', report)
    require(source.get('target_input') == 'supplied_public_point', 'fixture generation entered measured worker')
    wall = source['phase_wall_ns']
    require(type(wall) is dict and set(wall) == RAW_PHASES, 'missing/extra native phase interval')
    for phase, value in wall.items():
        natural(value, phase+' nanoseconds')
    natural(process_wall_ns, 'cold process nanoseconds', positive=True)
    snapshot = sum(wall.values())
    require(0 < snapshot <= process_wall_ns, 'worker intervals exceed native process boundary')
    if report['mode'] == 'ic':
        require(wall['reference_solve'] == 0, 'rho work in an IC interval')
        online_names = (*TARGET_CHILDREN, 'target_descent', 'recovery_check')
        cold = {p: wall[p] for p in PHASES if p != 'isogeny'}
        cold['target_descent'] += sum(wall[p] for p in TARGET_CHILDREN)
        cold['isogeny'] = 0
        cold['setup'] += process_wall_ns-snapshot
    else:
        require(all(wall[p] == 0 for p in RAW_PHASES-{'setup', 'reference_solve', 'recovery_check'}),
                'IC work in a rho interval')
        online_names = ('reference_solve', 'recovery_check')
        cold = {p: wall[p] for p in ('setup', 'reference_solve', 'recovery_check')}
        cold['setup'] += process_wall_ns-snapshot
    online = {p: wall[p] for p in online_names}
    natural(source['online_wall_ns'], 'online nanoseconds', positive=True)
    require(sum(online.values()) == source['online_wall_ns'], 'online phase intervals do not close')
    require(sum(cold.values()) == process_wall_ns, 'cold native phase intervals do not close')
    return {'schema_version': 1, 'unit': 'native_monotonic_ns', 'target_count': 1,
            'native_report_sha256': report_sha256(report),
            'online': {'phase_wall_ns': online, 'wall_ns': source['online_wall_ns'],
                'boundary': ('after certified reusable factor logs; target query through final scalar replay'
                             if report['mode'] == 'ic' else
                             'supplied public point; target-dependent rho setup, solve and final scalar replay'),
                'target_generation_included': False, 'scalar_replay_included': True},
            'cold': {'phase_wall_ns': cold, 'wall_ns': process_wall_ns,
                'worker_snapshot_wall_ns': snapshot, 'external_setup_remainder_ns': process_wall_ns-snapshot,
                'remainder_policy': 'setup charges process launch/input and report/exit tail outside worker snapshot',
                'unattributed_wall_ns': 0}}
