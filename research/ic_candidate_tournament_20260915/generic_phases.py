"""Exclusive generic-worker phase closure; not source/performance admission.

Legacy profiler labels are rejected. Unentered intervals stay null. A complete
clock partition still needs source/dispatch, group, query, base and matrix checks.
"""
import re
from pathlib import Path

from oracle import require

PHASES = ('setup', 'factor_base', 'precompute', 'queries', 'pdp', 'relation_check',
          'matrix_build', 'relation_la', 'target_query', 'target_pdp',
          'target_relation_check', 'target_descent', 'recovery_check', 'rho_solve')
ONLINE = ('target_query', 'target_pdp', 'target_relation_check', 'target_descent',
          'recovery_check', 'rho_solve')
REQUIRED = {'ic': set(PHASES) - {'rho_solve'},
            'rho': {'setup', 'precompute', 'rho_solve', 'recovery_check'}}


def natural(value, label, *, positive=False):
    require(type(value) is int and value >= (1 if positive else 0), 'invalid ' + label)
    return value


def verify_native(report, job, *, process_wall_ns=None):
    require(job.get('exclusive_phases') is True, 'exclusive phase mode not declared')
    mode = job['mode']
    require(mode in REQUIRED and report['mode'] == mode, 'phase mode mismatch')
    require(report['status'] in {'complete', 'incomplete'}, 'invalid terminal phase status')
    require(report.get('generic_phase_policy') == 'exclusive-owner-thread-v1' and
            report.get('online_timing_schema') == 2, 'unknown generic phase policy')
    require(job.get('public_targets') == report['fixture']['targets']
            and len(job['public_targets']) == 1, 'phase target mismatch')
    require(report['reusable_setup_excluded'] is True and
            report['target_input'] == 'supplied_public_point', 'unknown online accounting boundary')
    trace = report['generic_phase_timing']
    require(trace['schema_version'] == 1, 'unknown generic trace schema')
    cold, online = trace['phases_ns'], trace['online_phases_ns']
    require(type(cold) is dict and set(cold) == set(PHASES), 'missing generic phase slots')
    require(type(online) is dict and set(online) == set(ONLINE), 'missing online phase slots')
    for phase, cost in cold.items():
        if cost is not None:
            natural(cost, phase + ' wall ns')
        require(phase in REQUIRED[mode] or cost is None, 'inapplicable phase gained measured work')
    observed = natural(trace['observed_wall_ns'], 'observed wall ns', positive=True)
    require(cold['setup'] is not None, 'missing setup interval')
    require(sum(value for value in cold.values() if value is not None) == observed,
            'observed phase clocks do not close')
    online_ns = trace['online_wall_ns']
    require(online_ns == report['online_wall_ns'], 'headline online clock differs from phase clock')
    for phase, cost in online.items():
        if cost is not None:
            natural(cost, phase + ' online ns')
        require(cost == cold[phase], 'target-dependent work occurred outside the online interval')
    if online_ns is None:
        require(all(value is None for value in online.values()), 'unstarted online interval gained work')
        require(report['outer_online_wall_ns'] is None, 'unstarted online interval gained outer clock')
        require(report['scalar_replay_included'] is False, 'unstarted target claims scalar replay')
    else:
        natural(online_ns, 'online wall ns', positive=True)
        require(sum(value for value in online.values() if value is not None) == online_ns,
                'online phase clocks do not close')
        outer = natural(report['outer_online_wall_ns'], 'outer online clock', positive=True)
        require(outer <= online_ns <= observed, 'online interval exceeds its enclosing clocks')
    missing = sorted(phase for phase in REQUIRED[mode] if cold[phase] is None)
    if report['status'] == 'complete':
        require(not missing and online_ns is not None and report['scalar_replay_included'] is True,
                'complete result lacks required phase/replay coverage')
    process_phases = dict(cold)
    remainder = None
    if process_wall_ns is not None:
        natural(process_wall_ns, 'process wall ns', positive=True)
        require(process_wall_ns >= observed, 'worker clock exceeds whole process')
        remainder = process_wall_ns - observed
        process_phases['setup'] += remainder
    # A failed target never gains a complete scientific result, even if some
    # intervals are present. Its measured costs are retained as diagnostics.
    complete = report['status'] == 'complete' and not missing
    return dict(schema_version=1, unit='native_wall_ns', observed_phases_ns=dict(cold),
                observed_wall_ns=observed, online_phases_ns=dict(online), online_wall_ns=online_ns,
                process_phases_ns=process_phases, external_setup_remainder_ns=remainder,
                process_wall_ns=process_wall_ns, missing_required_phases=missing,
                complete_phase_coverage=complete,
                cold_wall_ns=process_wall_ns if complete else None,
                scope='exclusive clock closure only', promotion_eligible=False)


def parse_profiles(directory, report, job):
    """Check exclusive Callgrind intervals against both trace coverage and Ir checksum."""
    coverage = verify_native(report, job)
    directory = Path(directory)
    paths = sorted(directory.glob('callgrind.out*'))
    require(paths and all(p.is_file() and p.suffix != '.gz' for p in paths),
            'missing/unexpected generic instruction profiles')
    phases, parts, terminated = {}, set(), 0
    for path in paths:
        lines = path.read_text().splitlines()

        def one(prefix):
            fields = [line[len(prefix):].strip() for line in lines if line.startswith(prefix)]
            require(len(fields) == 1, 'missing/duplicate generic profile field ' + prefix)
            return fields[0]

        require(one('events:') == 'Ir', 'incompatible generic instruction unit')
        part = int(one('part:'))
        require(part >= 1 and part not in parts, 'duplicate/invalid generic profile part')
        parts.add(part)
        total = natural(int(one('summary:')), 'generic profile total')
        require(int(one('totals:')) == total, 'generic profile summary/totals mismatch')
        trigger = one('desc: Trigger:')
        if trigger == 'Program termination':
            phase = 'setup'
            terminated += 1
        else:
            prefix = 'Client Request: generic_ic_'
            require(trigger.startswith(prefix), 'legacy/unknown generic profiling boundary')
            phase = trigger[len(prefix):]
            require(phase in PHASES, 'unknown generic instruction phase')
        phases[phase] = phases.get(phase, 0) + total
    require(terminated == 1 and parts == set(range(1, len(parts) + 1)),
            'missing/duplicate final or interior generic profile interval')
    entered = {phase for phase, value in coverage['observed_phases_ns'].items() if value is not None}
    require(set(phases) == entered, 'native clock/profile phase coverage mismatch')
    collected = re.findall(r'Collected\s*:\s*(\d+)', (directory / 'stderr.txt').read_text())
    require(len(collected) == 1, 'missing whole-process generic instruction checksum')
    total = natural(int(collected[0]), 'whole-process instructions', positive=True)
    require(sum(phases.values()) == total, 'generic phase instructions do not close')
    values = {phase: phases.get(phase) for phase in PHASES}
    return dict(schema_version=1, unit='valgrind-3.22-amd64-Ir', phases=values,
                process_instructions=total,
                cold_instructions=total if coverage['complete_phase_coverage'] else None,
                online_instructions=(sum(phases.get(phase, 0) for phase in ONLINE)
                                     if report['online_wall_ns'] is not None else None),
                complete_phase_coverage=coverage['complete_phase_coverage'],
                scope='exclusive instruction closure only', promotion_eligible=False)
