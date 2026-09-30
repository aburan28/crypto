#!/usr/bin/env python3
"""Lock the v3 smoke registration and refuse to measure it.

Seed 2026093001 stays reserved. This checker does not generate targets,
restore archives, or start the tournament. A later commit may enable one
dispatch only after lost-v2-campaign-exposures.json is on main at the
sealed hash. Seeds 2026092901 and 2026092902 are not retried.
"""
import json
from pathlib import Path
import sys

from generic_solver_feasibility import assess
from oracle import require
from tournament import digest, read

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE / 'goal_20260924/generic-backend-qualification-v3'
PANEL = REGISTRATION / 'panel.json'
PANEL_SHA256 = 'df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4'
V2 = HERE / 'goal_20260924/generic-backend-qualification-v2'
LOST_V1 = V2 / 'lost-campaign-exposures.json'
LOST_V1_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'
LOST_V2 = V2 / 'lost-v2-campaign-exposures.json'
LOST_V2_SHA256 = '0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54'
# This file has no path that starts a campaign. Dispatch stays false until a
# later change adds one, after the v2 exposure corpus is merged.
DISPATCH_AUTHORIZED = False


def registration_check(panel):
    """Reject a panel that is not the frozen planning registration."""
    require(digest(PANEL) == PANEL_SHA256, 'changed premeasurement v3 panel bytes')
    require(panel['schema_version'] == 1 and panel['status'] == 'REGISTERED_PLANNING'
            and panel['seed'] == 2026093001
            and panel['comparison_kind'] == 'factor-base-policy'
            and panel['stages'] == ['aa', 'smoke']
            and panel['cells'] == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0']
            and panel['holdout_cells_excluded'] == ['n29a1']
            and panel['repetitions'] == 1 and panel['timeout_seconds'] == 300
            and panel['memory_bytes'] == 8 * 1024**3
            and panel['measure_timeout_minutes'] == 180
            and panel['pack_timeout_minutes'] == 35
            and panel['job_timeout_minutes'] == 240
            and panel['scheduled_pair_bound'] == 45 and panel['pair_cap'] == 60
            and panel['sat_in_scope'] is False and panel['promotion_eligible'] is False
            and panel['worker_commit'] == '765c3c5f19032bd852163805f257c56babef2040',
            'changed registered v3 schedule or resources')
    ids = [row['id'] for row in panel['candidates']]
    require(ids == ['incumbent', 'prepared_both', 'generic_pair_subspace_dense',
                    'generic_f4_subspace_dense', 'generic_f5_subspace_dense'],
            'changed registered v3 candidate ids')
    for row in panel['candidates'][2:]:
        base = row['config']['factor_base']
        require(row['adapter'] == 'generic-v1' and row['source'] == 'controlled-generic'
                and base == {'kind': 'standard_subspace', 'dimension': 6},
                'changed subspace arm '+row['id'])
    for row in panel['candidates'][3:]:
        config = row['config']
        require(config['summands'] == 3 and config['groebner_degree'] == 3
                and config['node_budget'] == 4096 and config['max_trials'] == 256
                and config['batch_trials'] == 8 and config['linear_algebra'] == 'dense',
                'changed algebraic budget '+row['id'])
    require(panel['rho_arms'] == [
        {'id': 'rho', 'role': 'cold_rho_reference', 'config': {'rho_parallel_walks': 8}},
        {'id': 'rho_online', 'role': 'online_rho_reference',
         'config': {'rho_parallel_walks': 16, 'linear_algebra': 'dense'}}],
        'changed registered v3 rho roles')
    require(digest(LOST_V1) == LOST_V1_SHA256, 'changed frozen first-run exposures')
    layout = assess(panel)
    require(layout['status'] == 'PASS_STATIC_LAYOUT_ONLY'
            and layout['impossible_algebraic_cells'] == 0,
            'v3 algebraic layout is outside the encoder cap')
    return layout


def exposure_block():
    """Why seed 2026093001 must not sample points on this checkout."""
    if not LOST_V2.is_file():
        return 'v2 exposure census is not on this checkout'
    if digest(LOST_V2) != LOST_V2_SHA256:
        return 'v2 exposure census hash does not match the sealed registration'
    lost = read(LOST_V2)
    if lost.get('seed') != 2026092902:
        return 'v2 exposure census is not seed 2026092902'
    return None


def status_report(panel):
    block = exposure_block()
    return {
        'seed': panel['seed'],
        'panel_sha256': PANEL_SHA256,
        'static_layout': 'PASS_STATIC_LAYOUT_ONLY',
        'v2_exposure_sha256': LOST_V2_SHA256,
        'v2_exposure_present': block is None,
        'dispatch_block': block,
        'dispatch_authorized': DISPATCH_AUTHORIZED,
        'measurement': 'not_run',
    }


def main():
    panel = read(PANEL)
    registration_check(panel)
    report = status_report(panel)
    require(report['dispatch_authorized'] is False, 'v3 dispatch flag was enabled in this checker')
    print(json.dumps(report, sort_keys=True), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
