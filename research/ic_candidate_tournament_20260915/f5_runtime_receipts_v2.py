"""Publication representation for raw resource diagnostics, never IC costs."""
import copy
import math

from oracle import require

SECONDS = ('wall_seconds', 'user_seconds', 'system_seconds',
           'total_core_seconds', 'single_core_seconds')


def publication_process(process):
    """Keep integer native clocks; encode finite raw seconds as decimal strings.

    The unchanged original metrics file remains in the archive and is hashed
    by the auditor. No rounding, unit conversion or performance admission is
    applied to the native_wall_ns clock or to any scientific phase.
    """
    result = copy.deepcopy(process)
    metrics = result.pop('metrics')
    require(set(metrics) == {*SECONDS, 'peak_rss_bytes', 'meter'}
            and type(metrics['peak_rss_bytes']) is int and metrics['peak_rss_bytes'] >= 0
            and type(metrics['meter']) is str,
            'unexpected native resource diagnostic schema')
    require(all(type(metrics[key]) in (int, float) and math.isfinite(metrics[key])
                and metrics[key] >= 0 for key in SECONDS),
            'invalid native resource seconds')
    result['resource_diagnostics'] = dict(peak_rss_bytes=metrics['peak_rss_bytes'],
        meter=metrics['meter'], seconds_decimal={key: repr(metrics[key]) for key in SECONDS},
        interpretation='raw uncalibrated process resource diagnostics; not scientific phase costs')
    return result
