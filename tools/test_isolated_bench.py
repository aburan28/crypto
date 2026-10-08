"""Tests for tools/isolated_bench.py. Run: python3 -m pytest tools/test_isolated_bench.py"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import isolated_bench as ib  # noqa: E402

TOOL = str(Path(__file__).with_name('isolated_bench.py'))


def spare_cpu():
    cpus = sorted(os.sched_getaffinity(0))
    if len(cpus) < 2:
        pytest.skip('needs at least two CPUs')
    return cpus[-1]


def test_parse_cpus():
    assert ib.parse_cpus('3') == {3}
    assert ib.parse_cpus('0-2,5') == {0, 1, 2, 5}
    with pytest.raises(ValueError):
        ib.parse_cpus('')


def test_other_use_excludes_and_ignores_idle():
    before = {1: ('a', 10), 2: ('b', 5), 3: ('c', 7)}
    after = {1: ('a', 10 + ib.TICK), 2: ('b', 5 + 2 * ib.TICK), 3: ('c', 7), 4: ('d', ib.TICK)}
    used = ib.other_use(before, after, exclude={1})
    assert used == {'b[2]': 2.0, 'd[4]': 1.0}


def test_monitor_excludes_descendants_from_its_cpu_snapshot():
    # The benchmark is present in the CPU scan, then exits before a separate
    # process-tree scan could find it. Parent IDs from the same scan retain it.
    snapshot = {
        10: ('taskset', 1, 1),
        11: ('python3', 5, 10),
        12: ('f4_f2_bench', 3 * ib.TICK, 11),
        13: ('background', ib.TICK, 1),
    }
    assert ib.descendants(10, snapshot) == {11, 12}
    now = {pid: (name, ticks) for pid, (name, ticks, _) in snapshot.items()}
    used = ib.other_use({}, now, exclude={10} | ib.descendants(10, snapshot))
    assert used == {'background[13]': 1.0}


def test_reserving_every_cpu_is_refused():
    with pytest.raises(SystemExit):
        ib.check_cpus(set(os.sched_getaffinity(0)))


def test_run_pins_records_and_restores(tmp_path):
    cpu = spare_cpu()
    out = tmp_path / 'r.jsonl'
    lock = tmp_path / 'lock'
    probe = 'import os,json;print(json.dumps(sorted(os.sched_getaffinity(0))))'
    result = subprocess.run([sys.executable, TOOL, 'run', '--cpus', str(cpu), '--settle', '0.5',
                             '--max-other-cpu', '0.5', '--lock', str(lock), '--out', str(out), '--',
                             sys.executable, '-c', probe], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [cpu]
    record = json.loads(out.read_text())
    run = record['run']
    assert run['cpus'] == [cpu] and run['exit_status'] == 0
    for key in ('wall_seconds', 'involuntary_switches', 'other_cpu_seconds', 'contended'):
        assert key in run
    assert record['left_on_reserved']['user_threads'] == []
    assert cpu in os.sched_getaffinity(os.getpid())


def test_an_inherited_narrow_mask_is_widened_and_recorded(tmp_path):
    # A harness forked while an isolated run had moved its parent off the
    # benchmark CPU inherits a mask without it; the run must still go ahead.
    cpu = spare_cpu()
    allowed = set(os.sched_getaffinity(0))
    out = tmp_path / 'r.jsonl'
    probe = 'import os,json;print(json.dumps(sorted(os.sched_getaffinity(0))))'
    result = subprocess.run([sys.executable, TOOL, 'run', '--cpus', str(cpu), '--settle', '0.5',
                             '--max-other-cpu', '0.5', '--lock', str(tmp_path / 'lock'), '--out', str(out),
                             '--', sys.executable, '-c', probe], capture_output=True, text=True,
                            preexec_fn=lambda: os.sched_setaffinity(0, allowed - {cpu}))
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [cpu]
    record = json.loads(out.read_text())
    assert record['affinity_widened_from'] == sorted(allowed - {cpu})
    assert cpu in os.sched_getaffinity(os.getpid())


def test_a_run_with_the_cpu_in_its_mask_records_no_widening(tmp_path):
    cpu = spare_cpu()
    out = tmp_path / 'r.jsonl'
    result = subprocess.run([sys.executable, TOOL, 'run', '--cpus', str(cpu), '--settle', '0.5',
                             '--max-other-cpu', '0.5', '--lock', str(tmp_path / 'lock'), '--out', str(out),
                             '--', 'true'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'affinity_widened_from' not in json.loads(out.read_text())


def test_busy_machine_is_refused(tmp_path):
    cpu = spare_cpu()
    hog = subprocess.Popen([sys.executable, '-c', 'while True: pass'])
    try:
        time.sleep(0.5)
        result = subprocess.run([sys.executable, TOOL, 'run', '--cpus', str(cpu), '--settle', '1',
                                 '--lock', str(tmp_path / 'lock'), '--out', str(tmp_path / 'r.jsonl'),
                                 '--', 'true'], capture_output=True, text=True)
    finally:
        hog.kill()
        hog.wait()
    assert result.returncode != 0
    assert 'machine is busy' in result.stderr
    assert not (tmp_path / 'r.jsonl').exists()
    assert cpu in os.sched_getaffinity(os.getpid())
