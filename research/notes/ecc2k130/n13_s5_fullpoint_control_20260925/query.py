#!/usr/bin/env python3
"""Exact Q_i+T_j full-point query export; construction only, no solver call."""
from __future__ import annotations

import sys
from pathlib import Path

from s5 import Circuit, N, MODULUS, build, read_factors, read_targets

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from export import write_cnf  # noqa: E402
sys.path.insert(0, str(HERE.parent / 'n19_s6_fullpoint_exporter_20260925'))
from admission import check_cnf  # noqa: E402


def circuit() -> Circuit:
    return build(N, MODULUS, read_factors())


def write_query(c: Circuit, label: str, path: Path, *, byte_cap: int) -> dict:
    matches = [row for row in read_targets() if row['id'] == label]
    if len(matches) != 1:
        raise ValueError('unknown or duplicate frozen branch label')
    target = tuple(matches[0]['full_target'])
    if path.exists():
        raise FileExistsError('refusing to overwrite an earlier exact query')
    units = c.target_units(target)
    meta = write_cnf(c, path, units=units, byte_cap=byte_cap)
    checked = check_cnf(c, path, units)
    if checked['sha256'] != meta['sha256']:
        raise AssertionError('query hash drift between writer and independent parser')
    return {'label': label, 'target': list(target), 'units': units,
            'cnf': meta, 'checked': checked}
