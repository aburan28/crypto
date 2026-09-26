#!/usr/bin/env python3
"""Producer-side archived n19 witness and signed-row construction."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CORPUS = NOTES / 'rotated_pdp_corpus_20260925/evidence/raw.tar.gz'
SWEEP = NOTES / 'rotated_beta_sweep_20260925/evidence/raw.tar.gz'
GATE = NOTES / 'rotated_subspace_support_20260925/gate.py'
CORPUS_SHA = '39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c'
SWEEP_SHA = 'fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140'
BETAS = (3, 338435, 303097, 464276)
KNOWN_SUPPORT = (62389, 66179, 66203, 59323)
Q = 130873
H = (385982, 301867)
LAM = 41811
TORSION = (None, (0, 1), (1, 0), (1, 1))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')


def point(value):
    return None if value is None else tuple(value)


def pjson(value):
    return None if value is None else list(value)


def arithmetic():
    spec = importlib.util.spec_from_file_location('joint_producer_arithmetic', GATE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Data:
    def __init__(self):
        assert sha(CORPUS) == CORPUS_SHA and sha(SWEEP) == SWEEP_SHA
        mod = arithmetic()
        self.field = mod.Field(19, [0, 1, 2, 5])
        self.field.rabin_prime_degree()
        self.curve = mod.Curve(self.field)
        assert self.field.poly == 0x80027
        assert mod.source_group_order(19) == 4 * Q and mod.is_prime_by_trial(Q)
        assert self.curve.on_curve(H) and self.curve.scalar(H, Q) is None
        assert self.curve.tau(H) == self.curve.scalar(H, LAM)
        assert (LAM * LAM + LAM + 2) % Q == 0
        assert tuple(self.curve.scalar(t, 4) for t in TORSION) == (None,) * 4
        self.arms = {}
        self.bytes_scanned = 0
        with tarfile.open(CORPUS, 'r:gz') as corpus, tarfile.open(SWEEP, 'r:gz') as sweep:
            for beta, known in zip(BETAS, KNOWN_SUPPORT):
                archive, prefix = ((corpus, 'raw/n19-m6/') if beta == 3
                                   else (sweep, f'raw/beta-{beta}/'))
                factors = [[point(p) for p in slot]
                           for slot in json.loads(self.read(archive, prefix + 'factors.json'))]
                assert len(factors) == 6 and all(len(slot) == 7 for slot in factors)
                assert all(self.curve.on_curve(p) for slot in factors for p in slot)
                assert all(factors[i + 1] == sorted(self.curve.tau(p) for p in factors[i])
                           for i in range(5))
                raw_hist = self.read(archive, prefix + 'projected_histogram.jsonl')
                hist = {}
                total = 0
                for line in raw_hist.splitlines():
                    entry = json.loads(line)
                    target = point(entry['point'])
                    assert target not in hist and entry['count'] > 0
                    assert len(entry['witness_indices']) == 6
                    assert all(0 <= j < 7 for j in entry['witness_indices'])
                    assert self.curve.on_curve(point(entry['full_sum']))
                    hist[target] = entry
                    total += entry['count']
                assert len(hist) == known and total == 7 ** 6
                terms = []
                all_cols = set()
                for i, slot in enumerate(factors):
                    slot_terms = []
                    for source in slot:
                        back = source
                        for _ in range((19 - i) % 19):
                            back = self.curve.tau(back)
                        projected = self.curve.scalar(back, 4)
                        if projected is None:
                            assert back in TORSION and source in TORSION
                            slot_terms.append((None, 0))
                            continue
                        canonical = min(projected, self.curve.neg(projected))
                        sign = 1 if canonical == projected else -1
                        coefficient = sign * pow(LAM, i, Q) % Q
                        assert self.curve.scalar(source, 4) == self.curve.scalar(canonical, coefficient)
                        all_cols.add(canonical)
                        slot_terms.append((canonical, coefficient))
                    terms.append(slot_terms)
                base_cols = {c for c, _ in terms[0] if c is not None}
                assert len(base_cols) == 3 and all_cols == base_cols
                self.arms[beta] = {'factors': factors, 'hist': hist, 'terms': terms,
                                   'columns': sorted(base_cols)}
        self.columns = sorted({c for arm in self.arms.values() for c in arm['columns']})

    def read(self, tar: tarfile.TarFile, name: str) -> str:
        member = tar.getmember(name)
        assert member.isfile() and member.size < 40_000_000
        handle = tar.extractfile(member)
        assert handle is not None
        raw = handle.read()
        assert len(raw) == member.size
        self.bytes_scanned += len(raw)
        return raw.decode('utf8')

    def row(self, beta: int, target_q):
        target_r = self.curve.scalar(target_q, 4)
        archived = self.arms[beta]['hist'].get(target_r)
        if archived is None:
            return None
        factors = self.arms[beta]['factors']
        indices = archived['witness_indices']
        sources = [factors[i][j] for i, j in enumerate(indices)]
        total = None
        for source in sources:
            total = self.curve.add(total, source)
        assert total == point(archived['full_sum'])
        assert self.curve.scalar(total, 4) == target_r
        torsion = self.curve.add(total, self.curve.neg(target_q))
        assert torsion in TORSION and self.curve.add(target_q, torsion) == total
        coefficients = {}
        for i, j in enumerate(indices):
            column, coeff = self.arms[beta]['terms'][i][j]
            if column is not None:
                coefficients[column] = (coefficients.get(column, 0) + coeff) % Q
        terms = [{'column': list(column), 'coefficient': value}
                 for column, value in sorted(coefficients.items()) if value]
        evaluated = None
        for term in terms:
            evaluated = self.curve.add(evaluated, self.curve.scalar(tuple(term['column']),
                                                                     term['coefficient']))
        assert evaluated == target_r
        return {'beta': beta, 'point': pjson(target_q), 'projected_point': pjson(target_r),
                'torsion_index': TORSION.index(torsion), 'row': terms,
                'source_witness_indices': indices, 'archived_full_sum': pjson(total)}


def vector(row, columns):
    by_col = {col: i for i, col in enumerate(columns)}
    result = [0] * len(columns)
    for term in row['row']:
        j = by_col[tuple(term['column'])]
        result[j] = (result[j] + term['coefficient']) % Q
    return result
