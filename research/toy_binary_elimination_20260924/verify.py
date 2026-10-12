"""Replay both frozen stages and independently check linear consequences."""
import hashlib
import json
from pathlib import Path
import run
import run_initial

ROOT = Path(__file__).resolve().parent


def assert_subset(a, b):
    if isinstance(a, dict):
        for k, v in a.items():
            assert_subset(v, b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            assert_subset(x, y)
    else:
        assert a == b


def dense_rank(rows):
    rows = [list(r) for r in rows]
    pivot = 0
    for col in range(7):
        found = next((i for i in range(pivot, len(rows)) if rows[i][col]), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        for i in range(pivot+1, len(rows)):
            if rows[i][col]:
                rows[i] = [a ^ b for a, b in zip(rows[i], rows[pivot])]
        pivot += 1
    return pivot


def main():
    original = json.loads((ROOT/'results.json').read_text())
    refined = json.loads((ROOT/'refined.json').read_text())
    assert json.loads(json.dumps(run_initial.run())) == original
    assert json.loads(json.dumps(run.run())) == refined
    assert_subset(original, refined)
    for record in refined['records']:
        zeros = record['polynomial_zeros']
        evaluation_rows = [[1] + [(z >> i) & 1 for i in range(6)] for z in zeros]
        expected = 7 - dense_rank(evaluation_rows)
        assert record['elimination']['stages'][-1]['linear_consequence_dimension'] == expected
    # Independent elementary ideals: zero ideal, unit ideal, and (b0).
    assert run.matrix_profile([0]*64, list(range(64)))['linear_completion_degree'] == 0
    assert run.matrix_profile([1]+[0]*63, [])['contradiction_degree'] == 0
    coordinate = [0]*64
    coordinate[1] = 1
    test = run.matrix_profile(coordinate, list(range(0,64,2)))
    assert test['contradiction_degree'] is None and test['linear_completion_degree'] == 1
    receipt_path = ROOT/'receipt.json'
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        for relative, expected in receipt['sha256'].items():
            assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() == expected
    print('PASS: two frozen replays, 106 independent linear-dimension checks, three elementary ideals.')


if __name__ == '__main__':
    main()
