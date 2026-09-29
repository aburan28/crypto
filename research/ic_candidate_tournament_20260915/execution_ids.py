"""Explicit execution numbering for diagnostic exports; no algorithm identity changes."""
import argparse

from identity import natural, run_id
from oracle import require

BLOCK_SIZE = 1000
LANES = ('controls', 'integration')


def ci_start(workflow_run, attempt, lane):
    """Injective run/attempt namespaces, with two disjoint bounded lanes."""
    natural(workflow_run, 'workflow run', positive=True)
    natural(attempt, 'workflow attempt', positive=True)
    require(lane in LANES, 'unknown execution lane')
    total = workflow_run + attempt
    namespace = total * (total + 1) // 2 + attempt  # Cantor pairing
    return (len(LANES) * namespace + LANES.index(lane)) * BLOCK_SIZE


def allocation(start, labels):
    natural(start, 'run number start')
    require(0 < len(labels) <= BLOCK_SIZE, 'execution allocation exceeds block capacity')
    require(all(type(label) is str and label for label in labels)
            and len(set(labels)) == len(labels), 'execution labels must be unique nonempty strings')
    return dict(schema_version=1, start=start, end_exclusive=start + len(labels),
                executions=[dict(label=label, number=start + i) for i, label in enumerate(labels)])


def audit_runs(records):
    keys = set()
    for record in records:
        cid, wid, rid = (record[key] for key in ('candidate_id', 'workload_id', 'run_id'))
        require(type(rid) is str and rid.rsplit('R', 1)[-1].isascii()
                and rid.rsplit('R', 1)[-1].isdigit(), 'invalid execution number')
        require(rid == run_id(cid, wid, int(rid.rsplit('R', 1)[-1])), 'run key identity mismatch')
        key = (cid, wid, rid)
        require(key not in keys, 'duplicate canonical run key')
        keys.add(key)
    return dict(status='PASS', run_records=len(records), unique_keys=len(keys))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-run', type=int, required=True)
    parser.add_argument('--attempt', type=int, required=True)
    parser.add_argument('--lane', choices=LANES, required=True)
    args = parser.parse_args()
    print(ci_start(args.workflow_run, args.attempt, args.lane))
