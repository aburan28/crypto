"""Versioned bounded-retention driver for the user-required legacy lab CLI.

Original runner and worker files remain unchanged. Every journal record and
hash-chain link is verified by the original canonical codec; retain only the
last state. Native worker timing/resource supervision is unchanged.
"""
from pathlib import Path
import hashlib
import json
import sys
from lab import runner as base

OriginalLaboratory = base.Laboratory

class JournalTail:
    def __init__(self, count, last):
        self.count, self.last = count, last
    def __len__(self):
        return self.count
    def __getitem__(self, index):
        if index != -1 or not self.count:
            raise IndexError('bounded journal view only exposes the verified tail')
        return self.last

class StreamingLaboratory(OriginalLaboratory):
    def __init__(self, ledger, process_runner=None):
        super().__init__(ledger, process_runner or self.native_process)
        self.driver_receipt = None
        self.active_experiment = None
    def journal(self):
        path = self.ledger / 'journal.jsonl'
        previous = '0' * 64
        count = 0
        last = record = None
        whole = hashlib.sha256()
        size = 0
        if path.exists():
            with path.open('rb') as handle:
                for number, line in enumerate(handle, 1):
                    last = record = None
                    whole.update(line)
                    size += len(line)
                    try:
                        record = base.strict_json(line)
                        checksum = record.pop('record_sha256')
                    except (ValueError, KeyError) as error:
                        raise base.LabError(f'journal record {number} is incomplete or invalid; preserve evidence') from error
                    if (record.get('sequence') != number or record.get('previous_sha256') != previous
                            or checksum != base.digest(base.canonical(record))):
                        raise base.LabError(f'journal record {number} fails integrity')
                    previous = checksum
                    record['record_sha256'] = checksum
                    last = record
                    count = number
        self.driver_receipt = dict(schema='rootfinding-streaming-driver/v1', journal_records=count,
            journal_bytes=size, journal_sha256=whole.hexdigest(), last_record_sha256=previous,
            streaming_driver_sha256=base.digest(Path(__file__).read_bytes()),
            original_runner_sha256=base.digest(Path(base.__file__).read_bytes()),
            original_runner_unmodified=True, worker_supervision='unchanged original run_cold_process')
        return JournalTail(count, last)
    def run_one(self, eid, timeout, **kwargs):
        self.active_experiment = eid
        return super().run_one(eid, timeout, **kwargs)
    def native_process(self, command, **kwargs):
        item = self.state['experiments'][self.active_experiment]
        directory = self.ledger / item['attempts'][-1]['directory']
        before = base.digest(Path(__file__).read_bytes())
        receipt = {**self.driver_receipt, 'streaming_driver_sha256':before,
                   'scientific_identity':item['registered_identity'], 'change':'journal memory retention only'}
        base.atomic_json(directory/'coordinator_driver.json', receipt, exclusive=True)
        result = base.run_cold_process(command, **kwargs)
        if base.digest(Path(__file__).read_bytes()) != before:
            raise base.LabError('streaming driver changed during worker execution')
        return result

if __name__ == '__main__':
    base.Laboratory = StreamingLaboratory
    raise SystemExit(base.cli())
