#!/usr/bin/env python3
"""Inventory historical seed identities; initialize only after live-owner review.

Run without --initialize first. The bootstrap is additive: existing identity
records are never overwritten, including their active owners or progress.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import re
import struct
import time

from seed_registry import CONFIG_KEY, PREFIX, modal_stream, run_key, slot_stream


def record_run_ids(body):
    offset, stride = (16, 72) if body.startswith(b'ECC2KDP2') else (0, 32)
    if len(body) < offset or (len(body) - offset) % stride:
        raise ValueError('unaligned historical corpus; identity audit incomplete')
    ids = {struct.unpack_from('<Q', body, i)[0] >> 48 for i in range(offset, len(body), stride)}
    for rid in ids:
        run_key(rid)
    return ids


def inventory(s3, bucket):
    def listing(prefix, delimiter=None):
        args = dict(Bucket=bucket, Prefix=prefix)
        if delimiter:
            args['Delimiter'] = delimiter
        return s3.get_paginator('list_objects_v2').paginate(**args)

    def document(key):
        return json.loads(s3.get_object(Bucket=bucket, Key=key)['Body'].read())

    now = time.time()
    prefixes = [''] + [p['Prefix'] for page in listing('campaigns/', '/')
                       for p in page.get('CommonPrefixes', [])]
    records, collisions, namespaces, unknown = {}, [], [], []
    for prefix in prefixes:
        try:
            config = document(prefix + 'campaign.json')
        except s3.exceptions.NoSuchKey:
            continue
        if config.get('curve') != 131:
            continue
        base = int(config.get('runIdBase', 1))
        namespaces.append(dict(prefix=prefix, run_id_base=base))
        slot_keys = [o['Key'] for page in listing(prefix + 'slots/')
                     for o in page.get('Contents', []) if o['Key'].endswith('.json')]
        with ThreadPoolExecutor(max_workers=12) as pool:
            slot_docs = dict(zip(slot_keys, pool.map(document, slot_keys)))
        checkpoints = {}
        for page in listing(prefix + 'ckpt/'):
            for obj in page.get('Contents', []):
                match = re.search(r'/slot-(\d+)(?:/[0-9a-f]{64})?\.ck$', '/' + obj['Key'])
                if not match:
                    continue
                slot = int(match[1])
                if slot not in checkpoints or obj['LastModified'] > checkpoints[slot]['LastModified']:
                    checkpoints[slot] = obj

        candidates = {}
        for key, doc in slot_docs.items():
            slot = int(re.search(r'slot-(\d+)\.json$', key)[1])
            candidates[slot] = doc
        for slot in checkpoints:
            candidates.setdefault(slot, {})
        for page in listing(prefix + 'dp/', '/'):
            for p in page.get('CommonPrefixes', []):
                match = re.search(r'/slot-(\d+)/$', '/' + p['Prefix'])
                if match:
                    candidates.setdefault(int(match[1]), {})

        def checkpoint(obj):
            h = s3.get_object(Bucket=bucket, Key=obj['Key'], Range='bytes=0-39')['Body'].read()
            if len(h) != 40 or h[:8] != b'ECC2K130':
                raise ValueError('unreadable checkpoint identity: ' + obj['Key'])
            version, curve, threads, batch, lanes, rid, iteration = struct.unpack('<6IQ', h[8:])
            if curve != 131 or not all((version, threads, batch, lanes)):
                raise ValueError('invalid checkpoint identity: ' + obj['Key'])
            return rid, iteration
        with ThreadPoolExecutor(max_workers=12) as pool:
            headers = dict(zip(checkpoints, pool.map(checkpoint, checkpoints.values())))
        for slot, doc in sorted(candidates.items()):
            modal = prefix == '' and slot >= 90000
            rid = slot - 90000 if modal else base + slot
            if not 1 <= rid <= 65535:
                # Old manual uploads need not use either slot mapping. Read
                # every record there, not a sample or a guessed arithmetic offset.
                objects = [o for page in listing(prefix + 'dp/slot-%05d/' % slot)
                           for o in page.get('Contents', []) if o['Key'].endswith('.bin')]
                if slot in headers or not objects:
                    raise ValueError('unmapped historical checkpoint or empty slot %s%d' % (prefix, slot))
                def ids_in_object(obj):
                    return record_run_ids(s3.get_object(Bucket=bucket, Key=obj['Key'])['Body'].read())
                with ThreadPoolExecutor(max_workers=8) as pool:
                    ids = set().union(*pool.map(ids_in_object, objects))
                unknown.append(dict(prefix=prefix, slot=slot, run_ids=sorted(ids),
                                    objects=len(objects), bytes=sum(o['Size'] for o in objects)))
                continue
            if slot in headers and headers[slot][0] != rid:
                raise ValueError('configuration disagrees with historical checkpoint run id')
            stream = modal_stream(rid) if modal else slot_stream(prefix, slot)
            active = doc.get('state') == 'active' and doc.get('leaseUntil', 0) > now
            row = dict(version=1, run_id=rid, stream=stream, started=True,
                       checkpoint_floor=max(int(doc.get('ckptIter', 0)), headers.get(slot, (rid, 0))[1], 0),
                       created_at=now, imported=True,
                       active_owner=('legacy:' + str(doc['owner'])) if active else None,
                       blocked=doc.get('state') in ('retired', 'error', 'solved'))
            if rid in records:
                prior = records[rid]
                collisions.append(dict(run_id=rid, canonical_stream=prior['stream'],
                                       conflicting_stream=stream, conflicting_active=active))
                # Only the already-investigated, stopped legacy Modal shadows
                # may alias the original root slots. Never pick a live winner.
                if not (modal and 1 <= rid <= 4 and not active and prior['stream'] == slot_stream('', rid - 1)):
                    raise ValueError('unresolved historical seed collision for run %d' % rid)
                prior['checkpoint_floor'] = max(prior['checkpoint_floor'], row['checkpoint_floor'])
                prior.setdefault('legacy_conflicts', []).append(stream)
            else:
                records[rid] = row
    for source in unknown:
        for rid in source['run_ids']:
            alias = 'legacy-dp:%s:slot%d' % (source['prefix'], source['slot'])
            if rid not in records:
                records[rid] = dict(version=1, run_id=rid, stream=alias, started=True,
                                    checkpoint_floor=0, created_at=now, imported=True,
                                    active_owner=None, blocked=True)
            row = records[rid]
            if row.get('active_owner'):
                raise ValueError('unmapped corpus overlaps active run %d; review before initializing' % rid)
            row['blocked'] = True
            row.setdefault('legacy_upload_aliases', []).append(alias)
    return dict(checked_at=now, bucket=bucket, namespaces=namespaces,
                records=list(records.values()), historical_conflicts=collisions,
                unmapped_corpora=unknown)


def initialize(s3, audit):
    bucket = audit['bucket']
    # Readiness is the last write. A partial bootstrap cannot authorize starts.
    for row in audit['records']:
        key = run_key(row['run_id'])
        try:
            s3.put_object(Bucket=bucket, Key=key, Body=json.dumps(row, sort_keys=True).encode(),
                          ContentType='application/json', IfNoneMatch='*')
        except Exception as exc:
            code = getattr(exc, 'response', {}).get('Error', {}).get('Code')
            if code not in ('PreconditionFailed', '412', 'ConditionalRequestConflict'):
                raise
            existing = json.loads(s3.get_object(Bucket=bucket, Key=key)['Body'].read())
            if existing.get('run_id') != row['run_id'] or existing.get('stream') != row['stream']:
                raise ValueError('existing registry binding differs; bootstrap refused')
    config = dict(version=1, state='ready', imported_records=len(audit['records']),
                  initialized_at=time.time(), automatic_owner_expiry=False)
    try:
        s3.put_object(Bucket=bucket, Key=CONFIG_KEY, Body=json.dumps(config).encode(),
                      ContentType='application/json', IfNoneMatch='*')
    except Exception as exc:
        if getattr(exc, 'response', {}).get('Error', {}).get('Code') not in ('PreconditionFailed', '412'):
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bucket', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--initialize', action='store_true')
    args = parser.parse_args()
    import boto3
    s3 = boto3.client('s3')
    audit = inventory(s3, args.bucket)
    with open(args.out, 'w') as output:
        json.dump(audit, output, indent=2, sort_keys=True)
    print(json.dumps(dict(run_ids=len(audit['records']),
                          active_owners=sum(bool(r.get('active_owner')) for r in audit['records']),
                          historical_conflicts=audit['historical_conflicts'])))
    if args.initialize:
        initialize(s3, audit)
        print('shared seed registry initialized; existing bindings preserved')


if __name__ == '__main__':
    main()
