"""Versioned JSON Redis artifacts for the ONB research harness.

Checksums detect accidental corruption; the Redis writer namespace is trusted.
No pickle or executable cache content. Credentials are read from the environment.
"""
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import time
import ir


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


class ArtifactCache:
    def __init__(self, client=None, ttl=3600, maxBytes=32*1024*1024):
        self.client, self.ttl, self.maxBytes = client, ttl, maxBytes
        if not 1 <= ttl <= 31536000 or maxBytes <= 0:
            raise ValueError('invalid cache limits')
        self.stats = Counter()
        root = Path(__file__).resolve().parent
        names = ['onb_artifacts.py', 'indexcalc_e2e.py', 'field.py', 'curves.py',
                 'decomp.py', 'build.py', 'ir.py', 'cnf.py', 'indexcalc.py']
        self.version = hashlib.sha256(b''.join((root/name).read_bytes() for name in names)).hexdigest()

    @classmethod
    def fromEnvironment(cls):
        import redis
        url = os.environ.get('ONB_REDIS_URL')
        if not url:
            raise ValueError('ONB_REDIS_URL is required for Redis mode')
        return cls(redis.Redis.from_url(url, socket_timeout=1, socket_connect_timeout=1,
                                       retry_on_timeout=False))

    def key(self, kind, parameters):
        return 'ic-onb-v1:' + kind + ':' + hashlib.sha256(encode([self.version, parameters])).hexdigest()

    def getOrBuild(self, kind, parameters, producer):
        key = self.key(kind, parameters)
        start = time.perf_counter_ns()
        if self.client is not None:
            try:
                # Bound response size atomically before GET; no oversized payload transfer.
                raw = self.client.eval("local n=redis.call('STRLEN',KEYS[1]); if n>tonumber(ARGV[1]) then return false end; return redis.call('GET',KEYS[1])", 1, key, self.maxBytes)
                if raw:
                    envelope = json.loads(raw)
                    payload = envelope['payload']
                    if envelope['key'] != key or envelope['sha256'] != hashlib.sha256(encode(payload)).hexdigest():
                        raise ValueError('cache integrity')
                    self.stats['hits'] += 1
                    self.stats['read_bytes'] += len(raw)
                    return payload
            except Exception:
                self.stats['read_errors'] += 1
            finally:
                self.stats['lookup_ns'] += time.perf_counter_ns() - start
        self.stats['misses'] += 1
        start = time.perf_counter_ns()
        payload = producer()
        self.stats['producer_ns'] += time.perf_counter_ns() - start
        self.stats['producers'] += 1
        raw = encode(dict(key=key, sha256=hashlib.sha256(encode(payload)).hexdigest(), payload=payload))
        if len(raw) > self.maxBytes:
            self.stats['oversized'] += 1
        elif self.client is not None:
            start = time.perf_counter_ns()
            try:
                self.client.set(key, raw, ex=self.ttl)
                self.stats['written_bytes'] += len(raw)
            except Exception:
                self.stats['write_errors'] += 1
            finally:
                self.stats['publish_ns'] += time.perf_counter_ns() - start
        return payload

    def circuit(self, m, nring, points, leaf, producer):
        def make():
            prog, roots = producer()
            return dict(ops=prog.ops, inputs=prog.inputRef, roots=roots, nInput=prog.nInput)
        data = self.getOrBuild('semaev-chain', [m, nring, points, leaf], make)
        prog = ir.Prog()
        prog.ops = [None if op is None else (op[0], tuple(op[1])) for op in data['ops']]
        prog.inputRef = [None if ref is None else tuple(ref) for ref in data['inputs']]
        prog.nInput = data['nInput']
        prog.hash = {op: i for i, op in enumerate(prog.ops) if op is not None}
        return prog, data['roots']

    def base(self, m, nring, prime, eigen, weight, producer):
        def make():
            reps, lookup = producer()
            return dict(reps=reps, lookup=[[list(p), col, coeff] for p, (col, coeff) in sorted(lookup.items())])
        data = self.getOrBuild('subgroup-base', [m, nring, str(prime), str(eigen), weight], make)
        return [tuple(p) for p in data['reps']], {tuple(p): (col, coeff) for p, col, coeff in data['lookup']}
