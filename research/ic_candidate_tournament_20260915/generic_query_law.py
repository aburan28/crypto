"""Independent query-law replay for the bounded generic worker, version 1.

The integer stream follows pinned rand 0.8.8, rand_core 0.6.4 and
rand_chacha 0.3.1: PCG seed expansion, ChaCha12, then u64 sample_single.
This is a reproducibility adapter, not a randomness or security claim.
It does not import or invoke the measured producer. Rust crate vector controls
and actual-worker records exercise the independently implemented stream.
"""
import hashlib
import json
from pathlib import Path

from generic_queries import verify_queries
from oracle import require

MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
STREAM_VERSION = 'rand-0.8.8/rand_core-0.6.4/rand_chacha-0.3.1/StdRng-u64-v1'
SOLVERS = {'pair_table', 'enumerate', 'f4', 'f5', 'inherited_f4', 'sat_xor', 'sat_cnf'}


def uint(value, bits, name):
    require(type(value) is int and 0 <= value < (1 << bits), f'invalid {name}')
    return value


def rotate(value, shift, bits):
    mask = (1 << bits) - 1
    shift %= bits
    return ((value << shift) | (value >> ((bits - shift) % bits))) & mask


class StdRng08:
    """Small scalar implementation of the exact pinned integer stream."""

    def __init__(self, words):
        require(len(words) == 8, 'invalid ChaCha seed length')
        self.key = [uint(word, 32, 'ChaCha seed word') for word in words]
        self.counter = 0
        self.buffer = iter(())

    @classmethod
    def seed_from_u64(cls, seed):
        state = uint(seed, 64, 'RNG seed')
        words = []
        for _ in range(8):
            state = (state * 6364136223846793005 + 11634580027462260723) & MASK64
            x = (((state >> 18) ^ state) >> 27) & MASK32
            words.append(rotate(x, -(state >> 59), 32))
        return cls(words)

    def _block(self):
        initial = [0x61707865, 0x3320646e, 0x79622d32, 0x6b206574] + self.key + [
            self.counter & MASK32, self.counter >> 32, 0, 0]
        words = initial.copy()

        def quarter(a, b, c, d):
            words[a] = (words[a] + words[b]) & MASK32
            words[d] = rotate(words[d] ^ words[a], 16, 32)
            words[c] = (words[c] + words[d]) & MASK32
            words[b] = rotate(words[b] ^ words[c], 12, 32)
            words[a] = (words[a] + words[b]) & MASK32
            words[d] = rotate(words[d] ^ words[a], 8, 32)
            words[c] = (words[c] + words[d]) & MASK32
            words[b] = rotate(words[b] ^ words[c], 7, 32)

        for _ in range(6):
            for indices in ((0, 4, 8, 12), (1, 5, 9, 13), (2, 6, 10, 14), (3, 7, 11, 15),
                            (0, 5, 10, 15), (1, 6, 11, 12), (2, 7, 8, 13), (3, 4, 9, 14)):
                quarter(*indices)
        self.counter = (self.counter + 1) & MASK64
        return iter((a + b) & MASK32 for a, b in zip(words, initial))

    def next_u32(self):
        value = next(self.buffer, None)
        if value is None:
            self.buffer = self._block()
            value = next(self.buffer)
        return value

    def next_u64(self):
        return self.next_u32() | (self.next_u32() << 32)

    def nonzero_below(self, order):
        uint(order, 64, 'subgroup order')
        require(order >= 2, 'degenerate subgroup order')
        width = order - 1
        zone = ((width << (64 - width.bit_length())) - 1) & MASK64
        while True:
            product = self.next_u64() * width
            if (product & MASK64) <= zone:
                return 1 + (product >> 64)


def probe_scalar(seed, trial, order):
    uint(seed, 64, 'algorithm seed')
    uint(trial, 64, 'trial')
    key = seed ^ 0x50524f4245534551 ^ rotate((trial * 0x9e3779b97f4a7c15) & MASK64, 17, 64)
    return StdRng08.seed_from_u64(key).nonzero_below(order)


def collection_coefficients(seed, order, walked):
    stride = (StdRng08.seed_from_u64(seed ^ 0x5354524944455f30).nonzero_below(order)
              if walked else None)
    trial = 0
    anchor = None
    while True:
        if walked:
            if trial % 64 == 0:
                anchor = probe_scalar(seed, trial // 64, order)
            a = (anchor + (trial % 64) * stride) % order
        else:
            a = probe_scalar(seed, trial, order)
        yield a, 0
        trial += 1


def descent_coefficients(seed, order, walked):
    rng = StdRng08.seed_from_u64(seed ^ (0x57414c4b44455300 if walked else 0x44455343454e5400))
    if walked:
        a0, b = rng.nonzero_below(order), rng.nonzero_below(order)
        stride = max(order // 64, 1 << 20)
        trial = 0
        while True:
            yield (a0 + (trial % 64) * stride + trial // 64) % order, b
            trial += 1
    else:
        while True:
            yield rng.nonzero_below(order), rng.nonzero_below(order)


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def verify_query_law(report, fixture, job):
    """Check all recorded queries against the externally frozen worker job.

    The ordinary group/ledger checker runs first. A passing receipt does not
    establish a rank-based stop condition or certify an implementation/source.
    """
    require(job['mode'] == 'ic', 'query-law adapter requires IC')
    degree = uint(job['degree'], 32, 'degree')
    require(5 <= degree <= 31 and degree % 2 == 1, 'query-law adapter degree unsupported')
    a = uint(job['curve_a'], 8, 'curve coefficient')
    require(a <= 1 and degree == fixture['degree'] and a == fixture['curve_a'],
            'query-law curve mismatch')
    require(job.get('public_targets') == fixture['targets'] and len(fixture['targets']) == 1,
            'query-law requires the same one supplied point')
    seed = uint(job['algorithm_seed'], 64, 'algorithm seed')
    config = job['config']
    solver = config.get('solver', 'pair_table')
    require(solver in SOLVERS, 'query-law solver unsupported')
    summands = uint(config.get('summands', 3), 64, 'summands')
    require(2 <= summands <= 4 and report['summands'] == summands, 'query-law summands mismatch')
    batch_size = uint(config.get('batch_trials', 64), 64, 'batch size')
    limit = uint(config.get('max_trials', 4096), 64, 'trial limit')
    require(1 <= batch_size <= min(limit, 4096) and limit <= 65536, 'query-law limits unsupported')
    window = config.get('collection_window')
    if window is not None:
        uint(window, 64, 'collection window')
    order = int(fixture['subgroup_order'])
    group_receipt = verify_queries(report, fixture, summands)
    # Every supported bounded field has the packed backend. A pair-table worker
    # cannot start without its pair table. The actual geometric point count
    # controls the window predicate (not orbit columns or the usable census B).
    walked_collection = (solver == 'pair_table' and summands == 3 and window is not None
                         and 0 < window < len(report['factor_base']))
    expected = collection_coefficients(seed, order, walked_collection)
    trial = 0
    for batch in report['collection_reports']:
        require(trial < limit and batch['trials'] == min(batch_size, limit - trial),
                'query-law batch partition mismatch')
        for attempt in batch['attempts']:
            require((attempt['a'], attempt['b']) == next(expected),
                    f'query-law collection coefficient mismatch at {trial}')
            trial += 1
    require(0 < trial <= limit, 'query-law collection budget mismatch')
    walked_descent = solver == 'pair_table'
    solutions = report.get('solutions', [])
    require(report['status'] in {'complete', 'incomplete'}, 'query-law invalid terminal status')
    require(len(solutions) <= 1, 'query-law unexpected target count')
    if report['status'] == 'complete':
        require(len(solutions) == 1, 'query-law completed target missing')
    elif not solutions:
        require(trial == limit, 'query-law failed preparation stopped early')
    for solution in solutions:
        require(solution['index'] == 0 and 0 < solution['trials'] <= limit,
                'query-law descent budget mismatch')
        expected = descent_coefficients(seed, order, walked_descent)
        for t, attempt in enumerate(solution['attempts']):
            require((attempt['a'], attempt['b']) == next(expected),
                    f'query-law descent coefficient mismatch at {t}')
        if solution['recovered'] is None:
            require(solution['trials'] == limit, 'query-law failed descent stopped early')
    return dict(schema_version=1, stream_version=STREAM_VERSION,
                collection_law='windowed-64' if walked_collection else 'sampled',
                descent_law='walked-64' if walked_descent else 'sampled',
                job_sha256=canonical_digest(job), query_sha256=group_receipt['query_sha256'],
                checker_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in ('generic_query_law.py', 'generic_queries.py', 'oracle.py')},
                collection_queries=trial, descent_queries=group_receipt['descent_queries'],
                scope='query law and group replay only', promotion_eligible=False)


def verify_rust_vectors(vectors):
    require(vectors['schema_version'] == 1, 'unknown RNG vector schema')
    seeds = (0, 1, 2026092556, MASK64)
    orders = (2, 31, 127, 65587, 1439393, MASK64)
    require([(row['seed'], row['order']) for row in vectors['streams']] ==
            [(seed, order) for seed in seeds for order in (None, *orders)], 'RNG stream panel changed')
    count = 0
    for row in vectors['streams']:
        require(len(row['values']) == 257, 'RNG stream length changed')
        rng = StdRng08.seed_from_u64(row['seed'])
        expected = [(rng.next_u64() if row['order'] is None else rng.nonzero_below(row['order']))
                    for _ in row['values']]
        require(expected == row['values'], 'independent RNG differs from Rust stream')
        count += len(expected)
    require([(row['seed'], row['order']) for row in vectors['probes']] ==
            [(seed, order) for seed in seeds for order in orders], 'probe vector panel changed')
    for row in vectors['probes']:
        require(row['trials'] == [0, 1, 7, 63, 64, 65, 127, 128, 65535, MASK64],
                'probe trials changed')
        expected = [probe_scalar(row['seed'], t, row['order']) for t in row['trials']]
        require(expected == row['values'], 'independent probe differs from Rust')
        count += len(expected)
    return dict(schema_version=1, status='PASS', values_verified=count,
                streams=len(vectors['streams']), probe_panels=len(vectors['probes']),
                vector_sha256=canonical_digest(vectors), promotion_eligible=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rust-vectors', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify_rust_vectors(json.loads(args.rust_vectors.read_text())), sort_keys=True))
