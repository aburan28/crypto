#!/usr/bin/env python3
"""Tests for the index-calculus artifact cache.

Runs without `redis` or `boto3` installed: every tier is exercised through a
fake, and the two real tiers are tested against injected clients.

    python3 cache/test_cache.py
"""

import json
import logging
from pathlib import Path
import subprocess
import sys
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from cache import artifacts, keys, store  # noqa: E402
from cache.keys import CacheKey, CacheKeyError, curve_fingerprint  # noqa: E402
from cache.store import (  # noqa: E402
    L1, L2, L3, Codec, ObjectLRU, RedisTier, S3Tier, TieredStore,
    bytes_codec, json_codec,
)

# `docs/ic/binary-example.json`, inline so a schema change shows up as a
# failure here and not as a silently different fingerprint.
K0N9 = {
    "field": {"kind": "binary", "degree": 9, "polynomial_terms": [0, 1, 9]},
    "a": "0",
    "b": "1",
}
PRIME_TOY = {"field": {"kind": "prime", "modulus": "97"}, "a": "2", "b": "3"}


def setUpModule():
    """Degradation is the expected path here; keep its warnings off the report.

    `test_a_tier_failure_is_reported` re-enables them to check they happen.
    """
    logging.getLogger("cache.store").setLevel(logging.ERROR)


# -- fakes ----------------------------------------------------------------

class FakeByteTier:
    """A byte tier that counts calls and can be made to fail or lie."""

    def __init__(self, name, *, fail_get=False, fail_put=False,
                 max_value_bytes=None):
        self.name = name
        self.data = {}
        self.fail_get = fail_get
        self.fail_put = fail_put
        self.max_value_bytes = max_value_bytes
        self.gets = 0
        self.puts = 0
        self.ttls = []

    def _addr(self, key):
        return key.redis_key() if self.name == L2 else key.s3_key("indexcalc")

    def get(self, key):
        self.gets += 1
        if self.fail_get:
            raise ConnectionError(f"{self.name} is down")
        return self.data.get(self._addr(key))

    def put(self, key, blob, ttl):
        self.puts += 1
        if self.fail_put:
            raise ConnectionError(f"{self.name} is down")
        if self.max_value_bytes is not None and len(blob) > self.max_value_bytes:
            return False
        self.data[self._addr(key)] = blob
        self.ttls.append(ttl)
        return True


class Counter:
    """A compute function that records how often it actually ran."""

    def __init__(self, value):
        self.value = value
        self.calls = 0

    def __call__(self):
        self.calls += 1
        return self.value


def a_key(curve=K0N9, **params):
    params.setdefault("n", 5)
    return CacheKey.build(keys.NS_SUMMATION, curve_fingerprint(curve), **params)


# -- key canonicalisation -------------------------------------------------

class CanonicalIntTests(unittest.TestCase):
    def test_reads_the_repository_integer_encoding(self):
        self.assertEqual(keys.canonical_int("0x160"), 0x160)
        self.assertEqual(keys.canonical_int("120"), 120)
        self.assertEqual(keys.canonical_int(7), 7)

    def test_zero_has_one_canonical_form(self):
        # docs/ic writes coefficients as "0" and coordinates as "0x...";
        # if these disagreed the same curve would never share an entry.
        self.assertEqual(keys.canonical_int("0"), keys.canonical_int("0x0"))

    def test_rejects_nonsense(self):
        for bad in ("", "  ", "zz", "0xzz", None, 1.5, True, -3):
            with self.assertRaises(CacheKeyError):
                keys.canonical_int(bad)


class FieldCanonicalisationTests(unittest.TestCase):
    def test_polynomial_term_order_does_not_matter(self):
        # params.rs builds [low_terms..., degree]; docs/ic writes ascending.
        ascending = {"kind": "binary", "degree": 97, "polynomial_terms": [0, 6, 97]}
        params_rs = {"kind": "binary", "degree": 97, "polynomial_terms": [6, 0, 97]}
        self.assertEqual(keys.canonical_field(ascending),
                         keys.canonical_field(params_rs))

    def test_terms_must_bracket_the_polynomial(self):
        for terms in ([0, 1], [1, 9], []):
            with self.assertRaises(CacheKeyError):
                keys.canonical_field({"kind": "binary", "degree": 9,
                                      "polynomial_terms": terms})

    def test_abstract_representation_is_kept_verbatim_modulo_whitespace(self):
        onb = {"kind": "binary_abstract", "degree": 131,
               "representation": "permuted type-II optimal normal basis"}
        spaced = dict(onb, representation="permuted  type-II\toptimal normal basis\n")
        self.assertEqual(keys.canonical_field(onb), keys.canonical_field(spaced))

    def test_unknown_kinds_are_refused(self):
        with self.assertRaises(CacheKeyError):
            keys.canonical_field({"kind": "ternary", "degree": 3})


class CurveFingerprintTests(unittest.TestCase):
    def test_same_curve_two_spellings_one_fingerprint(self):
        other = {
            "field": {"kind": "binary", "degree": 9, "polynomial_terms": [9, 1, 0]},
            "a": 0,
            "b": "0x1",
        }
        self.assertEqual(curve_fingerprint(K0N9), curve_fingerprint(other))

    def test_different_irreducible_polynomial_is_a_different_curve(self):
        # The hazard this fingerprint exists for: F_2[z]/(z^9+z+1) and
        # F_2[z]/(z^9+z^4+1) are isomorphic fields and unrelated byte strings.
        other = dict(K0N9, field={"kind": "binary", "degree": 9,
                                  "polynomial_terms": [0, 4, 9]})
        self.assertNotEqual(curve_fingerprint(K0N9), curve_fingerprint(other))

    def test_polynomial_basis_and_normal_basis_do_not_collide(self):
        # `ecc2k130/codegen/decomp.py` works in the ONB while
        # `scripts/ecc2k130_point_decomposition.py` uses a reduction
        # polynomial; params.rs calls these `binary` and `binary_abstract`.
        poly = {"field": {"kind": "binary", "degree": 131,
                          "polynomial_terms": [0, 1, 2, 13, 131]},
                "a": "0", "b": "1"}
        onb = {"field": {"kind": "binary_abstract", "degree": 131,
                         "representation": "permuted type-II optimal normal basis"},
               "a": "0", "b": "1"}
        self.assertNotEqual(curve_fingerprint(poly), curve_fingerprint(onb))

    def test_different_abstract_representations_do_not_collide(self):
        one = {"field": {"kind": "binary_abstract", "degree": 131,
                         "representation": "permuted type-II optimal normal basis"},
               "a": "0", "b": "1"}
        two = dict(one, field=dict(one["field"], representation="Gaussian normal basis"))
        self.assertNotEqual(curve_fingerprint(one), curve_fingerprint(two))

    def test_curve_coefficients_are_in_the_fingerprint(self):
        self.assertNotEqual(curve_fingerprint(K0N9),
                            curve_fingerprint(dict(K0N9, a="0x1")))

    def test_prime_coefficients_are_reduced(self):
        self.assertEqual(curve_fingerprint(PRIME_TOY),
                         curve_fingerprint(dict(PRIME_TOY, a="99")))

    def test_unreduced_binary_coefficients_are_refused(self):
        # Silently masking would name the wrong bytes; failing says so.
        with self.assertRaises(CacheKeyError):
            curve_fingerprint(dict(K0N9, a=str(1 << 9)))

    def test_instance_data_is_excluded(self):
        # A summation polynomial depends on the curve alone.  Folding the
        # subgroup, generator, target or name into the key would turn a
        # fleet-wide artifact into a per-run one.
        instanced = dict(K0N9, name="synthetic-k0-n9", subgroup_order="127",
                         cofactor="4", generator={"x": "0x160", "y": "0x1d2"},
                         point={"x": "0x11e", "y": "0x177"},
                         fixture={"known_log": "120", "seed": 42},
                         schema_version=1)
        self.assertEqual(curve_fingerprint(K0N9), curve_fingerprint(instanced))

    def test_missing_pieces_are_refused(self):
        for drop in ("field", "a", "b"):
            partial = {k: v for k, v in K0N9.items() if k != drop}
            with self.assertRaises(CacheKeyError):
                curve_fingerprint(partial)

    def test_schema_version_changes_every_fingerprint(self):
        before = curve_fingerprint(K0N9)
        original = keys.SCHEMA_VERSION
        try:
            keys.SCHEMA_VERSION = original + 1
            self.assertNotEqual(before, curve_fingerprint(K0N9))
        finally:
            keys.SCHEMA_VERSION = original
        self.assertEqual(before, curve_fingerprint(K0N9))


class InputsFingerprintTests(unittest.TestCase):
    """The specialization digest must not alias two different inputs."""

    def test_a_mapping_does_not_collide_with_a_list_of_its_pairs(self):
        # Same canonical shape, different inputs.  Sharing an entry here
        # hands back the wrong artifact, not merely a slow one.
        self.assertNotEqual(keys.inputs_fingerprint({"a": 1}),
                            keys.inputs_fingerprint([("a", 1)]))

    def test_a_mapping_does_not_collide_with_a_set(self):
        self.assertNotEqual(keys.inputs_fingerprint({"a": 1}),
                            keys.inputs_fingerprint({("a", 1)}))

    def test_equal_inputs_agree_regardless_of_order(self):
        self.assertEqual(keys.inputs_fingerprint({"a": 1, "b": 2}),
                         keys.inputs_fingerprint({"b": 2, "a": 1}))

    def test_nested_order_still_matters(self):
        self.assertNotEqual(keys.inputs_fingerprint({"fb": [1, 2]}),
                            keys.inputs_fingerprint({"fb": [2, 1]}))

    def test_objects_keyed_by_memory_address_are_refused(self):
        # A plain object's repr is its address: unstable across processes,
        # and reused within one, so two distinct objects can share a key.
        class Opaque:
            pass

        with self.assertRaises(CacheKeyError):
            keys.inputs_fingerprint(Opaque())

    def test_a_type_with_its_own_repr_is_accepted(self):
        class Described:
            def __repr__(self):
                return "Described(deg=9)"

        self.assertEqual(keys.inputs_fingerprint(Described()),
                         keys.inputs_fingerprint(Described()))


class RepositoryDocumentTests(unittest.TestCase):
    """The committed ic documents must stay fingerprintable as written."""

    def test_docs_ic_examples_fingerprint(self):
        for name in ("binary-example.json", "prime-example.json"):
            path = REPO / "docs" / "ic" / name
            with self.subTest(document=name):
                document = json.loads(path.read_text())
                self.assertEqual(len(curve_fingerprint(document)),
                                 keys.FINGERPRINT_BITS // 4)

    def test_the_two_examples_are_different_curves(self):
        binary = json.loads((REPO / "docs/ic/binary-example.json").read_text())
        prime = json.loads((REPO / "docs/ic/prime-example.json").read_text())
        self.assertNotEqual(curve_fingerprint(binary), curve_fingerprint(prime))


# -- key rendering --------------------------------------------------------

class CacheKeyTests(unittest.TestCase):
    def test_parameter_order_does_not_matter(self):
        one = CacheKey.build("sumpoly", "ff", n=5, form="symmetric", codec="j")
        two = CacheKey.build("sumpoly", "ff", codec="j", form="symmetric", n=5)
        self.assertEqual(one, two)
        self.assertEqual(one.redis_key(), two.redis_key())

    def test_redis_key_carries_exactly_one_hash_tag(self):
        rendered = a_key().redis_key()
        self.assertEqual(rendered.count("{"), 1)
        self.assertEqual(rendered.count("}"), 1)
        tag = rendered[rendered.index("{") + 1:rendered.index("}")]
        self.assertEqual(tag, curve_fingerprint(K0N9))

    def test_one_curve_lands_on_one_shard(self):
        # Cluster mode hashes the tag, so every artifact for a curve shares a
        # slot and a warming worker does not fan out across the cluster.
        def tag(key):
            return key.redis_key().split("{")[1].split("}")[0]
        self.assertEqual(tag(a_key(n=5)), tag(a_key(n=6)))
        self.assertEqual(tag(a_key(n=5)),
                         tag(CacheKey.build(keys.NS_GROEBNER,
                                            curve_fingerprint(K0N9), m=3)))

    def test_s3_key_has_no_braces_and_honours_the_prefix(self):
        rendered = a_key().s3_key("indexcalc")
        self.assertNotIn("{", rendered)
        self.assertNotIn("}", rendered)
        self.assertTrue(rendered.startswith("indexcalc/v%d/" % keys.SCHEMA_VERSION))
        self.assertEqual(a_key().s3_key("/indexcalc/"), rendered)

    def test_schema_version_appears_in_both_renderings(self):
        self.assertIn("v%d" % keys.SCHEMA_VERSION, a_key().redis_key())
        self.assertIn("v%d" % keys.SCHEMA_VERSION, a_key().s3_key("p"))

    def test_segments_that_would_break_a_hash_tag_are_refused(self):
        for bad in ("a{b", "a}b", "a/b", "a b", ""):
            with self.assertRaises(CacheKeyError):
                CacheKey.build("sumpoly", "ff", n=bad)
            with self.assertRaises(CacheKeyError):
                CacheKey.build(bad, "ff", n=5)


# -- tiering --------------------------------------------------------------

class ReadThroughTests(unittest.TestCase):
    def setUp(self):
        self.l1 = ObjectLRU(max_entries=8)
        self.l2 = FakeByteTier(L2)
        self.l3 = FakeByteTier(L3)
        self.store = TieredStore(l1=self.l1, l2=self.l2, l3=self.l3)
        self.codec = json_codec()
        self.key = a_key()

    def get(self, compute):
        return self.store.get_or_compute(self.key, compute, self.codec)

    def test_a_miss_computes_once_and_fills_every_tier(self):
        compute = Counter({"terms": [1, 2, 3]})
        self.assertEqual(self.get(compute), {"terms": [1, 2, 3]})
        self.assertEqual(compute.calls, 1)
        self.assertEqual(len(self.l1), 1)
        self.assertEqual(len(self.l2.data), 1)
        self.assertEqual(len(self.l3.data), 1)

    def test_l1_hit_touches_no_network_tier(self):
        self.get(Counter({"terms": []}))
        gets_before = (self.l2.gets, self.l3.gets)
        compute = Counter({"never": True})
        self.assertEqual(self.get(compute), {"terms": []})
        self.assertEqual(compute.calls, 0)
        self.assertEqual((self.l2.gets, self.l3.gets), gets_before)

    def test_l2_hit_backfills_l1(self):
        self.get(Counter({"terms": [7]}))
        self.l1.clear()
        compute = Counter({"never": True})
        self.assertEqual(self.get(compute), {"terms": [7]})
        self.assertEqual(compute.calls, 0)
        self.assertEqual(self.store.stats()[L2]["hit"], 1)
        self.assertEqual(len(self.l1), 1)

    def test_l3_hit_backfills_l2_with_the_same_bytes(self):
        self.get(Counter({"terms": [7]}))
        stored = dict(self.l3.data)
        self.l1.clear()
        self.l2.data.clear()
        compute = Counter({"never": True})
        self.assertEqual(self.get(compute), {"terms": [7]})
        self.assertEqual(compute.calls, 0)
        self.assertEqual(self.store.stats()[L3]["hit"], 1)
        self.assertEqual(list(self.l2.data.values()), list(stored.values()))

    def test_redis_flush_refills_from_s3_rather_than_recomputing(self):
        # The load-bearing decision: Redis is a cache, S3 is the authority.
        compute = Counter({"terms": [7]})
        self.get(compute)
        self.l1.clear()
        self.l2.data.clear()
        self.get(compute)
        self.assertEqual(compute.calls, 1)

    def test_artifact_writes_carry_the_thirty_day_ttl(self):
        self.get(Counter({"terms": []}))
        self.assertEqual(self.l2.ttls, [store.DEFAULT_TTL_SECONDS])
        self.assertEqual(self.l3.ttls, [store.DEFAULT_TTL_SECONDS])

    def test_distinct_keys_do_not_share_entries(self):
        first = self.store.get_or_compute(a_key(n=5), Counter("s5"), self.codec)
        second = self.store.get_or_compute(a_key(n=6), Counter("s6"), self.codec)
        self.assertEqual((first, second), ("s5", "s6"))


class DegradationTests(unittest.TestCase):
    """Every tier failure degrades to recomputation, never to a wrong answer."""

    def setUp(self):
        self.codec = json_codec()
        self.key = a_key()

    def test_a_failing_read_still_returns_the_right_value(self):
        st = TieredStore(l1=ObjectLRU(4), l2=FakeByteTier(L2, fail_get=True),
                         l3=FakeByteTier(L3, fail_get=True))
        compute = Counter({"terms": [1]})
        self.assertEqual(st.get_or_compute(self.key, compute, self.codec),
                         {"terms": [1]})
        self.assertEqual(compute.calls, 1)
        self.assertEqual(st.stats()[L2]["error"], 1)
        self.assertEqual(st.stats()[L3]["error"], 1)

    def test_a_failing_write_still_returns_the_right_value(self):
        l2 = FakeByteTier(L2, fail_put=True)
        st = TieredStore(l1=ObjectLRU(4), l2=l2, l3=FakeByteTier(L3))
        self.assertEqual(st.get_or_compute(self.key, Counter("v"), self.codec), "v")
        self.assertEqual(st.stats()[L2]["error"], 1)
        self.assertEqual(st.stats()[L3]["write"], 1)

    def test_an_undecodable_value_is_a_miss(self):
        l2 = FakeByteTier(L2)
        st = TieredStore(l2=l2)
        st.get_or_compute(self.key, Counter({"terms": []}), self.codec)
        for addr in list(l2.data):
            l2.data[addr] = b"\xff not json at all"
        compute = Counter({"terms": [9]})
        self.assertEqual(st.get_or_compute(self.key, compute, self.codec),
                         {"terms": [9]})
        self.assertEqual(compute.calls, 1)
        self.assertEqual(st.stats()[L2]["error"], 1)

    def test_an_unencodable_value_is_returned_but_not_cached(self):
        l2 = FakeByteTier(L2)
        st = TieredStore(l1=ObjectLRU(4), l2=l2)
        sentinel = object()
        self.assertIs(st.get_or_compute(self.key, Counter(sentinel), self.codec),
                      sentinel)
        self.assertEqual(l2.data, {})
        self.assertEqual(st.stats()[L2]["error"], 1)

    def test_an_oversized_value_skips_redis_and_still_reaches_s3(self):
        l2 = FakeByteTier(L2, max_value_bytes=64)
        l3 = FakeByteTier(L3)
        st = TieredStore(l2=l2, l3=l3)
        big = "x" * 4096
        self.assertEqual(st.get_or_compute(self.key, Counter(big), self.codec), big)
        self.assertEqual(l2.data, {})
        self.assertEqual(st.stats()[L2]["skip"], 1)
        self.assertEqual(len(l3.data), 1)

    def test_a_dead_tier_stops_being_consulted(self):
        # A slow cache must never be slower than the work it replaces: 0.5 s
        # per relation attempt against a dead cluster is exactly that.
        l2 = FakeByteTier(L2, fail_get=True, fail_put=True)
        st = TieredStore(l2=l2)
        for index in range(12):
            st.get_or_compute(a_key(n=5, seq=index), Counter("v"), self.codec)
        self.assertLessEqual(l2.gets, store._BREAKER_THRESHOLD)
        self.assertGreater(st.stats()[L2]["skip"], 0)

    def test_the_breaker_reopens_after_its_cooldown(self):
        now = [1000.0]
        breaker = store._Breaker(threshold=2, cooldown=30.0, clock=lambda: now[0])
        self.assertTrue(breaker.closed())
        breaker.record_failure()
        breaker.record_failure()
        self.assertFalse(breaker.closed())
        now[0] += 31.0
        self.assertTrue(breaker.closed())
        breaker.record_success()
        self.assertTrue(breaker.closed())

    def test_a_tier_failure_is_reported(self):
        # Silent degradation is the other way to get this wrong: an operator
        # must be able to see that the cache stopped working.
        st = TieredStore(l2=FakeByteTier(L2, fail_get=True))
        with self.assertLogs("cache.store", level="WARNING") as captured:
            st.get_or_compute(self.key, Counter("v"), self.codec)
        self.assertTrue(any("l2" in line for line in captured.output), captured.output)

    def test_a_disabled_store_just_computes(self):
        st = TieredStore(enabled=False)
        compute = Counter("v")
        self.assertEqual(st.get_or_compute(self.key, compute, self.codec), "v")
        self.assertEqual(st.get_or_compute(self.key, compute, self.codec), "v")
        self.assertEqual(compute.calls, 2)

    def test_a_store_with_no_tiers_just_computes(self):
        st = TieredStore()
        compute = Counter("v")
        st.get_or_compute(self.key, compute, self.codec)
        st.get_or_compute(self.key, compute, self.codec)
        self.assertEqual(compute.calls, 2)


class ObjectLRUTests(unittest.TestCase):
    def test_eviction_is_bounded_and_least_recently_used(self):
        lru = ObjectLRU(max_entries=2)
        first, second, third = a_key(n=3), a_key(n=4), a_key(n=5)
        lru.put(first, "a")
        lru.put(second, "b")
        self.assertEqual(lru.get(first), "a")   # first is now most recent
        lru.put(third, "c")
        self.assertEqual(len(lru), 2)
        self.assertIsNone(lru.get(second))
        self.assertEqual(lru.get(first), "a")

    def test_zero_capacity_is_refused(self):
        with self.assertRaises(ValueError):
            ObjectLRU(max_entries=0)


# -- the real tiers, against injected clients -----------------------------

class FakeRedisClient:
    def __init__(self):
        self.data = {}
        self.sets = []

    def get(self, name):
        return self.data.get(name)

    def set(self, name, value, ex=None):
        self.data[name] = value
        self.sets.append((name, len(value), ex))


class NoSuchKey(Exception):
    def __init__(self):
        super().__init__("NoSuchKey")
        self.response = {"Error": {"Code": "NoSuchKey"}}


class Throttled(Exception):
    def __init__(self):
        super().__init__("SlowDown")
        self.response = {"Error": {"Code": "SlowDown"}}


class FakeS3Client:
    def __init__(self, *, error=None):
        self.data = {}
        self.error = error

    def get_object(self, Bucket, Key):  # noqa: N803 - boto3's own spelling
        if self.error is not None:
            raise self.error()
        if (Bucket, Key) not in self.data:
            raise NoSuchKey()
        return {"Body": _Body(self.data[(Bucket, Key)])}

    def put_object(self, Bucket, Key, Body):  # noqa: N803
        self.data[(Bucket, Key)] = Body


class _Body:
    def __init__(self, blob):
        self.blob = blob

    def read(self):
        return self.blob


class RedisTierTests(unittest.TestCase):
    def test_round_trip_and_ttl(self):
        client = FakeRedisClient()
        tier = RedisTier("redis://unused", client=client)
        key = a_key()
        self.assertIsNone(tier.get(key))
        self.assertTrue(tier.put(key, b"bytes", 42))
        self.assertEqual(tier.get(key), b"bytes")
        self.assertEqual(client.sets[-1][2], 42)
        self.assertIn("{", client.sets[-1][0])

    def test_the_value_cap_is_a_refusal_not_an_error(self):
        client = FakeRedisClient()
        tier = RedisTier("redis://unused", client=client, max_value_bytes=8)
        self.assertFalse(tier.put(a_key(), b"x" * 9, 42))
        self.assertEqual(client.data, {})

    def test_the_default_cap_is_eight_megabytes(self):
        self.assertEqual(store.MAX_REDIS_VALUE_BYTES, 8 * 1024 * 1024)
        self.assertEqual(store.REDIS_TIMEOUT_SECONDS, 0.5)


class S3TierTests(unittest.TestCase):
    def test_round_trip_under_the_prefix(self):
        client = FakeS3Client()
        tier = S3Tier("bucket", "indexcalc", client=client)
        key = a_key()
        self.assertIsNone(tier.get(key))
        tier.put(key, b"bytes", None)
        self.assertEqual(tier.get(key), b"bytes")
        self.assertIn(("bucket", key.s3_key("indexcalc")), client.data)

    def test_a_missing_object_is_a_miss(self):
        tier = S3Tier("bucket", "indexcalc", client=FakeS3Client(error=NoSuchKey))
        self.assertIsNone(tier.get(a_key()))

    def test_anything_else_is_a_tier_failure(self):
        # Throttling must reach the store so the breaker can open; swallowing
        # it as a miss would quietly turn S3 trouble into recomputation.
        tier = S3Tier("bucket", "indexcalc", client=FakeS3Client(error=Throttled))
        with self.assertRaises(Throttled):
            tier.get(a_key())


class FromEnvTests(unittest.TestCase):
    def test_absent_configuration_means_absent_tiers(self):
        st = TieredStore.from_env({})
        self.assertTrue(st.enabled)
        self.assertIsNotNone(st.tier(L1))
        self.assertIsNone(st.tier(L2))
        self.assertIsNone(st.tier(L3))

    def test_the_documented_settings_are_read(self):
        st = TieredStore.from_env({
            "INDEXCALC_REDIS_URL": "rediss://primary:6379",
            "INDEXCALC_S3_BUCKET": "artifacts",
            "INDEXCALC_S3_PREFIX": "indexcalc",
            "INDEXCALC_L1_ENTRIES": "7",
        })
        self.assertEqual(st.tier(L1).max_entries, 7)
        self.assertEqual(st.tier(L2).url, "rediss://primary:6379")
        self.assertEqual(st.tier(L3).bucket, "artifacts")
        self.assertEqual(st.tier(L3).prefix, "indexcalc")

    def test_clients_are_built_lazily(self):
        # `redis` and `boto3` are not installed in CI; constructing a store
        # from configuration must not need them.
        st = TieredStore.from_env({"INDEXCALC_REDIS_URL": "rediss://primary:6379",
                                   "INDEXCALC_S3_BUCKET": "artifacts"})
        self.assertIsNone(st.tier(L2)._client)
        self.assertIsNone(st.tier(L3)._client)

    def test_the_kill_switch(self):
        for value in ("0", "false", "no", ""):
            with self.subTest(value=value):
                st = TieredStore.from_env({"INDEXCALC_CACHE": value,
                                           "INDEXCALC_S3_BUCKET": "artifacts"})
                self.assertFalse(st.enabled)
                self.assertIsNone(st.tier(L3))

    def test_a_malformed_setting_falls_back_rather_than_raising(self):
        # A cache setting must not be able to kill a multi-hour run.
        st = TieredStore.from_env({"INDEXCALC_L1_ENTRIES": "sixty-four"})
        self.assertEqual(st.tier(L1).max_entries, 64)

    def test_l1_can_be_switched_off_alone(self):
        st = TieredStore.from_env({"INDEXCALC_L1_ENTRIES": "0"})
        self.assertIsNone(st.tier(L1))


# -- accessors ------------------------------------------------------------

class LazyImportTests(unittest.TestCase):
    """`redis` and `boto3` must stay optional, whatever the machine ships.

    Asserting they are *absent* is not a check this can own -- the GitHub
    runner image ships boto3, and what is installed will drift.  The property
    that actually matters is that importing this package, and building a
    store from configuration that names both tiers, pulls in neither module.
    A clean subprocess is the only place that can be observed honestly, since
    by the time a test runs some other import may already have loaded them.
    """

    # Note what is and is not asserted.  Importing the package and building
    # both network tiers from configuration must not load either module.
    # *Using* one legitimately does -- that is where the lazy import lives --
    # so the round trip below runs against an L1-only store instead.
    PROBE = """
import sys
sys.path.insert(0, {repo!r})
import cache

store = cache.TieredStore.from_env({{
    "INDEXCALC_REDIS_URL": "rediss://primary:6379",
    "INDEXCALC_S3_BUCKET": "artifacts",
}})
assert store.tier("l2") is not None and store.tier("l3") is not None
eager = [name for name in ("redis", "boto3") if name in sys.modules]
assert not eager, f"eagerly imported: {{eager}}"

local = cache.TieredStore.from_env({{}})
key = cache.CacheKey.build("sumpoly", cache.curve_fingerprint(
    {{"field": {{"kind": "binary", "degree": 9, "polynomial_terms": [0, 1, 9]}},
      "a": "0", "b": "1"}}), n=5)
assert local.get_or_compute(key, lambda: "v", cache.json_codec()) == "v"
assert local.get_or_compute(key, lambda: "w", cache.json_codec()) == "v"
print("clean")
"""

    def test_importing_the_package_imports_neither_optional_dependency(self):
        probe = self.PROBE.format(repo=str(REPO))
        done = subprocess.run([sys.executable, "-c", probe],
                              capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(done.stdout.strip(), "clean")


class AccessorTests(unittest.TestCase):
    def setUp(self):
        self.l2 = FakeByteTier(L2)
        self.l3 = FakeByteTier(L3)
        self.store = TieredStore(l1=ObjectLRU(8), l2=self.l2, l3=self.l3)
        self.codec = json_codec()

    def test_summation_polynomial_is_computed_once_per_curve(self):
        compute = Counter({"e": [1, 2, 3]})
        for _ in range(4):
            value = artifacts.summation_polynomial(self.store, K0N9, 5,
                                                   compute=compute, codec=self.codec)
        self.assertEqual(value, {"e": [1, 2, 3]})
        self.assertEqual(compute.calls, 1)

    def test_a_different_curve_is_a_different_artifact(self):
        other = dict(K0N9, field={"kind": "binary", "degree": 9,
                                  "polynomial_terms": [0, 4, 9]})
        first = artifacts.summation_polynomial(self.store, K0N9, 5,
                                               compute=Counter("A"), codec=self.codec)
        second = artifacts.summation_polynomial(self.store, other, 5,
                                                compute=Counter("B"), codec=self.codec)
        self.assertEqual((first, second), ("A", "B"))

    def test_the_symmetrised_and_expanded_forms_are_separate_entries(self):
        sym = artifacts.summation_polynomial(self.store, K0N9, 5, form=artifacts.SYMMETRIC,
                                             compute=Counter("sym"), codec=self.codec)
        exp = artifacts.summation_polynomial(self.store, K0N9, 5, form=artifacts.EXPANDED,
                                             compute=Counter("exp"), codec=self.codec)
        self.assertEqual((sym, exp), ("sym", "exp"))

    def test_the_codec_is_part_of_the_key(self):
        # Two serializers produce different bytes for one object, so they
        # must not read each other's entries.
        other = Codec(name="other", encode=self.codec.encode, decode=self.codec.decode)
        first = artifacts.summation_polynomial(self.store, K0N9, 5,
                                               compute=Counter("A"), codec=self.codec)
        second = artifacts.summation_polynomial(self.store, K0N9, 5,
                                                compute=Counter("B"), codec=other)
        self.assertEqual((first, second), ("A", "B"))

    def test_the_summation_index_is_validated(self):
        for bad in (1, 0, -1, "5", 2.0, True):
            with self.assertRaises(CacheKeyError):
                artifacts.summation_polynomial_key(K0N9, bad, "json")

    def test_the_form_is_validated(self):
        with self.assertRaises(CacheKeyError):
            artifacts.summation_polynomial_key(K0N9, 5, "json", form="compact")

    def test_groebner_structure_keys_separate_every_shape(self):
        base = dict(summands=3, factor_base_dim=12, monomial_order="grevlex",
                    algorithm="f4")
        key = artifacts.groebner_structure_key(K0N9, "json", **base)
        for field, value in (("summands", 4), ("factor_base_dim", 13),
                             ("monomial_order", "lex"), ("algorithm", "f5")):
            with self.subTest(changed=field):
                self.assertNotEqual(
                    key, artifacts.groebner_structure_key(K0N9, "json",
                                                          **dict(base, **{field: value})))

    def test_groebner_extra_options_reach_the_key(self):
        base = dict(summands=3, factor_base_dim=12, monomial_order="grevlex",
                    algorithm="f4")
        plain = artifacts.groebner_structure_key(K0N9, "json", **base)
        tuned = artifacts.groebner_structure_key(K0N9, "json", extra={"threads": 8}, **base)
        other = artifacts.groebner_structure_key(K0N9, "json", extra={"threads": 4}, **base)
        self.assertNotEqual(plain, tuned)
        self.assertNotEqual(tuned, other)

    def test_groebner_extra_options_are_order_insensitive(self):
        base = dict(summands=3, factor_base_dim=12, monomial_order="grevlex",
                    algorithm="f4")
        one = artifacts.groebner_structure_key(K0N9, "json", extra={"a": 1, "b": 2}, **base)
        two = artifacts.groebner_structure_key(K0N9, "json", extra={"b": 2, "a": 1}, **base)
        self.assertEqual(one, two)

    def test_groebner_shape_is_validated(self):
        base = dict(summands=3, factor_base_dim=12, monomial_order="grevlex",
                    algorithm="f4")
        for field, value in (("summands", 1), ("factor_base_dim", 0),
                             ("summands", "3")):
            with self.subTest(field=field, value=value):
                with self.assertRaises(CacheKeyError):
                    artifacts.groebner_structure_key(K0N9, "json",
                                                     **dict(base, **{field: value}))

    def test_specialization_is_redis_only_and_short_lived(self):
        compute = Counter([1, 2])
        for _ in range(3):
            artifacts.specialization(self.store, K0N9, inputs={"x": 5, "y": 9},
                                     compute=compute, codec=self.codec)
        self.assertEqual(compute.calls, 1)
        self.assertEqual(len(self.l2.data), 1)
        self.assertEqual(self.l3.data, {}, "T3 churn must not reach S3")
        self.assertEqual(self.l2.ttls, [store.SPECIALIZATION_TTL_SECONDS])
        self.assertLess(store.SPECIALIZATION_TTL_SECONDS, store.DEFAULT_TTL_SECONDS)

    def test_specialization_inputs_are_fingerprinted_order_insensitively(self):
        one = artifacts.specialization_key(K0N9, "json", {"x": 5, "y": 9})
        two = artifacts.specialization_key(K0N9, "json", {"y": 9, "x": 5})
        three = artifacts.specialization_key(K0N9, "json", {"x": 5, "y": 10})
        self.assertEqual(one, two)
        self.assertNotEqual(one, three)

    def test_specialization_distinguishes_nested_inputs(self):
        one = artifacts.specialization_key(K0N9, "json", {"fb": [1, 2], "t": 3})
        two = artifacts.specialization_key(K0N9, "json", {"fb": [2, 1], "t": 3})
        self.assertNotEqual(one, two)


class RealComputeTests(unittest.TestCase):
    """The injected-compute contract, against real code from this repository.

    `scripts/semaev_symbolic.py` builds S_4 as `{exponent tuple: 1}`.  Tuple
    keys are not JSON object keys, which is the point: the caller supplies a
    codec that knows its own representation, and this package never learns
    what a summation polynomial is.
    """

    @staticmethod
    def _codec():
        return Codec(
            name="semaev-exponents-v1",
            encode=lambda poly: json.dumps(sorted(list(e) for e in poly)).encode(),
            decode=lambda blob: {tuple(e): 1 for e in json.loads(blob.decode())},
        )

    def test_the_repository_s4_round_trips_through_every_tier(self):
        sys.path.insert(0, str(REPO / "scripts"))
        try:
            import semaev_symbolic
        finally:
            sys.path.pop(0)

        # This routine keeps a6 symbolic, so its result is shared more
        # narrowly than it strictly needs to be -- a caller-side choice about
        # what `compute` returns, not something the key schema can fix.
        curve = {"field": {"kind": "binary", "degree": 131,
                           "polynomial_terms": [0, 1, 2, 13, 131]},
                 "a": "0", "b": "1"}
        expected = semaev_symbolic.build_s4()
        compute = Counter(expected)
        l2, l3 = FakeByteTier(L2), FakeByteTier(L3)
        st = TieredStore(l1=ObjectLRU(4), l2=l2, l3=l3)

        first = artifacts.summation_polynomial(st, curve, 4, compute=compute,
                                               codec=self._codec())
        self.assertEqual(first, expected)

        # Cold process, flushed Redis: S3 must still spare the recomputation.
        cold = TieredStore(l1=ObjectLRU(4), l2=FakeByteTier(L2), l3=l3)
        never = Counter(None)
        second = artifacts.summation_polynomial(cold, curve, 4, compute=never,
                                                codec=self._codec())
        self.assertEqual(second, expected)
        self.assertEqual(compute.calls, 1)
        self.assertEqual(never.calls, 0)


class CodecTests(unittest.TestCase):
    def test_bytes_codec_round_trips(self):
        codec = bytes_codec()
        self.assertEqual(codec.decode(codec.encode(b"abc")), b"abc")

    def test_bytes_codec_refuses_non_bytes(self):
        with self.assertRaises(TypeError):
            bytes_codec().encode({"not": "bytes"})

    def test_json_codec_is_deterministic(self):
        codec = json_codec()
        self.assertEqual(codec.encode({"b": 1, "a": 2}), codec.encode({"a": 2, "b": 1}))


if __name__ == "__main__":
    unittest.main(verbosity=2)
