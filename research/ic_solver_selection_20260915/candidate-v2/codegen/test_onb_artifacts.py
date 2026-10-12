import json
import unittest
import indexcalc_e2e as e
from onb_artifacts import ArtifactCache


class MemoryRedis:
    def __init__(self):
        self.data = {}
    def eval(self, script, count, key, size):
        value = self.data.get(key)
        return value if value is None or len(value) <= size else None
    def set(self, key, value, ex):
        self.data[key] = value


class ArtifactTests(unittest.TestCase):
    def test_cold_warm_roundtrip_and_parameter_isolation(self):
        store = MemoryRedis()
        cold = ArtifactCache(store)
        a = e.setup(5,4,3,e.Ledger(),cold)
        warm = ArtifactCache(store)
        b = e.setup(5,4,3,e.Ledger(),warm)
        self.assertEqual(a[2:7],b[2:7])
        self.assertEqual(a[7].ops,b[7].ops)
        self.assertEqual(a[7].inputRef,b[7].inputRef)
        self.assertEqual(a[8],b[8])
        self.assertEqual(warm.stats['hits'],2)
        self.assertNotEqual(warm.key('x',[131,2]),warm.key('x',[131,4]))
        self.assertNotEqual(warm.key('x',[131,2]),warm.key('y',[131,2]))

    def test_corrupt_and_unavailable_recompute(self):
        store = MemoryRedis()
        cache = ArtifactCache(store)
        cache.getOrBuild('x',[],lambda:[1])
        key = cache.key('x',[])
        value = json.loads(store.data[key])
        value['payload'] = [2]
        store.data[key] = json.dumps(value).encode()
        self.assertEqual(cache.getOrBuild('x',[],lambda:[3]),[3])
        self.assertEqual(cache.stats['read_errors'],1)
        class Down:
            def eval(self,*args):
                raise OSError('offline')
            def set(self,*args,**kwargs):
                raise OSError('offline')
        cache = ArtifactCache(Down())
        self.assertEqual(cache.getOrBuild('x',[],lambda:[4]),[4])
        self.assertEqual(cache.stats['read_errors'],1)
        self.assertEqual(cache.stats['write_errors'],1)

    def test_oversized_not_published(self):
        store = MemoryRedis()
        cache = ArtifactCache(store,maxBytes=10)
        self.assertEqual(cache.getOrBuild('x',[],lambda:[1]),[1])
        self.assertFalse(store.data)
        self.assertEqual(cache.stats['oversized'],1)


if __name__ == '__main__':
    unittest.main()
