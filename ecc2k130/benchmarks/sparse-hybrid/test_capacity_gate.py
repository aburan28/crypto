"""Reproduction must reject changed source, instructions or resources before GPU timing."""
from pathlib import Path
import json,re,tempfile,unittest
import capacity_gate

ROOT=Path(__file__).resolve().parent

class CapacityGateTests(unittest.TestCase):
    def check(self,mutation=None):
        actual=json.loads((ROOT/'evidence/capacity-compile-result.json').read_text())
        if mutation:mutation(actual)
        with tempfile.TemporaryDirectory() as directory:
            raw=Path(directory)/'compile.json';review=Path(directory)/'review.json'
            raw.write_text(json.dumps(actual,indent=2)+'\n')
            return capacity_gate.bind_review(raw,review)

    def test_exact_reproduction(self):
        self.assertTrue(self.check()['valid'])

    def test_changed_instruction(self):
        def mutate(actual):
            row=next(x for x in actual['commands'] if x['label']=='SASS')
            row['output']=re.sub(r'0x[0-9a-f]{16}',lambda m:'0x'+format(int(m[0],16)^1,'016x'),row['output'],count=1)
        with self.assertRaises(AssertionError):self.check(mutate)

    def test_changed_resources(self):
        def mutate(actual):
            row=next(x for x in actual['commands'] if x['label']=='resources')
            row['output']=row['output'].replace('REG:40','REG:41',1)
        with self.assertRaises(AssertionError):self.check(mutate)

    def test_changed_source(self):
        def mutate(actual):actual['sourceHashes']['bench.cu']='0'*64
        with self.assertRaises(AssertionError):self.check(mutate)

if __name__=='__main__':unittest.main()
