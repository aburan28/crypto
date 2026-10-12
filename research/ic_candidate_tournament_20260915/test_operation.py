import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from oracle import InvalidEvidence, verify
from tournament import candidates, decision, proposals_from_previous, snapshot_build, write


FIXTURE = {'cofactor':'4','curve_a':0,'degree':13,'generator':['4793','7108'],
    'group_order':'8012','irreducible':{'degree':13,'low_terms':[0,1,3,4]},
    'lambda':'89','subgroup_order':'2003','target_scalar_constructed':False,
    'target_seeds':[1009],'targets':[['2168','3628']]}


class CertificateTests(unittest.TestCase):
    def report(self):
        return {'schema_version':1,'status':'complete','mode':'rho','fixture':copy.deepcopy(FIXTURE),
                'solutions':[{'index':0,'recovered':'1830'}],'automorphism_order':26}

    def test_real_public_fixture_recovers_and_wrong_scalar_fails(self):
        r=self.report()
        self.assertEqual(verify(r,FIXTURE,expected_mode='rho')['verified_targets'],1)
        r['solutions'][0]['recovered']='1831'
        with self.assertRaises(InvalidEvidence):verify(r,FIXTURE,expected_mode='rho')

    def test_changed_target_and_duplicate_output_rejected(self):
        r=self.report();r['fixture']['targets'][0][0]='2169'
        with self.assertRaises(InvalidEvidence):verify(r,FIXTURE,expected_mode='rho')
        r=self.report();r['solutions']*=2
        with self.assertRaises(InvalidEvidence):verify(r,FIXTURE,expected_mode='rho')


class OperationTests(unittest.TestCase):
    def test_literal_native_dependencies_are_frozen(self):
        # Reproduces the observed include_str snapshot failure without compiling
        # the complete library in every test run.
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/'input';dest=Path(d)/'frozen'
            source.mkdir();dest.mkdir()
            (source/'src').mkdir();(source/'examples').mkdir();(source/'gpu').mkdir()
            (source/'Cargo.toml').write_text('[package]\nname="test"\n')
            (source/'Cargo.lock').write_text('lock')
            (source/'gpu/reduce.cu').write_text('original native source')
            (source/'src/lib.rs').write_text('const X:&str=include_str!("../gpu/reduce.cu");')
            (source/'examples/ic_tournament_worker.rs').write_text('fn main() {}')
            def fake_build(*args,**kwargs):
                output=dest/'build/release/examples/ic_tournament_worker'
                output.parent.mkdir(parents=True)
                output.write_bytes(b'unit-test executable placeholder')
            with patch('tournament.subprocess.run',fake_build):
                binary,manifest=snapshot_build(source,dest)
            self.assertTrue(binary.exists())
            self.assertIn('gpu/reduce.cu',manifest)
            (source/'gpu/reduce.cu').write_text('later source edit')
            self.assertEqual((dest/'source/gpu/reduce.cu').read_text(),'original native source')

    def test_successor_keeps_selected_baseline_and_avoids_previous_configs(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            write(root/'candidates.json',candidates())
            write(root/'decision.json',{'status':'promoted','winner':'batch16'})
            proposed,source=proposals_from_previous(root)
            self.assertEqual(proposed[0]['config']['batch_trials'],16)
            old={json.dumps(a['config'],sort_keys=True) for a in candidates()}
            self.assertTrue(all(json.dumps(a['config'],sort_keys=True) not in old for a in proposed[1:]))
            self.assertEqual(source,root/'source')

    def test_incomplete_incumbent_is_not_retained_as_a_verified_winner(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            write(root/'summaries/selection.json',{'provisional_challenger':None})
            for stage in ('confirmation','replay'):
                write(root/f'summaries/{stage}.json',{'comparisons':[],'rho_over_incumbent':None})
                write(root/f'runs/{stage}/case/incumbent/rep-0/receipt.json',
                      {'status':'TIMEOUT','total_operations':None})
            result=decision(root,{'repetitions':1,'unit':'test','evidence_scope':'test'},
                            {'confirmation':[{'id':'case'}],'replay':[{'id':'case'}]},[],save=False)
            self.assertIsNone(result['winner'])
            self.assertEqual(result['status'],'inconclusive')


if __name__=='__main__':unittest.main()
