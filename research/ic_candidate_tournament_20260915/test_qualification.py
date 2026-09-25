import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from oracle import InvalidEvidence
from qualification import reference_report
from tournament import (comparison, gate, is_rho, qualification_references,
                        run_campaign, stage_arms)
from test_tournament import rows


class ReferenceQualificationTests(unittest.TestCase):
    def test_widths_cover_each_source_once_and_preserve_source_binding(self):
        base = dict(id='incumbent', config={'summands':3}, source_manifest_sha256='a', binary_relative='worker')
        same = dict(base, id='same_source')
        other = dict(base, id='other', source_manifest_sha256='b', binary_relative='other/worker')
        refs = qualification_references([base,same,other])
        self.assertEqual(len(refs),6)
        self.assertTrue(all(is_rho(r) for r in refs))
        self.assertEqual({r['config']['rho_parallel_walks'] for r in refs},{1,8,32})
        self.assertEqual({r['binary_relative'] for r in refs[3:]},{'other/worker'})
        self.assertNotIn('rho_parallel_walks',base['config'])

    def test_named_rho_arms_do_not_require_an_ic_factor_base(self):
        data = rows(.7)
        for row in data:
            if row['arm']=='candidate':
                row.update(arm='rho_scaled_8', mode='rho', certificate={})
        self.assertTrue(comparison(data,'rho_scaled_8',draws=20)['eligible'])

    def test_reference_selection_never_promotes_even_with_qualification_field(self):
        contract=dict(purpose='reference-qualification', reference_qualification={'complete':True},
                      confirmation_ratio=.8,max_cell_ratio=1.1)
        self.assertFalse(gate(comparison(rows(.1),'candidate',draws=20),contract))

    def test_development_retains_failed_references_and_every_width(self):
        arms=[dict(id='incumbent'),dict(id='challenger')]
        refs=[dict(id='rho_incumbent_1',kind='rho-reference')]
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'summaries').mkdir()
            (root/'contract.json').write_text(json.dumps(dict(purpose='reference-qualification',reference_arms=refs)))
            (root/'summaries/smoke.json').write_text(json.dumps(dict(failures=[dict(arm='challenger')])))
            self.assertEqual(stage_arms(root,'development',arms),arms+refs)

    def test_online_and_cold_leaders_remain_separate_and_smoke_failure_disqualifies(self):
        cases=[dict(id='n13a0-000',cell='n13a0',fixture={'subgroup_order':'2003'})]
        ic=[dict(id=name,source_manifest_sha256=name,config={'summands':3})
            for name in ('incumbent','cold','online')]
        refs=[dict(id='rho_incumbent_1',kind='rho-reference',source_manifest_sha256='rho',
                   config={'summands':3,'rho_parallel_walks':1})]
        arms=ic+refs
        values={'incumbent':(1000,1000,100),'cold':(500,500,90),
                'online':(900,900,50),'rho_incumbent_1':(2000,2000,200)}
        data=[]
        for arm in arms:
            instructions,cold,online=values[arm['id']]
            for rep in range(3):
                data.append(dict(arm=arm['id'],case='n13a0-000',cell='n13a0',repetition=rep,
                    case_sha256='fixed',status='VERIFIED',mode='rho' if is_rho(arm) else 'ic',
                    total_operations=instructions,native_process={'process_wall_seconds':cold/1e9},
                    certificate={'factor_base_sha256':'fixed'},measurement={'native_timing':{
                        'cold':{'wall_ns':cold},'online':{'wall_ns':online}}}))
        contract=dict(purpose='reference-qualification',reference_arms=refs,repetitions=3,
                      comparison_kind='fixed-support',bootstrap_draws=20,cells=['n13a0'],
                      stages=['aa','smoke','development'])
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'summaries').mkdir()
            (root/'summaries/aa.json').write_text('{"passed":true}')
            for rep in range(3):
                for mode in ('native','profile'):
                    path=root/'runs/development/n13a0-000/rho_incumbent_1'/f'rep-{rep}'/mode
                    path.mkdir(parents=True)
                    (path/'stdout.json').write_text('{"solutions":[{"effective_walks":1}]}')
            fixtures={'development':cases,'smoke':cases}
            with patch('tournament.load_stage',return_value=data):
                report=reference_report(root,contract,fixtures,ic)
            self.assertEqual(report['selected_ic_cold'],'cold')
            self.assertEqual(report['selected_ic_online'],'online')
            self.assertFalse(report['promotion_eligible'])
            smoke=copy.deepcopy(data)
            next(row for row in smoke if row['arm']=='cold')['status']='TIMEOUT'
            with patch('tournament.load_stage',side_effect=[data,smoke]):
                report=reference_report(root,contract,fixtures,ic)
            self.assertEqual(report['selected_ic_cold'],'online')
            failed=next(row for row in report['table'] if row['alias']=='cold')
            self.assertFalse(failed['qualified'])
            self.assertEqual(len(failed['smoke_failures']),1)

    def test_confirmation_is_not_an_executable_qualification_stage(self):
        with tempfile.TemporaryDirectory() as directory:
            contract=dict(stages=['aa','smoke','development'],purpose='reference-qualification')
            args=type('Args',(),dict(round=Path(directory),stage='confirmation'))()
            with patch('tournament.frozen_inputs',return_value=(contract,{},[])):
                with self.assertRaisesRegex(InvalidEvidence,'stage is not in this frozen campaign'):
                    run_campaign(args)


if __name__=='__main__':
    unittest.main()
