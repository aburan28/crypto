import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/'ic_candidate_tournament_20260915'))
import agent_loop
import probe

class FreeModelTests(unittest.TestCase):
    def setUp(self):
        self.name=agent_loop.PREFERRED[0]
        self.catalog={'opencode':{'models':{self.name:{'tool_call':True,'cost':{'input':0,'output':0,'cache_read':0}}}}}
        self.active={'data':[{'id':self.name}]}
    def test_zero_prices_select(self):
        self.assertEqual(agent_loop.free_models(self.catalog,self.active),[self.name])
    def test_paid_cache_rejected(self):
        self.catalog['opencode']['models'][self.name]['cost']['cache_read']=0.01
        with self.assertRaises(RuntimeError):agent_loop.free_models(self.catalog,self.active)
    def test_missing_price_rejected(self):
        del self.catalog['opencode']['models'][self.name]['cost']['input']
        with self.assertRaises(RuntimeError):agent_loop.free_models(self.catalog,self.active)
    def test_inactive_rejected(self):
        with self.assertRaises(RuntimeError):agent_loop.free_models(self.catalog,{'data':[]})
    def test_deprecated_rejected(self):
        self.catalog['opencode']['models'][self.name]['status']='deprecated'
        with self.assertRaises(RuntimeError):agent_loop.free_models(self.catalog,self.active)
    def test_boolean_cost_is_not_a_price(self):
        self.catalog['opencode']['models'][self.name]['cost']['input']=False
        with self.assertRaises(RuntimeError):agent_loop.free_models(self.catalog,self.active)
    def test_single_provider_no_paid_fallback(self):
        c=agent_loop.config(self.name)
        self.assertEqual(c['enabled_providers'],['opencode'])
        self.assertEqual(c['provider']['opencode']['whitelist'],[self.name])
        self.assertEqual(c['permission']['task'],'deny')

class SubmissionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.base=self.root/'base';self.source=self.root/'source';self.built=self.root/'built'
        self.built.mkdir()
        for root in (self.base,self.source):
            (root/'src').mkdir(parents=True)
            (root/'src/binary_ecc.rs').write_text('pub fn add() {}\n')
            (root/'Cargo.toml').write_text('[package]\nname="test"\n')
        manifest={name:probe.t.digest(self.base/name) for name in ('src/binary_ecc.rs','Cargo.toml')}
        (self.built/'source-manifest.json').write_text(json.dumps(manifest))
        self.patch1=patch.object(probe,'BASE',self.base);self.patch2=patch.object(probe,'BASE_BUILT',self.built)
        self.patch1.start();self.patch2.start();self.addCleanup(self.patch1.stop);self.addCleanup(self.patch2.stop)
    def test_unchanged_valid(self):self.assertEqual(probe.validate(self.source,probe.CONFIG),[])
    def test_allowed_implementation_edit(self):
        (self.source/'src/binary_ecc.rs').write_text('pub fn add() { let _ = 1; }\n')
        self.assertEqual(probe.validate(self.source,probe.CONFIG),['src/binary_ecc.rs'])
    def test_frozen_build_rejected(self):
        (self.source/'Cargo.toml').write_text('changed')
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,probe.CONFIG)
    def test_added_counter_control_rejected(self):
        (self.source/'src/binary_ecc.rs').write_text('valgrind_counter_off();\n')
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,probe.CONFIG)
    def test_added_source_rejected(self):
        (self.source/'src/hidden.rs').write_text('bad')
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,probe.CONFIG)
    def test_symlink_rejected(self):
        p=self.source/'src/binary_ecc.rs';p.unlink();p.symlink_to(self.base/'src/binary_ecc.rs')
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,probe.CONFIG)
    def test_changed_summands_rejected(self):
        c={**probe.CONFIG,'summands':2}
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,c)
    def test_unknown_config_rejected(self):
        c={**probe.CONFIG,'free_oracle':True}
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,c)
    def test_trial_budget_rejected(self):
        c={**probe.CONFIG,'max_trials':65536}
        with self.assertRaises(probe.t.InvalidEvidence):probe.validate(self.source,c)

if __name__=='__main__':unittest.main()
