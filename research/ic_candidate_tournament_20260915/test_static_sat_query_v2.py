"""Preserve a conflict-capped CryptoMiniSat attempt as censored evidence."""
import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from static_sat_query_v2 import one_query


class StaticSatQueryV2Tests(unittest.TestCase):
    def test_conflict_cap_is_inconclusive_and_verification_is_timed(self):
        panel=dict(export_nonce=7,cms_conflict_budget=1000000,
                   export_timeout_seconds=60,cms_timeout_seconds=120)
        item=dict(trial=0,probe_scalar=3,point=[1,2])
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            def fake_meter(_,directory,name,__):
                if name=='export':
                    instance=directory/'instance'
                    instance.mkdir()
                    (instance/'manifest.json').write_text(json.dumps({}))
                    return dict(returncode=0,timed_out=False)
                (directory/'cms.stdout').write_text('s INDETERMINATE\n')
                return dict(returncode=15,timed_out=False,
                            command=['cms','--maxconfl','1000000'])
            started=time.monotonic_ns()
            with patch('static_sat_query_v2.meter',side_effect=fake_meter),\
                 patch('static_sat_query_v2.validate_export',return_value={}):
                row=one_query(panel,item,'exporter','cms',None,[],root)
            wall=time.monotonic_ns()-started
        self.assertEqual(row['status'],'CONFLICT_BUDGET_INCONCLUSIVE')
        self.assertIsNone(row['point_witness'])
        self.assertLessEqual(0,row['verification_wall_ns'])
        self.assertLessEqual(row['verification_wall_ns'],wall)


if __name__=='__main__':
    unittest.main()
