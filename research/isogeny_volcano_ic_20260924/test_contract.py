#!/usr/bin/env python3
import json,re
from pathlib import Path
R=Path(__file__).resolve().parent
c=json.loads((R/'experiments.json').read_text())
s=json.loads((R/'identity_schema.json').read_text())
assert c['schema_version']==1 and c['status']=='planned'
assert c['boundary_schema'].endswith('schema_version=2')
assert len(c['worker_lanes'])==8
assert sum(x['count'] for x in c['families'])>=80
assert c['common']['holdout_classes_per_field']>=2
assert len(c['common']['seeds'])>=5
assert 'ffd_or_degree_of_regularity' in c['required_measurements']
assert s['curve_id']['canonical'].startswith('ICV1:')
assert s['candidate_id']['canonical'].startswith('ICCAN1/')
assert any('j or order alone' in x for x in s['curve_id']['rules'])
print('PASS: campaign contract, >=80 experiment cells, 8 worker lanes, ICV1/ICCAN1 identity schema')
