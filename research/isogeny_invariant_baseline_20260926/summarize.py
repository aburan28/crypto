"""Audit exact incidence, emit a deterministic stage-only summary (no DLP speed claim)."""
import argparse,gzip,hashlib,json
from pathlib import Path
P=Path(__file__).parent
EXPECTED=[
 ('y^2=x^3+1x+1',64,[(5,7),(21,63),(30,63)]),
 ('y^2=x^3+1x+0',64,[(5,7),(21,63),(3,3)]),
 ('y^2=x^3+1x+2',64,[(5,7),(21,63),(30,63)]),
 ('y^2=x^3+2x+0',64,[(51,63),(3,3),(21,63)]),
 ('y^2=x^3+2x+1',91,[(30,90),(30,90),(30,90)]),
 ('y^2=x^3+2x+2',91,[(30,90),(30,90),(30,90)]),
]
def audit(raw):
 cases=[raw['tiny']]+raw['other_five_curves'];assert len(cases)==6
 assert len({c['curve'] for c in cases})==6
 rows=[]
 for c,(curve,order,expected) in zip(cases,EXPECTED):
  assert c['curve']==curve and c['curve_order']==order
  assert len(c['points'])==order and c['points'][0] is None
  assert c['affine_target_count']==order-1 and len(c['all_40_hyperplanes'])==40
  linear=[v for v in c['all_40_hyperplanes'] if v['frobenius_invariant']]
  assert len(linear)==2
  assert [v['functional'] for v in linear]==[[0,1,1,2],[1,1,1,1]]
  variants=[c['fraction']]+linear
  assert [(v['affine_base_points'],v['covered_affine_targets']) for v in variants]==expected
  for v in [c['fraction']]+c['all_40_hyperplanes']:
   assert v['field_set_size']==27
   counts=v['counts_in_point_order'];distinct=v['distinct_x_counts_in_point_order']
   assert len(counts)==len(distinct)==order
   assert sum(counts)==v['affine_base_points']**2
   assert all(0<=b<=a for a,b in zip(counts,distinct))
   assert sum(x>0 for x in counts[1:])==v['covered_affine_targets']
   assert sum(x>0 for x in distinct[1:])==v['covered_affine_targets_distinct_x']
   assert sum(int(k)*n for k,n in v['ordered_pair_count_histogram_affine'].items())==sum(counts[1:])
   if v['frobenius_invariant']:
    assert sum(int(k)*n for k,n in v['field_orbits'].items())==27
    assert sum(int(k)*n for k,n in v['point_orbits'].items())==v['affine_base_points']
  assert c['f3_equivalence_checks']==len(_ys_from_points(c))**3
  for v in variants:
   rows.append({'curve':curve,'variant':v['name'],'field_set_size':27,'affine_base_points':v['affine_base_points'],'covered_affine_targets':v['covered_affine_targets'],'affine_targets':order-1,'distinct_x_covered_targets':v['covered_affine_targets_distinct_x'],'stage':'two_summand_incidence','total_common_operations_over_sqrt_subgroup_order':None,'matched_rho_operations':None,'speedup':None,'verified':True})
 p=raw['published_13_7']
 assert p['q']==13 and p['n']==7 and p['set_size']==2197
 assert p['frobenius_orbits']=={'1':13,'7':312}
 assert p['representation_multiplicity_histogram']=={'1':2184,'14':13}
 assert sum(c['f3_equivalence_checks'] for c in cases)==325998
 return {'schema_version':1,'evidence_type':'exact_small_field_stage_diagnostic','status':'verified','operation_unit':'affine_point_count_and_affine_target_coverage','end_to_end_cost':None,'runtime_comparison':None,'baseline_comparator':'all_40_three_dimensional_F3_subspaces_with_two_invariant','cases':len(cases),'summation_lifting_checks':325998,'published_construction':{'q':13,'degree':7,'set_size':2197,'orbits':p['frobenius_orbits']},'source_sha256':hashlib.sha256((P/'benchmark.py').read_bytes()).hexdigest(),'rows':rows}
def _ys_from_points(c):
 return {x for p in c['points'] if p is not None for x in [p[0]]}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('raw');ap.add_argument('--output',required=True);args=ap.parse_args()
 raw_path=Path(args.raw)
 raw_text=gzip.open(raw_path,'rt').read() if raw_path.suffix=='.gz' else raw_path.read_text()
 result=audit(json.loads(raw_text))
 with open(args.output,'x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(f"Audited {result['cases']} curves, {len(result['rows'])} rows and {result['summation_lifting_checks']} lifting checks; S and speedup unset")
