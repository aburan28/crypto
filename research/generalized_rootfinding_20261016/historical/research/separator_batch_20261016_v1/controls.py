import json
from .core import BatchSeparator
from rootfinder.fields import PrimeField,ExtensionField
from rootfinder.binaryfield import BinaryField
from research.advance.test_general_chain import full_oracle
from research.advance.general_chain import _Unresolved
cases=0
for f,c in [(PrimeField(3),{'a2':0,'a4':1,'a6':1}),(PrimeField(5),{'a2':0,'a4':1,'a6':1}),(ExtensionField(3,2,(1,0,1)),{'a2':3,'a4':1,'a6':4}),(BinaryField(3,11),{'a2':0,'b':1}),(BinaryField(4,19),{'a2':5,'b':7})]:
 batch=BatchSeparator(f,(0,1),c)
 for t in (0,1,2):
  result=batch.solve(t);expected=full_oracle(f,(0,1),t,c)
  assert result['search_complete'] and result['output_complete'] and {tuple(r) for r in result['solutions']}==expected and result['solution_count_exact']==len(expected);cases+=1
 truncated=batch.solve(0,max_output=0)
 assert not truncated['output_complete'] or truncated['solution_count_exact']==0
 try:BatchSeparator(f,(0,1),c,max_left_records=0)
 except _Unresolved:pass
 else:raise AssertionError('unboundedleftaccepted')
 empty=BatchSeparator(f,(),c);assert empty.solve(0)['solution_count_exact']==0
print(json.dumps({'status':'passed','independent_exhaustive_D5K3_cases':cases,'empty_and_truncated_and_preflight_refusal':True,'actualwide_runs':0,'external_independent_review':False}))
