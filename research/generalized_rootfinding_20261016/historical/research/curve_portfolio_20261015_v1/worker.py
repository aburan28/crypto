import argparse,json,time
from .fields import BudgetedBinaryField
from .adapter import solve_chain,curve_capability
from rootfinder.fields import PrimeField,ExtensionField
from rootfinder.binaryfield import BinaryField
from rootfinder import polynomial as poly
from research.advance.general_chain import GeneralS3Fiber
from research.advance.test_general_chain import full_oracle,independent_s3
from lab.research.nist_fields import CertifiedPrimeField,preset

def run(started):
 rows=[]
 small=[(PrimeField(3),'prime',{'a2':0,'a4':1,'a6':1}),(PrimeField(5),'prime',{'a2':0,'a4':1,'a6':1}),(ExtensionField(3,2,(1,0,1)),'odd_extension',{'a2':3,'a4':1,'a6':4}),(BinaryField(3,11),'binary',{'a2':3,'b':3}),(BinaryField(4,19),'binary',{'a2':5,'b':7}),(BinaryField(3,11),'koblitz',{'a2':0,'b':1}),(BinaryField(3,11),'koblitz',{'a2':1,'b':1})]
 for f,family,curve in small:
  target=min(2,f.order-1);expected=full_oracle(f,(0,1),target,curve);r=solve_chain(f,(0,1),target,family=family,curve=curve,max_output=100000,max_seconds=20)
  assert r['complete'] and {tuple(x) for x in r['roots']}==expected and r['root_count_exact']==len(expected)
  rows.append(dict(family=family,field_bits=f.order.bit_length()-1,field_degree=f.degree,field_order=str(f.order),full_small_ordered_chain_root_count=len(expected),oracle='independent_D5_times_K3_exhaustion'))
 wide=[]
 for degree,mod in [(83,0x800000000200000000007),(163,(1<<163)|(1<<7)|(1<<6)|(1<<3)|1),(233,(1<<233)|(1<<74)|1),(283,(1<<283)|(1<<12)|(1<<7)|(1<<5)|1),(409,(1<<409)|(1<<87)|1),(571,(1<<571)|(1<<10)|(1<<5)|(1<<2)|1)]:
  f=BudgetedBinaryField(degree,mod,validation_seconds=15)
  for family,curve in [('koblitz',{'a2':1,'b':1}),('binary',{'a2':1<<(degree-3),'b':(1<<(degree-2))|1})]:
   cap=curve_capability(f,family,curve);fiber=GeneralS3Fiber(f,curve);x=(1<<(degree-2))|7;y=(1<<(degree-5))|11
   roots=fiber.roots(x,y);cs=fiber.coefficients(x,y);P=poly.trim(cs,f)
   frobenius=poly.powmod_poly((0,1),f.order,P,f);restricted=poly.gcd_poly(P,poly.sub(frobenius,(0,1),f),f)
   assert len(roots)==len(restricted)-1 and all(independent_s3(f,x,y,z,curve)==0 for z in roots)
   wide.append(dict(family=family,bits=degree,exact_quadratic_fiber_roots=len(roots),completeness='degree_gcd(P,X^q-X)',full_system_search=False))
 for name in ['P192','P224','P256','P384','P521']:
  rec=preset(name);f=CertifiedPrimeField(rec['p']);curve=rec['curve'];curve_capability(f,'prime',curve);fiber=GeneralS3Fiber(f,curve);x=(1<<(f.order.bit_length()-3))+7;y=x+11;roots=fiber.roots(x,y);P=poly.trim(fiber.coefficients(x,y),f);restricted=poly.gcd_poly(P,poly.sub(poly.powmod_poly((0,1),f.order,P,f),(0,1),f),f)
  assert len(roots)==len(restricted)-1 and all(independent_s3(f,x,y,z,curve)==0 for z in roots)
  wide.append(dict(family='prime',curve=name,bits=f.order.bit_length(),exact_quadratic_fiber_roots=len(roots),completeness='degree_gcd(P,X^q-X)',full_system_search=False))
 return dict(small_complete_allroot_oracles=rows,wide_complete_fiber_certificates=wide,arbitrary_hard_systems_solved=False,uniform_all_scales_speed_claim=False,all_index_calculus_models_supported=False,unsupported=['supersingular_binary_models','unconverted_general_Weierstrass_models','per_variable_domains_in_general_polynomial_adapter'],release_gate_claim=False)
def main():
 t=time.monotonic();p=argparse.ArgumentParser();p.add_argument('--experiment-id',required=True);a=p.parse_args()
 try:e=run(t);status='passed';summary='Prime, ordinary binary, Koblitz and odd-extension adapters matched small all-root oracles; wide finite fibers certified through571bits.'
 except Exception as x:e={'release_gate_claim':False,'hard_allroot_count_exact':None};status='inconclusive';summary=type(x).__name__+': '+str(x)
 r=dict(schema='scientific_result_v1',experiment_id=a.experiment_id,status=status,summary=summary,evidence=e,next_actions=['Promote each family/scale only after complete matching hard allroots benchmarks; finite-fiber and small-oracle coverage do not transfer the binary support optimization.']);s=json.dumps(r,sort_keys=True)
 if time.monotonic()-t>=75:r.update(status='inconclusive',summary='Late portfolio serialization.',evidence={'release_gate_claim':False,'hard_allroot_count_exact':None});s=json.dumps(r,sort_keys=True)
 print(s)
if __name__=='__main__':main()
