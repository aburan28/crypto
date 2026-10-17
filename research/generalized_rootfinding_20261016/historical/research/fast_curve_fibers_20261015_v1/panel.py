from rootfinder.fields import PrimeField,ExtensionField
from research.curve_portfolio_20261015_v1.fields import BudgetedBinaryField
from lab.research.nist_fields import CertifiedPrimeField,preset

def cases():
 for p in (3,5,7,11):yield f'prime-{p}',PrimeField(p),{'a2':0,'a4':1,'a6':1}
 for name in ('P192','P224','P256','P384','P521'):
  r=preset(name);yield name,CertifiedPrimeField(r['p']),r['curve']
 for n,m in [(2,7),(3,11),(4,19),(83,0x800000000200000000007),(163,(1<<163)|(1<<7)|(1<<6)|(1<<3)|1),(233,(1<<233)|(1<<74)|1),(283,(1<<283)|(1<<12)|(1<<7)|(1<<5)|1),(409,(1<<409)|(1<<87)|1),(571,(1<<571)|(1<<10)|(1<<5)|(1<<2)|1)]:
  f=BudgetedBinaryField(n,m)
  for a in (0,1):yield f'Koblitz-{n}-a{a}',f,{'a2':a,'b':1}
  yield f'binary-{n}',f,{'a2':1<<(n-1),'b':(1<<(n-2))|1}
 yield 'oddext-3-2',ExtensionField(3,2,(1,0,1)),{'a2':3,'a4':1,'a6':4}
 yield 'oddext-5-2',ExtensionField(5,2,(2,0,1)),{'a2':6,'a4':1,'a6':1}
