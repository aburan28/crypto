from research.advance.general_chain import GeneralS3Fiber,solve_general_chain
from research.advance.fiber_solver import solve_fiber_dag

def curve_capability(field,family,curve):
 if family not in ('prime','binary','koblitz','odd_extension'):raise ValueError('explicit supported curve family required')
 if family=='prime' and (field.characteristic==2 or field.degree!=1):raise ValueError('odd prime field required')
 if family=='odd_extension' and (field.characteristic==2 or field.degree<2):raise ValueError('odd extension required')
 if family in ('binary','koblitz') and field.characteristic!=2:raise ValueError('characteristic two required')
 f=GeneralS3Fiber(field,curve)
 if family=='koblitz' and (f.curve['b']!=1 or f.curve['a2'] not in (0,1)):raise ValueError('Koblitz form requires a2 in F2 and b=1')
 return dict(curve_family=family,field_characteristic=field.characteristic,field_degree=field.degree,field_order=str(field.order),curve=f.curve,all_field_sizes_performance_claim=False,koblitz_frobenius_quotient_enabled=False,contracts=['five-source-S3-full-field-auxiliary-chain','general-polynomial-system-common-explicit-domain'])

def solve_chain(field,domain,target,*,family,curve,**budgets):
 cap=curve_capability(field,family,curve);raw=solve_general_chain(field,domain,target,curve=curve,**budgets)
 complete=bool(raw['search_complete'] and raw['output_complete'])
 return dict(contract='all_ordered_original_chain_roots',complete=complete,status='complete' if complete else 'unresolved_or_truncated',root_count_exact=raw['solution_count_exact'],roots=raw['solutions'],capability=cap,raw=raw)

def solve_polynomials(field,equations,n_variables,domain,**budgets):
 raw=solve_fiber_dag(equations,n_variables,domain,field,**budgets);complete=bool(raw['count_complete'] and raw['output_complete'])
 return dict(contract='all_original_polynomial_roots_common_explicit_domain',complete=complete,status='complete' if complete else 'unresolved_or_truncated',root_count_exact=raw['solution_count_exact'],roots=raw['solutions'],raw=raw)
