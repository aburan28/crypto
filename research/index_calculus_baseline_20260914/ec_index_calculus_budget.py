#!/usr/bin/env python3
"""Auditable EC index-calculus budget calculator; estimates are not measurements.

Defaults describe an illustrative order approximately 2^60 subgroup, not a named
curve. Use --order for its actual prime order and --success for measured yield.
M counts signed allowed summands; B counts independent factor-base columns.
"""
import argparse
import json
import math

def budget(order, signed_base, columns, m, rho_orbit, rho_rate,
           success=None, usable=1.0, overhead=1.1, other_s=0.0,
           attempt_s=None, row_weight=None, matvec_factor=2.0):
    for name, value in [('order', order), ('signed_base', signed_base),
                        ('columns', columns), ('m', m), ('rho_orbit', rho_orbit)]:
        if type(value) is not int:
            raise ValueError(f'{name} must be an integer')
    for value in (rho_rate, success, usable, overhead, other_s, attempt_s,
                  row_weight, matvec_factor):
        if value is not None and not math.isfinite(value):
            raise ValueError('all supplied measurements and assumptions must be finite')
    if not(order>1 and signed_base>=1 and columns>=1 and m>=2 and rho_orbit>=1 and rho_rate>0):
        raise ValueError('positive order, base sizes, orbit and rate required; m >= 2')
    if not(0 < usable <= 1 and overhead>=1 and other_s>=0 and matvec_factor>0):
        raise ValueError('invalid usable, overhead, other_s or matvec_factor')
    if signed_base>order or columns>signed_base or rho_orbit>order:
        raise ValueError('base sizes or orbit inconsistent with group size')
    if attempt_s is not None and attempt_s<0:raise ValueError('negative attempt time')
    # Counting bound includes repeated summands. It does not presume uniform sums.
    tuples=math.comb(signed_base+m-1,m)
    p_upper=min(1.0,tuples/order)
    lam=math.exp(m*math.log(signed_base)-math.lgamma(m+1)-math.log(order))
    heuristic_p=-math.expm1(-lam)
    p=heuristic_p if success is None else success
    if not(0<=p<=1):raise ValueError('success must be in [0, 1]')
    rows=math.ceil(overhead*columns)
    rho_steps=math.sqrt(math.pi*order/(2*rho_orbit))
    rho_s=rho_steps/rho_rate
    attempts=rows/(p*usable) if p else None
    available=max(0.0,rho_s-other_s)
    cap=available/attempts if attempts else None
    w=m if row_weight is None else row_weight
    if w<=0:raise ValueError('row weight must be positive')
    # Stylized O(B) sparse matrix-vector sequence, not an exact Wiedemann count.
    matrix_mac=matvec_factor*columns*rows*w
    result={
        'evidence':'calculated_model_not_measured',
        'order':order,'signed_factor_base_M':signed_base,'independent_columns_B':columns,
        'm':m,'rho_orbit':rho_orbit,'rho_rate_steps_per_s':rho_rate,
        'lambda_heuristic':lam,'p_heuristic':heuristic_p,
        'p_counting_upper_bound_uniform_target':p_upper,
        'p_used':p,'p_source':'heuristic' if success is None else 'user_supplied',
        'usable_fraction':usable,'rows_target_not_rank_guarantee':rows,
        'expected_attempts':attempts,'rho_expected_steps':rho_steps,
        'attempt_count_floor_first_relation_mode':rows/p_upper,
        'count_floor_scope':'uniform independent targets; at most one fresh row per attempt',
        'rho_expected_s':rho_s,'other_attack_cost_s':other_s,
        'other_costs_already_exhaust_rho_budget':other_s>=rho_s,
        'full_dlp_S':None,'measured_rho_ratio':None,
        'attempt_budget_s':cap,'attempt_budget_rho_step_equivalents':cap*rho_rate if cap is not None else None,
        'assumed_relation_matrix_row_weight':w,'assumed_matvec_factor':matvec_factor,
        'stylized_relation_matrix_mod_r_coefficient_accumulations':matrix_mac,
        'matrix_operations_are_not_GF2_operations':True,
        'warnings':['Poisson success is a heuristic; record empirical success and rank.',
                    'Matrix update estimate is not automatically included in other_s; supply measured matrix time.',
                    'A coefficient accumulation can be an addition for coefficients +/-1; it is not always a full modular multiplication.',
                    'All rho rates must match the curve, quotient and hardware being compared.']
    }
    if attempt_s is not None:
        result['modeled_full_attack_s']=other_s+attempts*attempt_s if attempts else None
        result['modeled_rho_over_ic']=rho_s/result['modeled_full_attack_s'] if result['modeled_full_attack_s'] else None
    return result

def monomials(v,d,boolean=True):
    if v<1 or d<0:raise ValueError('invalid variable count or degree')
    c=sum(math.comb(v,i) for i in range(min(d,v)+1)) if boolean else math.comb(v+d,d)
    return {'variables':v,'degree':d,'boolean':boolean,'column_universe':c,
            'hypothetical_square_bit_matrix_GiB':c*c/(8*2**30) if boolean else None,
            'actual_F4_memory_prediction':False}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--order',type=int,default=2**60)
    p.add_argument('--signed-base',type=int,default=2**20)
    p.add_argument('--columns',type=int,default=2**19)
    p.add_argument('--m',type=int,default=3)
    p.add_argument('--rho-orbit',type=int,default=2)
    p.add_argument('--rho-rate',type=float,default=1e6)
    p.add_argument('--success',type=float)
    p.add_argument('--usable',type=float,default=1.0)
    p.add_argument('--overhead',type=float,default=1.1)
    p.add_argument('--other-s',type=float,default=0.0)
    p.add_argument('--attempt-s',type=float)
    a=p.parse_args()
    try:
        r=budget(a.order,a.signed_base,a.columns,a.m,a.rho_orbit,a.rho_rate,a.success,a.usable,a.overhead,a.other_s,a.attempt_s)
    except (ValueError, OverflowError) as exc:
        p.error(str(exc))
    r['macaulay_examples']=[monomials(30,d) for d in [4,5,6]]
    print(json.dumps(r,indent=2,allow_nan=False))

if __name__=='__main__':main()
