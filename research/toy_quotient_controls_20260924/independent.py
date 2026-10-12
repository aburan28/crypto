"""Independent ordinary-polynomial F5B/Buchberger audit with Boolean equations."""
import hashlib
import inspect
import signal
import sys
import gc


def audit(n, packed, expected_roots, timeout=15, call_cap=100000):
    import sympy
    from sympy.polys import groebnertools as gt
    from sympy.polys.rings import ring
    from sympy.polys.domains import GF
    if sympy.__version__!='1.14.0':
        raise RuntimeError('frozen independent audit requires SymPy 1.14.0')
    if not 1<=n<=6:
        raise ValueError('toy audit is limited to six variables')
    R, *variables = ring(','.join(f'x{i}' for i in range(n)), GF(2), order='grevlex')
    masks = sorted(range(1<<n),key=lambda m:(m.bit_count(),m))
    F = [R.from_dict({tuple((mask>>j)&1 for j in range(n)): R.domain.one
                     for i,mask in enumerate(masks) if row>>i&1}) for row in packed if row]
    F = [x*x+x for x in variables] + F
    original = gt.f5_reduce
    original_profile = sys.getprofile()
    original_trace = sys.gettrace()
    gc_enabled = gc.isenabled()
    if original_profile is not None or original_trace is not None:
        raise RuntimeError('independent audit requires no existing profile hook')
    calls = 0
    class BudgetExceeded(Exception):
        pass
    def budget(frame,event,arg):
        nonlocal calls
        if event=='call':
            calls += 1
    def checkpoint(frame,event,arg):
        if frame.f_code is gt._f5b.__code__:
            if event=='line' and calls>call_cap:
                raise BudgetExceeded()
            return checkpoint
        return None
    degree_inputs, degree_outputs = [], []

    def degree(p):
        return max((sum(m) for m in p),default=0)

    def observed(f, basis):
        degree_inputs.append(degree(gt.Polyn(f)))
        out = original(f,basis)
        degree_outputs.append(degree(gt.Polyn(out)))
        return out

    def expired(signum, frame):
        raise TimeoutError('independent audit exceeded frozen timeout')

    old_handler = signal.signal(signal.SIGALRM, expired)
    gt.f5_reduce = observed
    signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        gt.f5_reduce = original
        buchberger = gt.groebner(F,R,method='buchberger')
        f5_status = 'VERIFIED'
        gt.f5_reduce = observed
        try:
            gc.disable()
            sys.setprofile(budget)
            sys.settrace(checkpoint)
            f5 = gt.groebner(F,R,method='f5b')
        except BudgetExceeded:
            f5_status = 'BUDGET_EXCEEDED'
            f5 = buchberger
        finally:
            sys.settrace(original_trace)
            sys.setprofile(original_profile)
            if gc_enabled:
                gc.enable()
            gt.f5_reduce = original
        if f5!=buchberger:
            raise AssertionError('independent reduced bases disagree')
        if any(f.rem(f5) for f in F):
            raise AssertionError('independent basis misses an input equation')
        for i,f in enumerate(f5):
            for g in f5[:i]:
                if gt.spoly(f,g,R).rem(f5):
                    raise AssertionError('independent basis fails Buchberger certificate')
        roots = [a for a in range(1<<n) if all(
            sum(all(not exponent or a>>j&1 for j,exponent in enumerate(m))
                for m,c in f.items() if int(c)&1)%2==0 for f in f5)]
        if roots!=expected_roots:
            raise AssertionError('independent roots differ')
        return dict(status='VERIFIED', sympy=sympy.__version__,
            f5b_status=f5_status, f5b_python_calls=calls, f5b_call_cap=call_cap,
            f5_reduce_source_sha256=hashlib.sha256(inspect.getsource(original).encode()).hexdigest(),
            f5b_source_sha256=hashlib.sha256(inspect.getsource(gt._f5b).encode()).hexdigest(),
            reduced_basis=[[[list(m),int(c)] for m,c in sorted(f.items())] for f in f5],
            roots=roots, f5b_signature_reduction_calls=len(degree_inputs),
            observed_f5b_reduction_input_degree=max(degree_inputs,default=None),
            observed_f5b_reduction_output_degree=max(degree_outputs,default=None),
            reduced_basis_max_degree=max(map(degree,f5),default=0),
            degree_scope='observed signature-reduction boundaries only; excludes preprocessing/intermediate reductions; not intrinsic regularity')
    except TimeoutError as error:
        return dict(status='TIMEOUT',message=str(error))
    finally:
        sys.settrace(original_trace)
        sys.setprofile(original_profile)
        if gc_enabled:
            gc.enable()
        signal.setitimer(signal.ITIMER_REAL,0)
        gt.f5_reduce = original
        signal.signal(signal.SIGALRM,old_handler)
