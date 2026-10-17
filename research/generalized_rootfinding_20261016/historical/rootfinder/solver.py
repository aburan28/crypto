"""A bounded exact solver for the generalized finite-field root problem."""

from .fields import make_field
from .composition import composition_roots
from .domain import linearized_roots, polynomial_domain
from .search import solve_on_domain, _components
from .native_boolean import solve_boolean_mq, supports_boolean_mq


DEFAULTS = dict(max_work=1_000_000, max_table=100_000, max_solutions=10_000,
                max_domain=4096, max_degree=256, max_trials=2048, seed=0)


def element(value, field):
    if isinstance(value, dict):
        if set(value) != {"scalar"} or type(value["scalar"]) is not int:
            raise ValueError("A scalar coefficient must be {scalar: integer}")
        return value["scalar"] % field.characteristic
    if isinstance(value, str):
        value = int(value, 16) if value.lower().startswith("0x") else int(value, 10)
    if type(value) is not int:
        raise ValueError("Field elements must be integer or hexadecimal polynomial-basis encodings")
    field.validate(value)
    return value


def coefficients(values, field):
    if not isinstance(values, list) or not values:
        raise ValueError("Coefficients must be a nonempty ascending-degree list")
    result = [element(x, field) for x in values]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def normalize_equations(raw, n, field):
    if not isinstance(raw, list) or not raw:
        raise ValueError("Supply at least one polynomial equation")
    equations = []
    for equation in raw:
        if not isinstance(equation, list):
            raise ValueError("An equation must be a list of sparse terms")
        merged = {}
        for term in equation:
            powers = term["powers"]
            if (not isinstance(powers, list) or len(powers) != n or
                    any(type(e) is not int or e < 0 or e > 1_000_000 for e in powers)):
                raise ValueError("Each powers list must match variables and contain exponents from 0 to 1000000")
            powers = tuple(powers)
            c = element(term["coefficient"], field)
            merged[powers] = field.add(merged.get(powers, 0), c)
        equations.append([(c, exps) for exps, c in sorted(merged.items()) if c])
    return equations


def solve_problem(problem, **overrides):
    if not isinstance(problem, dict):
        raise ValueError("The problem must be a JSON object")
    unknown = set(overrides) - set(DEFAULTS)
    if unknown:
        raise ValueError("Unknown limit option(s): " + ", ".join(sorted(unknown)))
    limits = {**DEFAULTS, **overrides}
    for key in DEFAULTS:
        if type(limits[key]) is not int or limits[key] < 0:
            raise ValueError(f"{key} must be a nonnegative integer")
    field_spec = problem.get("field", {})
    if not isinstance(field_spec, dict):
        raise ValueError("field must be a JSON object")
    field = make_field(field_spec)
    variables = problem["variables"]
    if (not isinstance(variables, list) or not variables or
            any(not isinstance(x, str) or not x for x in variables) or
            len(set(variables)) != len(variables)):
        raise ValueError("variables must be a nonempty list of distinct names")
    n = len(variables)
    raw = problem.get("equations")
    if raw is None:
        raw = [problem["f"]]
    equations = normalize_equations(raw, n, field)
    l_spec = problem["L"]
    if not isinstance(l_spec, dict):
        raise ValueError("L must be a JSON object")
    representations = sum(name in l_spec for name in ("linearized", "composition", "coefficients"))
    if representations != 1:
        raise ValueError("L must specify exactly one of coefficients, linearized, or composition")
    if "composition" in l_spec:
        if "denominator" in l_spec:
            raise ValueError("Composition denominators belong to individual maps")
        raw_maps = l_spec["composition"]
        if not isinstance(raw_maps, list) or not raw_maps:
            raise ValueError("L.composition must be a nonempty inner-to-outer map list")
        maps = []
        rational = False
        for raw_map in raw_maps:
            if not isinstance(raw_map, dict):
                raise ValueError("Each composition map must be a JSON object")
            numerator = coefficients(raw_map["coefficients"], field)
            denominator = coefficients(raw_map.get("denominator", [1]), field)
            rational |= len(denominator) > 1
            maps.append({"coefficients": numerator, "denominator": denominator})
        if rational and l_spec.get("semantics") != "partial_circuit":
            raise ValueError("Rational compositions require explicit L.semantics='partial_circuit' to define intermediate poles")
        if "semantics" in l_spec and l_spec["semantics"] != "partial_circuit":
            raise ValueError("Unsupported composition semantics")
        domain = composition_roots(maps, field, max_domain=limits["max_domain"],
                                   max_degree=limits["max_degree"], max_work=limits["max_work"],
                                   max_trials=limits["max_trials"], seed=limits["seed"])
    elif "linearized" in l_spec:
        if "denominator" in l_spec or "coefficients" in l_spec:
            raise ValueError("A linearized L must not also specify dense coefficients or a denominator")
        linear = l_spec["linearized"]
        if not isinstance(linear, dict):
            raise ValueError("L.linearized must be a JSON object")
        terms = [(t["power"], element(t["coefficient"], field)) for t in linear["terms"]]
        constant = element(linear.get("constant", 0), field)
        domain = linearized_roots(terms, constant, field, limits["max_domain"])
    else:
        numerator = coefficients(l_spec["coefficients"], field)
        denominator = coefficients(l_spec.get("denominator", [1]), field)
        domain = polynomial_domain(numerator, denominator, field,
                                  max_domain=limits["max_domain"], max_degree=limits["max_degree"],
                                  seed=limits["seed"], max_trials=limits["max_trials"])
    report = {
        "field": {"kind": {"PrimeField": "prime", "BinaryField": "binary", "ExtensionField": "extension"}.get(type(field).__name__, "finite"), "degree": field.degree,
                  "characteristic": str(field.characteristic), "order": str(field.order),
                  "modulus": hex(field.modulus) if isinstance(field.modulus, int) else field.modulus,
                  "field_verified": True,
                  "element_encoding": "canonical residue" if field.degree == 1 else "ascending base-characteristic polynomial-basis digits encoded as an integer"},
        "variables": variables, "limits": limits, "domain": domain,
        "certificate": {"arithmetic": "exact", "domain_complete": domain["status"] == "complete",
                        "search_complete": False, "all_solutions_returned": False,
                        "returned_tuples_verified": False},
    }
    if domain["status"] != "complete":
        report.update(status=domain["status"], solutions=[], solution_count_exact=None,
                      reason=domain.get("reason", "Domain construction did not complete"))
        return report
    requested = problem.get("method", "auto")
    native_eligible = n <= 62 and supports_boolean_mq(equations, field, domain["roots"])
    choose_native = requested == "native_boolean" or (
        requested == "auto" and native_eligible and n >= 12 and len(_components(equations, n)) == 1)
    if choose_native:
        if not native_eligible:
            raise ValueError("Native Boolean MQ requires characteristic two, domain {0,1}, and at most two distinct variables per monomial")
        search = solve_boolean_mq(equations, n, field, max_work=limits["max_work"],
                                 max_solutions=limits["max_solutions"])
    else:
        search = solve_on_domain(equations, n, domain["roots"], field,
                                 method=requested, max_work=limits["max_work"],
                                 max_table=limits["max_table"], max_solutions=limits["max_solutions"])
    report.update(search)
    counted = search["status"] in ("complete", "output_truncated")
    report["certificate"].update(search_complete=counted,
                                 all_solutions_returned=search["status"] == "complete",
                                 returned_tuples_verified=counted)
    return report


def json_report(report):
    """Encode field elements as hex so JSON consumers never round 83-bit values."""
    result = dict(report)
    result["solutions"] = [[hex(x) for x in solution] for solution in report["solutions"]]
    result["domain"] = dict(report["domain"])
    for name in ("roots", "affine_basis"):
        if name in result["domain"]:
            result["domain"][name] = [hex(x) for x in result["domain"][name]]
    if "affine_origin" in result["domain"]:
        result["domain"]["affine_origin"] = hex(result["domain"]["affine_origin"])
    return result
