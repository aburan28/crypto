"""Reverse exact preimages of compact partial rational-map compositions.

Maps are supplied inner-to-outer. Every literal denominator must be nonzero
at the point where its map is evaluated. Common factors are deliberately not
cancelled: this API describes a partial arithmetic circuit rather than a
single element of K(X). Local degrees stay bounded even when an expanded
composition would have an enormous degree.
"""

from .polynomial import field_roots, evaluate_poly, trim


def composition_roots(
    maps,
    field,
    *,
    max_domain=4096,
    max_degree=256,
    max_work=1_000_000,
    max_trials=2048,
    seed=0,
):
    """Certify all inputs whose defined composition has terminal value zero.

    ``max_work`` caps inverse problems, not field operations. ``max_trials``
    is shared across all randomized root splits. ``max_domain`` caps each
    intermediate preimage set; refusal can occur even if an earlier map would
    subsequently shrink that set. A zero inverse polynomial is enumerated
    only when the entire field fits ``max_domain``.
    """
    for name, value in (
        ("max_domain", max_domain), ("max_degree", max_degree),
        ("max_work", max_work), ("max_trials", max_trials),
    ):
        if type(value) is not int or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if type(seed) is not int:
        raise ValueError("seed must be an integer")
    if not isinstance(maps, (list, tuple)) or not maps:
        raise ValueError("A composition must contain a nonempty list of maps")

    prepared = []
    for index, mapping in enumerate(maps):
        if not isinstance(mapping, dict) or "coefficients" not in mapping:
            raise ValueError("Every map must supply ascending numerator coefficients")
        numerator = trim(mapping["coefficients"], field)
        denominator = trim(mapping.get("denominator", (1,)), field)
        if not denominator:
            raise ValueError(f"Map {index} has a zero denominator polynomial")
        prepared.append((numerator, denominator))

    result = {
        "roots": [],
        "status": "complete",
        "method": "reverse_preimages",
        "inverse_problems": 0,
        "intermediate_sizes": [],
        "split_trials": 0,
        "excluded_pole_preimages": 0,
        "rational_semantics": "partial composition circuit: literal denominators excluded at every stage; no cancellation",
        "composition_order": "maps inner-to-outer; intermediate_sizes outer-to-inner",
        "reason": "Every inverse set is certified exactly and every final root is verified through the original circuit",
    }

    def stopped(status, reason):
        result.update(status=status, roots=[], reason=reason)
        return result

    for index, (numerator, denominator) in enumerate(prepared):
        if max(len(numerator), len(denominator)) - 1 > max_degree:
            return stopped("budget_exceeded", f"Map {index} exceeds the local degree cap {max_degree}")

    targets = {0}
    for stage in range(len(prepared) - 1, -1, -1):
        numerator, denominator = prepared[stage]
        preimages = set()
        for target in sorted(targets):
            if result["inverse_problems"] >= max_work:
                return stopped("budget_exceeded", f"The inverse-problem cap {max_work} was reached before processing map {stage}")
            result["inverse_problems"] += 1
            # A/B=y is equivalent to A-y*B=0 in every characteristic.
            inverse = trim(tuple(
                field.sub(numerator[i] if i < len(numerator) else 0,
                          field.mul(target, denominator[i] if i < len(denominator) else 0))
                for i in range(max(len(numerator), len(denominator)))
            ), field)
            if not inverse:
                # Avoid the root isolator's own unrestricted enumeration,
                # which has a separate, potentially much larger threshold.
                if field.order > max_domain:
                    return stopped(
                        "budget_exceeded",
                        f"An unrestricted inverse at map {stage} requires field enumeration ({field.order} elements), exceeding max_domain={max_domain}",
                    )
                candidates = range(field.order)
            else:
                isolated = field_roots(
                    inverse, field,
                    seed=seed + result["inverse_problems"] - 1,
                    max_trials=max_trials - result["split_trials"],
                )
                result["split_trials"] += isolated["split_trials"]
                if isolated["status"] != "complete":
                    return stopped("unresolved", f"Inverse root splitting at map {stage} is unresolved: {isolated['reason']}")
                candidates = isolated["roots"]
            for candidate in candidates:
                if evaluate_poly(denominator, candidate, field) == 0:
                    result["excluded_pole_preimages"] += 1
                    continue
                preimages.add(candidate)
                if len(preimages) > max_domain:
                    return stopped(
                        "budget_exceeded",
                        f"The intermediate preimage set at map {stage} exceeds max_domain={max_domain}",
                    )
        targets = preimages
        result["intermediate_sizes"].append(len(targets))

    roots = sorted(targets)
    for root in roots:
        value = root
        for numerator, denominator in prepared:
            divisor = evaluate_poly(denominator, value, field)
            if divisor == 0:
                raise ArithmeticError("A composition root encountered a literal circuit pole")
            value = field.mul(field.inv(divisor), evaluate_poly(numerator, value, field))
        if value != 0:
            raise ArithmeticError("Exact forward composition verification failed")
    result.update(roots=roots, domain_size=len(roots), exact_domain_count=True)
    return result
