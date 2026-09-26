#!/usr/bin/env python3
"""Read-only necessary target-support/rank gate for a proposed PDP run.

The bound assumes every target probe is marginally uniform on the declared
group/coset, irrespective of dependence between probes. It never predicts an
observed hit rate or certifies UNSAT. See PROTOCOL.md for cofactor and row
multiplicity rules.
"""

import argparse
from fractions import Fraction
from math import comb, prod
import json


def positive(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    bases = parser.add_mutually_exclusive_group(required=True)
    bases.add_argument("--unordered-size", type=positive,
                       help="physical points in one common, unordered base")
    bases.add_argument("--ordered-slot-sizes", type=positive, nargs="+",
                       help="physical point count in each distinct ordered slot")
    parser.add_argument("--arity", type=positive,
                        help="summand count (required with --unordered-size)")
    parser.add_argument("--target-space-size", type=positive, required=True,
                        help="cardinality of the exact uniform target group/coset")
    parser.add_argument("--targets", type=positive, required=True,
                        help="all probed targets, including torsion translates")
    parser.add_argument("--required-rank", type=positive, required=True,
                        help="independent rows needed for the frozen solve")
    parser.add_argument("--rows-per-target-max", type=positive, default=1,
                        help="enforced maximum emitted relation rows per probe")
    args = parser.parse_args()

    if args.unordered_size is not None:
        if args.arity is None:
            parser.error("--arity is required with --unordered-size")
        count = comb(args.unordered_size + args.arity - 1, args.arity)
        base = {"unordered_size": args.unordered_size, "arity": args.arity}
    else:
        if args.arity is not None and args.arity != len(args.ordered_slot_sizes):
            parser.error("--arity must equal the number of ordered slots")
        count = prod(args.ordered_slot_sizes)
        base = {"ordered_slot_sizes": args.ordered_slot_sizes,
                "arity": len(args.ordered_slot_sizes)}

    hit = min(Fraction(1), Fraction(count, args.target_space_size))
    expected = args.targets * hit
    row_cap = args.targets * args.rows_per_target_max
    rank = (Fraction(0) if row_cap < args.required_rank else
            min(Fraction(1), expected * args.rows_per_target_max /
                args.required_rank))
    output = {
        "classification": "NECESSARY_BOUND_ONLY",
        "base": base,
        "support_cardinality_upper": str(count),
        "target_space_size": str(args.target_space_size),
        "targets": args.targets,
        "rows_per_target_max": args.rows_per_target_max,
        "deterministic_row_cap": row_cap,
        "required_rank": args.required_rank,
        "hit_probability_upper_fraction": f"{hit.numerator}/{hit.denominator}",
        "expected_supported_probes_upper_fraction":
            f"{expected.numerator}/{expected.denominator}",
        "full_rank_probability_upper_fraction":
            f"{rank.numerator}/{rank.denominator}",
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
