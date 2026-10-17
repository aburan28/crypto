import argparse
import json
import sys

from .solver import DEFAULTS, json_report, solve_problem


def main():
    parser = argparse.ArgumentParser(description="Exact bounded generalized root finding over finite fields, without Gröbner bases")
    parser.add_argument("input", help="JSON problem file")
    parser.add_argument("--output", help="Write JSON result to this path")
    for name, default in DEFAULTS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=int, default=default)
    args = parser.parse_args()
    try:
        with open(args.input, encoding="utf-8") as handle:
            problem = json.load(handle)
        report = solve_problem(problem, **{name: getattr(args, name) for name in DEFAULTS})
        result = json_report(report)
        code = 0 if result["status"] == "complete" else 3
    except (ValueError, TypeError, KeyError, ArithmeticError) as error:
        result = {"status": "invalid_input", "reason": str(error)}
        code = 2
    except RuntimeError as error:
        result = {"status": "unresolved", "reason": str(error)}
        code = 3
    output = json.dumps(result, indent=2) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(output)
    else:
        sys.stdout.write(output)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
