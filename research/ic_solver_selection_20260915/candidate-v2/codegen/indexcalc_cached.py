"""Matched ONB cache integration entry point; both arms use candidate algebra."""
import argparse
import json
import indexcalc_e2e as e
from onb_artifacts import ArtifactCache


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--m', type=int, choices=(5, 9, 131), required=True)
    p.add_argument('--variant', choices=('reference', 'candidate'), required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--workload', choices=('stage', 'dlp'), default='stage')
    p.add_argument('--trials', type=int, default=20)
    args = p.parse_args()
    if args.m == 131 and args.workload == 'dlp':
        p.error('degree-131 recovery is not implemented')
    cache = ArtifactCache.fromEnvironment() if args.variant == 'candidate' else None
    report = e.experiment(args.m, 'candidate', args.seed, args.workload, args.trials,
                          weight=4 if args.m == 5 else 2, artifactCache=cache)
    report['algebra_variant'] = 'candidate'
    report['variant'] = args.variant
    report['cache_mode'] = 'redis' if cache is not None else 'disabled'
    print(json.dumps(report))


if __name__ == '__main__':
    main()
