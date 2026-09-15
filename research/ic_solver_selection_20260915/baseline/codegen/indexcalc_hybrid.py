"""ONB Redis + bounded CPU/CUDA/host-RDMA XOR preprocessing entry point."""
import argparse
import json
import os
import indexcalc_e2e as e
from onb_artifacts import ArtifactCache


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--m', type=int, choices=(5,9,131), required=True)
    p.add_argument('--variant', choices=('reference','candidate'), required=True)
    p.add_argument('--backend', choices=('auto-cuda','auto-rdma-cpu','auto-rdma-cuda'), default='auto-cuda')
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--workload', choices=('stage','dlp'), default='stage')
    p.add_argument('--trials', type=int, default=20)
    p.add_argument('--points', type=int, choices=range(3,17), default=3)
    args = p.parse_args()
    if args.points != 3 and args.m != 131:
        p.error('longer-chain mode is a degree-131 stage diagnostic')
    if args.m == 131 and args.workload == 'dlp':
        p.error('degree-131 scalar recovery is not implemented')
    os.environ['ONB_F2_BACKEND'] = 'cpu' if args.variant == 'reference' else args.backend
    cache = ArtifactCache.fromEnvironment() if os.environ.get('ONB_REDIS_URL') else None
    report = e.experiment(args.m,'candidate',args.seed,args.workload,args.trials,
                          weight=4 if args.m==5 else 2,points=args.points,artifactCache=cache)
    report.update(variant=args.variant,algebra_variant='candidate',xor_backend=os.environ['ONB_F2_BACKEND'],
                  cache_mode='redis' if cache else 'disabled')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
