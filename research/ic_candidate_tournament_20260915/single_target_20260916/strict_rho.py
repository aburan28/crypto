"""Report predeclared strict rho beating, separately from frozen IC promotion."""
import argparse
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--round', required=True, type=Path)
    args = parser.parse_args()
    root = args.round.resolve()
    sys.path.insert(0, str(root/'evaluator'))
    import tournament as frozen
    assert frozen.read(root/'audit.json')['status'] == 'VERIFIED'
    contract, fixtures, arms = frozen.frozen_inputs(root)
    decision = frozen.read(root/'decision.json')
    candidates = list(dict.fromkeys([decision['winner'], decision['provisional_challenger']]))
    results = {}
    for candidate in candidates:
        if candidate is None:
            continue
        comparisons = {}
        for stage in ('confirmation', 'replay'):
            active = frozen.stage_arms(root, stage, arms)
            rows = frozen.load_stage(root, stage, fixtures[stage], active, contract['repetitions'])
            comparisons[stage] = frozen.comparison(rows, candidate, baseline='rho',
                                                    draws=contract['bootstrap_draws'])
        passed = all(p.get('eligible') and p['ci95'][1] < 1
                     and p['native_wall_ci95'][1] < 1
                     and max(p['per_cell'].values()) < 1
                     and max(p['native_wall_per_cell'].values()) < 1
                     for p in comparisons.values())
        results[candidate] = dict(strict_beats_rho=passed, comparisons=comparisons)
    result = dict(round=root.name, target_count=contract.get('target_count', 1),
                  definition='Both metric upper paired 95% limits and every cell ratio <1, confirmation AND replay.',
                  scope='Same-host bounded engineering measurement on five small Koblitz cells; no complexity claim.',
                  promotion_winner=decision['winner'], candidates=results)
    output = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    path = root/'strict-rho.json'
    if path.exists():
        assert path.read_text() == output
    else:
        path.write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
