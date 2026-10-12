"""Frozen one-process F4/F5 and SAT qualification gate.

This gate is deliberately separate from the three-repetition first registration.
Its input must have passed the retained tournament verifier and natural-query
auditor before this read-only decision is used as evidence.
"""
import argparse
from pathlib import Path

from generic_backend_gate import evaluate
from oracle import require
from tournament import read, write

PANEL_SHA256 = 'd283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8'
LOST_EXPOSURES_SHA256 = 'a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41'


def evaluate_v2(summary, qualification, natural):
    require(summary.get('registration_panel_sha256') == PANEL_SHA256
            and summary.get('prior_censored_exposures_sha256') == LOST_EXPOSURES_SHA256,
            'not the fresh source-bound registration')
    result = evaluate(summary, qualification, natural, repetitions=1)
    result['registration_panel_sha256'] = PANEL_SHA256
    result['prior_censored_exposures_sha256'] = LOST_EXPOSURES_SHA256
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    result = evaluate_v2(read(bundle/'summary.json'),
                         read(bundle/'tournament/qualification.json'),
                         read(bundle/'natural-yield.json'))
    write(args.out, result, exclusive=True)
    print(result['status'])


if __name__ == '__main__':
    main()
