"""Frozen-source complete-workload adapter; parent measures process lifetime."""
import argparse
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True)
    parser.add_argument('--params', required=True)
    parser.add_argument('--dir', required=True)
    parser.add_argument('--model')
    args = parser.parse_args()
    sys.path.insert(0, str(Path(args.source).resolve()))
    import indexcalc_fixed as fixed
    document = fixed.readJson(args.params)
    options = {} if args.model is None else {'solver': 'learned', 'selectorModel': args.model}
    with fixed.Campaign(document, args.dir, **options) as campaign:
        report = campaign.run(attempts=512, pairBudget=100000, seconds=.25)
        encode = lambda q: None if q is None else [campaign.params.onb.toCoords(x) for x in q]
        fixture = {'parameters': document, 'generator_onb': encode(campaign.params.generator),
                   'prime': str(campaign.params.prime), 'eigen': str(campaign.params.eigen),
                   'base': [encode(q) for q in campaign.points],
                   'representatives': [encode(q) for q in campaign.reps]}
        relations = []
        for stream, number, witness, row, a, b in campaign.db.execute(
                'SELECT stream,number,witness,coefficients,a,b FROM relations ORDER BY stream,number'):
            relations.append({'stream': stream, 'number': number, 'witness': json.loads(witness),
                              'row': json.loads(row), 'a': a, 'b': b,
                              'target': encode(campaign.probe(stream, number)[2])})
        result = {'report': report, 'fixture': fixture, 'relations': relations,
                  'logs': campaign.meta('logs'),
                  'attempts': [{'stream': stream, 'number': number, 'result': json.loads(details)}
                               for stream, number, details in campaign.db.execute(
                                   'SELECT stream,number,result FROM attempts ORDER BY stream,number')]}
    print(json.dumps(result, sort_keys=True))
    return 0 if report['status'] == 'complete' else 2


if __name__ == '__main__':
    sys.exit(main())
