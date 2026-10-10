"""Create isolated, reviewable candidates from the frozen previous winner."""
import difflib
import json
from pathlib import Path
import re
import shutil

WORK = Path(__file__).resolve().parent
ROOT = WORK.parent
BASE = ROOT / 'runs/round-0002/source'
REL = Path('src/cryptanalysis/koblitz_index_calculus.rs')


def change(text, kind):
    if kind == 'verify':
        needle = '        let g = kc.generator();\n        self.columns.iter().all(|(point, log)| {'
        replacement = '''        let g = kc.generator();
        if let Some(fc) = FastCurve::new(&kc.curve) {
            let g = fc.lift(g);
            return self.columns.iter().all(|(point, log)| {
                *point != BinaryPoint::Infinity && log < &kc.subgroup_order
                    && fc.mul(g, log) == fc.lift(point)
            });
        }
        self.columns.iter().all(|(point, log)| {'''
        assert text.count(needle) == 1
        text = text.replace(needle, replacement, 1)
        a = text.index('pub fn verify_collected_relation(')
        i = text.index('    let mut sum = BinaryPoint::Infinity;', a)
        text = text[:i] + '''    if let Some(fc) = FastCurve::new(&kc.curve) {
        let sum = rel.points.iter().fold(FastPoint::INFINITY,
            |sum, &i| fc.add(sum, fc.lift(&fb.points[i])));
        return !sum.infinity && sum == fc.mul_u64(fc.lift(kc.generator()), rel.a);
    }
''' + text[i:]
    elif kind == 'cofactor':
        a = text.index('    pub fn distinct_cofactor_classes(')
        i = text.index('        let mut classes = HashMap::new();', a)
        text = text[:i] + '        let fast = FastCurve::new(&kc.curve);\n' + text[i:]
        needle = '            let mut current = kc.mul(&self.points[representative], &kc.subgroup_order);'
        replacement = '''            let mut current = if let Some(fc) = &fast {
                fc.lower(fc.mul(fc.lift(&self.points[representative]), &kc.subgroup_order))
            } else {
                kc.mul(&self.points[representative], &kc.subgroup_order)
            };'''
        assert text.count(needle) == 1
        text = text.replace(needle, replacement)
    elif kind == 'serial':
        a = text.index('    pub fn build_within(')
        following = re.search(r'\n    (?:pub )?fn ', text[a + 1:])
        b = a + 1 + following.start()
        block = text[a:b]
        assert '.into_par_iter()' in block and '.par_sort_unstable()' in block
        block = block.replace('.into_par_iter()', '.into_iter()').replace('.par_sort_unstable()', '.sort_unstable()').replace('.flat_map_iter(', '.flat_map(')
        text = text[:a] + block + text[b:]
    else:
        raise ValueError(kind)
    return text


def main():
    original = (BASE / REL).read_text()
    config = json.loads((ROOT / 'runs/round-0002/winner-config.json').read_text())['config']
    registry = [{'id': 'incumbent', 'parent': 'round-0002/batch16',
                 'hypothesis': 'Retain frozen previous winner.', 'config': config},
                {'id': 'batch8', 'parent': 'incumbent',
                 'hypothesis': 'Reduce surplus verification by halving the batch.',
                 'config': dict(config, batch_trials=8)}]
    variants = [('serial_pairs', ['serial'], 16), ('fast_verify', ['verify'], 16),
                ('fast_cofactor', ['cofactor'], 16),
                ('combined', ['serial', 'verify', 'cofactor'], 16),
                ('combined_batch8', ['serial', 'verify', 'cofactor'], 8)]
    hypotheses = {
        'serial_pairs': 'Avoid Rayon scheduling and spin work for a one-CPU pair table.',
        'fast_verify': 'Verify every relation and log with equivalent single-word arithmetic.',
        'fast_cofactor': 'Compute cofactor classes with equivalent single-word arithmetic.',
        'combined': 'Measure the three source changes together, retaining every check.',
        'combined_batch8': 'Measure the combined source with fewer surplus relations.',
    }
    for name, kinds, batch in variants:
        source = WORK / 'sources' / name
        if not source.exists():
            shutil.copytree(BASE, source)
        assert (source / REL).read_text() == original, 'Do not replace an existing candidate'
        changed = original
        for kind in kinds:
            changed = change(changed, kind)
        (source / REL).write_text(changed)
        (WORK / (name + '.patch')).write_text(''.join(difflib.unified_diff(
            original.splitlines(True), changed.splitlines(True),
            fromfile='a/' + str(REL), tofile='b/' + str(REL))))
        registry.append({'id': name, 'source_root': str(source), 'parent': 'incumbent',
                         'config': dict(config, batch_trials=batch), 'hypothesis': hypotheses[name],
                         'falsification': 'Wrong certificate or changed base rejects. Promotion requires >=20% lower instructions and native wall, upper95<1 and every cell<=1.10, confirmation and replay.'})
    (WORK / 'round-0003-candidates.json').write_text(json.dumps(registry, indent=2) + '\n')
    print('Created six concrete challengers.')


if __name__ == '__main__':
    main()
