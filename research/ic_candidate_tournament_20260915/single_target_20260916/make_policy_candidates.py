"""Build fixed-support and separately declared factor-base-policy candidates."""
import copy
import difflib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE = ROOT/'runs/round-0004/source_candidates/folded_lift/source'
COMBINED = HERE/'sources/combined'
REL = Path('src/cryptanalysis/koblitz_index_calculus.rs')
WORKER = ROOT.parents[1]/'examples/ic_tournament_worker.rs'


def cached_lifting(text):
    marker = 'fn factor_base_points_with_x(curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {'
    assert text.count(marker) == 1
    replacement = '''fn factor_base_points_with_x(curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    let fast = FastCurve::new(curve);
    factor_base_points_with_x_cached(curve, x, fast.as_ref())
}

fn factor_base_points_with_x_cached(
    curve: &BinaryCurve, x: &F2mElement, fast: Option<&FastCurve>,
) -> Vec<BinaryPoint> {'''
    text = text.replace(marker, replacement)
    start = text.index('fn factor_base_points_with_x_cached(')
    end = text.index('#[cfg(test)]', start)
    block = text[start:end]
    assert block.count('let Some(fc) = FastCurve::new(curve) else') == 1
    text = text[:start] + block.replace('let Some(fc) = FastCurve::new(curve) else',
                                      'let Some(fc) = fast else') + text[end:]
    text = text.replace('let Some(point) = factor_base_points_with_x(&kc.curve, &x).into_iter().next()',
                        'let Some(point) = factor_base_points_with_x_cached(&kc.curve, &x, Some(&curve)).into_iter().next()')
    old = '    let mut points: Vec<BinaryPoint> = Vec::new();\n    for x in &subspace {\n        for p in factor_base_points_with_x(&kc.curve, x) {'
    new = '    let fast = FastCurve::new(&kc.curve);\n    let mut points: Vec<BinaryPoint> = Vec::new();\n    for x in &subspace {\n        for p in factor_base_points_with_x_cached(&kc.curve, x, fast.as_ref()) {'
    assert text.count(old) == 1
    return text.replace(old, new)


def main():
    roots = {}
    for name, parent in [('incumbent', BASE), ('cold_context', COMBINED), ('policy', COMBINED)]:
        target = HERE/'policy-sources'/name
        shutil.copytree(parent, target)
        shutil.copy2(WORKER, target/'examples/ic_tournament_worker.rs')
        original = (parent/REL).read_text()
        changed = original if name == 'incumbent' else cached_lifting(original)
        if name == 'policy':
            assert changed.count('let batch = 8usize;') == 1
            changed = changed.replace('let batch = 8usize;', 'let batch = 1usize;')
        (target/REL).write_text(changed)
        (HERE/('policy-'+name+'.patch')).write_text(''.join(difflib.unified_diff(
            original.splitlines(True), changed.splitlines(True),
            fromfile='a/'+str(REL), tofile='b/'+str(REL))))
        roots[name] = target
    config = json.loads((ROOT/'runs/round-0004/winner-config.json').read_text())['config']
    registry = [{'id': 'incumbent', 'config': config,
                 'hypothesis': 'Last promoted single-target implementation; default base recipe unchanged.'}]
    fast_config = dict(config, batch_trials=1)
    registry.append({'id': 'cold_context', 'source_root': str(roots['cold_context']),
                     'config': fast_config, 'hypothesis': 'Reuse immutable reduction tables during all base lifts; original support.'})
    for orbits in [1, 2, 3, 4]:
        registry.append({'id': f'orbits{orbits}', 'source_root': str(roots['policy']),
                         'config': dict(fast_config, factor_base_orbits=orbits),
                         'hypothesis': f'Pay for {orbits} sampled signed orbit(s) before complete single-target recovery.'})
    registry.append({'id': 'cube_root', 'source_root': str(roots['policy']),
                     'config': dict(fast_config, factor_base_cube_root=True),
                     'hypothesis': 'Predeclared cube-root size rule balances folded setup and collection costs.'})
    (HERE/'policy-candidates.json').write_text(json.dumps(registry, indent=2)+'\n')


if __name__ == '__main__':
    main()
