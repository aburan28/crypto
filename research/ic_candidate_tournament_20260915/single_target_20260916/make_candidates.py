"""Exact, isolated derivatives for the cold single-target continuation."""
import copy
import difflib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE = ROOT / 'runs/round-0004/source_candidates/folded_lift/source'
DESCENT = ROOT / 'runs/round-0005-batch16/source_candidates/combined_descent/source'
REL = Path('src/cryptanalysis/koblitz_index_calculus.rs')


def lazy_field(text):
    assert text.count('    field: FieldStructure,') == 2
    text = text.replace('    field: FieldStructure,', '    field: Option<FieldStructure>,')
    text = text.replace('field: FieldStructure::new(kc.n, &kc.curve.irreducible),',
        'field: matches!(opts.strategy, DecompositionStrategy::Groebner | DecompositionStrategy::Sat)\n'
        '                .then(|| FieldStructure::new(kc.n, &kc.curve.irreducible)),')
    start = text.index('fn decompose_once(')
    end = text.index('// ── Relation collection', start)
    part = text[start:end].replace('field: &FieldStructure,', 'field: Option<&FieldStructure>,')
    assert part.count('                field,') == 2
    part = part.replace('                field,', '                field.expect("algebraic strategy requires field structure"),')
    text = text[:start] + part + text[end:]
    assert text.count('&self.field,') == 3
    return text.replace('&self.field,', 'self.field.as_ref(),')


def fast_orbits(text):
    start = text.index('    // Index points for the orbit walk', text.index('fn finish_factor_base_domain('))
    end = text.index('    Some(FrobeniusFactorBase {', start)
    body = text[start:end]
    text = text[:start] + '''    let (orbits, orbit_of, signed_orbits, signed_orbit_of) =
        factor_base_orbits_fast(kc, &points)
            .or_else(|| factor_base_orbits_general(kc, &points))?;

''' + text[end:]
    fast = body.replace('HashMap<(BigUint, BigUint), usize>', 'HashMap<u64, usize>')
    for old, new in [('point_key(p)', 'p.pack()'), ('point_key(&cur)', 'cur.pack()'),
                     ('point_key(&point)', 'point.pack()'),
                     ('kc.frobenius(&cur)', 'fc.frobenius_k(cur, kc.k)'),
                     ('point_neg(&current)', 'fc.neg(current)'),
                     ('kc.frobenius(&current)', 'fc.frobenius_k(current, kc.k)'),
                     ('current != points[start]', 'current.pack() != points[start].pack()')]:
        fast = fast.replace(old, new)
    definitions = '''type FactorBaseOrbitTables = (
    Vec<Vec<usize>>, Vec<(usize, u32)>, Vec<Vec<usize>>, Vec<(usize, u32, bool)>,
);

fn factor_base_orbits_general(kc: &KoblitzCurve, points: &[BinaryPoint]) -> Option<FactorBaseOrbitTables> {
''' + body + '''    Some((orbits, orbit_of, signed_orbits, signed_orbit_of))
}

fn factor_base_orbits_fast(kc: &KoblitzCurve, points: &[BinaryPoint]) -> Option<FactorBaseOrbitTables> {
    let fc = FastCurve::new(&kc.curve)?;
    let points: Vec<_> = points.iter().map(|point| fc.lift(point)).collect();
''' + fast + '''    Some((orbits, orbit_of, signed_orbits, signed_orbit_of))
}

#[cfg(test)]
mod single_target_orbit_equivalence {
    use super::*;
    #[test]
    fn packed_orbit_tables_preserve_order_and_rejection() {
        for (n, a) in [(13, 0), (17, 1), (19, 0), (19, 1), (23, 0)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            for seed in [7, 43, 97] {
                let fb = build_subgroup_orbit_factor_base(&kc, seed, 6 * n as usize).unwrap();
                assert_eq!(factor_base_orbits_fast(&kc, &fb.points),
                           factor_base_orbits_general(&kc, &fb.points));
                let broken = &fb.points[..fb.points.len()-1];
                assert_eq!(factor_base_orbits_fast(&kc, broken),
                           factor_base_orbits_general(&kc, broken));
                assert!(factor_base_orbits_fast(&kc, broken).is_none());
            }
        }
    }
}

'''
    location = text.index('/// Hashable identity of a point')
    return text[:location] + definitions + text[location:]


def main():
    config = json.loads((ROOT/'runs/round-0004/winner-config.json').read_text())['config']
    registry = [{'id': 'incumbent', 'config': config,
                 'hypothesis': 'Frozen single-target incumbent folded_lift_batch4.'}]
    original = (BASE/REL).read_text()
    variants = [('descent', DESCENT, lambda t: t), ('lazy_field', BASE, lazy_field),
                ('fast_orbits', BASE, fast_orbits),
                ('combined', DESCENT, lambda t: lazy_field(fast_orbits(t)))]
    for name, parent, transform in variants:
        target = HERE/'sources'/name
        shutil.copytree(parent, target)
        changed = transform((parent/REL).read_text())
        (target/REL).write_text(changed)
        (HERE/(name+'.patch')).write_text(''.join(difflib.unified_diff(
            original.splitlines(True), changed.splitlines(True),
            fromfile='a/'+str(REL), tofile='b/'+str(REL))))
        registry.append({'id': name, 'source_root': str(target), 'config': copy.deepcopy(config),
                         'hypothesis': name+' as pre-registered in PLAN.md',
                         'falsification': 'Any mismatched certificate or failure of the two-metric gate.'})
    batch1 = copy.deepcopy(registry[-1])
    batch1['id'] = 'combined_batch1'
    batch1['config']['batch_trials'] = 1
    batch1['hypothesis'] = 'Reduce surplus relation work after setup and descent optimizations.'
    registry.append(batch1)
    (HERE/'candidates.json').write_text(json.dumps(registry, indent=2)+'\n')


if __name__ == '__main__':
    main()
