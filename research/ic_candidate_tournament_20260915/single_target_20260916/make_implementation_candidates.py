"""Construct exact-support implementation ablations from the two-orbit source."""
import difflib
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
REL = Path('src/cryptanalysis/koblitz_index_calculus.rs')
BASE = HERE/'policy-sources/policy'

HELPER = '''fn project_by_signed_orbit(
    kc: &KoblitzCurve, fb: &FrobeniusFactorBase, fc: &FastCurve,
) -> Option<Vec<FastPoint>> {
    if fb.signed_orbit_of.len() != fb.points.len() { return None; }
    let mut conjugates = Vec::with_capacity(fb.signed_orbits.len());
    for orbit in &fb.signed_orbits {
        let representative = fb.points.get(*orbit.first()?)?;
        let mut current = fc.mul(fc.lift(representative), &kc.cofactor);
        let mut powers = Vec::with_capacity(kc.n as usize);
        for _ in 0..kc.n {
            powers.push(current);
            current = fc.frobenius_k(current, kc.k);
        }
        conjugates.push(powers);
    }
    fb.signed_orbit_of.iter().map(|&(orbit, k, negated)| {
        let point = *conjugates.get(orbit)?.get(k as usize)?;
        Some(if negated { fc.neg(point) } else { point })
    }).collect()
}

#[cfg(test)]
mod single_target_projection_equivalence {
    use super::*;
    #[test]
    fn signed_orbit_projection_matches_every_point() {
        for (n, a) in [(13, 0), (17, 1), (19, 0), (19, 1), (23, 0)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fc = FastCurve::new(&kc.curve).unwrap();
            for seed in [7, 43, 97] {
                for count in [4, 6] {
                    let mut fb = build_subgroup_orbit_factor_base(&kc, seed, count*n as usize).unwrap();
                    let expected: Vec<_> = fb.points.iter()
                        .map(|p| fc.mul(fc.lift(p), &kc.cofactor).pack()).collect();
                    let actual: Vec<_> = project_by_signed_orbit(&kc, &fb, &fc).unwrap()
                        .iter().map(|p| p.pack()).collect();
                    assert_eq!(actual, expected);
                    fb.signed_orbit_of[0].0 = usize::MAX;
                    assert!(project_by_signed_orbit(&kc, &fb, &fc).is_none());
                }
            }
        }
    }
}

'''


def projection(text):
    old = '''    let projected: Vec<FastPoint> = fb
        .points
        .par_iter()
        .map(|point| fc.mul(fc.lift(point), &kc.cofactor))
        .collect();'''
    assert text.count(old) == 1
    text = text.replace(old, '''    let projected = project_by_signed_orbit(kc, fb, &fc).unwrap_or_else(|| {
        fb.points.iter().map(|point| fc.mul(fc.lift(point), &kc.cofactor)).collect()
    });''')
    return text.replace('fn projected_signed_orbit_map_fast(', HELPER+'fn projected_signed_orbit_map_fast(', 1)


def serial(text):
    old = 'let mut relations: Vec<CollectedRelation> = (unit.start..end)\n            .into_par_iter()'
    assert text.count(old) == 1
    return text.replace(old, 'let mut relations: Vec<CollectedRelation> = (unit.start..end)\n            .into_iter()')


def main():
    original = (BASE/REL).read_text()
    config = dict(batch_trials=1, linear_algebra='sparse', max_trials=4096,
                  solver='pair_table', summands=3, factor_base_orbits=2)
    registry = [dict(id='incumbent', config=config,
                     hypothesis='Promoted two-orbit policy, unchanged complete cold cost.')]
    for name, changed in [('orbit_projection', projection(original)),
                          ('serial_collection', serial(original)),
                          ('combined', serial(projection(original)))]:
        target = HERE/'implementation-sources'/name
        shutil.copytree(BASE, target)
        (target/REL).write_text(changed)
        (HERE/('implementation-'+name+'.patch')).write_text(''.join(difflib.unified_diff(
            original.splitlines(True), changed.splitlines(True),
            fromfile='a/'+str(REL), tofile='b/'+str(REL))))
        registry.append(dict(id=name, config=config, source_root=str(target),
                             hypothesis='Fixed-support engineering: '+name+'. See IMPLEMENTATION_PLAN.md.'))
    registry.append(dict(id='combined_dense', config=dict(config, linear_algebra='dense'),
                         source_root=str(HERE/'implementation-sources/combined'),
                         hypothesis='Combined implementation with dense two-column elimination.'))
    (HERE/'implementation-candidates.json').write_text(json.dumps(registry, indent=2)+'\n')


if __name__ == '__main__':
    main()
