//! Grouping curves by keys, and ranking curves by likeness to one.
//!
//! Both read only the precomputed [`CurveTraits::keys`] and two
//! size-free reals, so they run on the committed trait file without
//! recomputing anything.  A key whose value is `unknown` on either side
//! never counts as a match: two curves whose discriminants did not
//! factor are not thereby alike.

use std::collections::BTreeMap;

use super::CurveTraits;

/// The default weights of [`similar`]: what has to agree for two curves
/// to be alike, most structural first.
pub const DEFAULT_WEIGHTS: &[(&str, u32)] = &[
    ("char", 3),
    ("ordinary", 3),
    ("cm", 4),
    ("descent", 4),
    ("jfield", 2),
    ("cofactor", 2),
    ("twist_cofactor", 1),
    ("embedding", 1),
    ("conductor", 1),
    ("split", 1),
    ("depth", 1),
];

fn value<'a>(c: &'a CurveTraits, key: &str) -> &'a str {
    c.keys.get(key).map_or("?", String::as_str)
}

fn known(v: &str) -> bool {
    v != "unknown" && v != "?" && !v.ends_with("=?")
}

/// Curves sharing one value of every key.
#[derive(Debug)]
pub struct Group<'a> {
    pub values: Vec<String>,
    /// Ascending by field size, then slug.
    pub members: Vec<&'a CurveTraits>,
}

impl Group<'_> {
    /// The distinct field sizes the group spans, ascending.
    pub fn sizes(&self) -> Vec<u64> {
        let mut s: Vec<u64> = self.members.iter().map(|c| c.field.size_bits).collect();
        s.dedup();
        s
    }
}

/// Partition `curves` by their values of `keys`; largest groups first.
pub fn group<'a>(curves: &'a [CurveTraits], keys: &[&str]) -> Vec<Group<'a>> {
    let mut map: BTreeMap<Vec<String>, Vec<&CurveTraits>> = BTreeMap::new();
    for c in curves {
        let values = keys.iter().map(|k| value(c, k).to_string()).collect();
        map.entry(values).or_default().push(c);
    }
    let mut groups: Vec<Group> = map
        .into_iter()
        .map(|(values, mut members)| {
            members.sort_by(|a, b| (a.field.size_bits, &a.slug).cmp(&(b.field.size_bits, &b.slug)));
            Group { values, members }
        })
        .collect();
    groups.sort_by(|a, b| {
        b.members
            .len()
            .cmp(&a.members.len())
            .then_with(|| a.values.cmp(&b.values))
    });
    groups
}

/// One candidate's likeness to the target.
#[derive(Debug)]
pub struct Match<'a> {
    pub curve: &'a CurveTraits,
    /// Sum of the weights of the keys that differ (or are unknown).
    pub mismatch: u32,
    /// `|Δ trace_ratio| + |Δ conductor_fraction|`, the second taken as 1
    /// when either side is unknown.  Breaks ties among equal mismatch.
    pub distance: f64,
    pub differs: Vec<String>,
}

/// Every other curve ranked by likeness to `target`: least mismatch
/// weight first, then nearest in the two size-free reals.  With
/// `other_sizes`, curves over a field of the target's size are left out.
pub fn similar<'a>(
    curves: &'a [CurveTraits],
    target: &CurveTraits,
    weights: &[(&str, u32)],
    other_sizes: bool,
) -> Vec<Match<'a>> {
    let mut out: Vec<Match> = curves
        .iter()
        .filter(|c| c.slug != target.slug)
        .filter(|c| {
            !other_sizes
                || c.field.kind != target.field.kind
                || c.field.size_bits != target.field.size_bits
        })
        .map(|c| {
            let mut mismatch = 0;
            let mut differs = Vec::new();
            for (key, w) in weights {
                let (a, b) = (value(c, key), value(target, key));
                if a != b || !known(a) {
                    mismatch += w;
                    differs.push(if known(a) && known(b) {
                        key.to_string()
                    } else {
                        format!("{key}?")
                    });
                }
            }
            let cf = match (
                c.frobenius.conductor_fraction,
                target.frobenius.conductor_fraction,
            ) {
                (Some(x), Some(y)) => (x - y).abs(),
                _ => 1.0,
            };
            Match {
                curve: c,
                mismatch,
                distance: (c.trace_ratio - target.trace_ratio).abs() + cf,
                differs,
            }
        })
        .collect();
    out.sort_by(|a, b| {
        a.mismatch
            .cmp(&b.mismatch)
            .then(a.distance.total_cmp(&b.distance))
            .then(a.curve.slug.cmp(&b.curve.slug))
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::curve_traits::TraitFile;

    fn committed() -> TraitFile {
        serde_json::from_str(include_str!("../../../docs/curves/traits.json")).unwrap()
    }

    #[test]
    fn the_koblitz_family_is_one_group_across_sizes() {
        let file = committed();
        let groups = group(&file.curves, &["cm", "descent"]);
        let koblitz = groups
            .iter()
            .find(|g| g.values == ["-7", "k1|t|=1"])
            .expect("a Koblitz group");
        assert!(koblitz.members.iter().all(|c| c.family == "koblitz"));
        assert!(koblitz.sizes().len() >= 30);
        let sizes = koblitz.sizes();
        assert!(sizes.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn neighbours_of_ecc2k_130_share_its_cm_field_and_descent() {
        let file = committed();
        let target = file.find("ECC2K-130").expect("registered");
        let ranked = similar(&file.curves, target, DEFAULT_WEIGHTS, true);
        assert!(ranked
            .iter()
            .all(|m| m.curve.field.size_bits != target.field.size_bits));
        for m in ranked.iter().take(5) {
            assert_eq!(m.curve.keys["cm"], "-7");
            assert_eq!(m.curve.keys["descent"], "k1|t|=1");
        }
        assert!(ranked.windows(2).all(|w| w[0].mismatch <= w[1].mismatch));
    }

    #[test]
    fn unknown_values_never_match() {
        let file = committed();
        let mut a = file.curves[0].clone();
        let mut b = file.curves[1].clone();
        a.keys.insert("cm".into(), "unknown".into());
        b.keys.insert("cm".into(), "unknown".into());
        let pair = [b];
        let m = &similar(&pair, &a, &[("cm", 1)], false)[0];
        assert_eq!(m.mismatch, 1);
        assert_eq!(m.differs, ["cm?"]);
    }
}
