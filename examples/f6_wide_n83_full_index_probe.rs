//! Frozen full-base probe: research/f6_wide_n83_20261004/FULL_INDEX_PROTOCOL.md.
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::f6_wide_geometry::F6WidePairIndex;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use serde_json::json;

fn peak_rss_bytes() -> Option<u64> {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        // SAFETY: getrusage fills the entire rusage output on success.
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
            return None;
        }
        let usage = unsafe { usage.assume_init() };
        #[cfg(target_os = "macos")]
        return Some(usage.ru_maxrss as u64);
        #[cfg(target_os = "linux")]
        return Some((usage.ru_maxrss as u64) * 1024);
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    None
}

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K_0 curve");
    let target = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&target));
    assert_eq!(kc.mul(&target, &kc.subgroup_order), BinaryPoint::Infinity);

    let parent = build_standard_subspace_factor_base(&kc, 12).expect("dimension-12 base");
    let base = cofactor_project_factor_base(&kc, &parent).expect("usable subgroup base");
    assert_eq!(base.points.len(), 4_054);
    let pairs = base.points.len() * (base.points.len() + 1) / 2;
    assert_eq!(pairs, 8_219_485);

    let start = Instant::now();
    let index = F6WidePairIndex::new(&kc.curve, &base.points, pairs).expect("full pair index");
    let build_ns = start.elapsed().as_nanos();
    assert_eq!(index.pair_count(), pairs);

    let start = Instant::now();
    let answer = index.solve4(&target);
    let query_ns = start.elapsed().as_nanos();
    if let Some(indices) = answer {
        let replay = indices.iter().fold(BinaryPoint::Infinity, |sum, &i| {
            kc.add(&sum, &base.points[i])
        });
        assert_eq!(replay, target);
    }
    println!(
        "{}",
        json!({
            "curve_id": kc.label(),
            "ell": 12,
            "usable_points": base.points.len(),
            "signed_columns": base.signed_orbits.len(),
            "pairs": pairs,
            "build_ns": build_ns,
            "query_ns": query_ns,
            "peak_rss_bytes": peak_rss_bytes(),
            "outcome": if answer.is_some() { "witness" } else { "no_witness" },
            "correctness": "PASS",
            "claim_scope": "exploratory_stage_only"
        })
    );
}
