//! Attest the actual Gf2 arithmetic dispatch under the inherited environment.

use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

fn main() {
    for n in [41, 53] {
        let curve = KoblitzCurve::new(0, n).expect("frozen Koblitz curve");
        let field = Gf2::new(&curve.curve.irreducible);
        println!("n={n} kernel={} folds={}", field.kernel_name(), field.folds());
    }
}
