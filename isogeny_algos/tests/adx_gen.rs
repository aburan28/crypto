/// The generator reproduces the assembly embedded in `fpm::adx` byte for byte.
#[cfg(target_arch = "x86_64")]
#[test]
fn generator_matches_embedded_assembly() {
    use isogeny_algos::adx_gen::mont_adx;
    use isogeny_algos::fpm::adx::{MONT_ADX_4, MONT_ADX_8};
    let text = |n| {
        mont_adx(n)
            .iter()
            .map(|l| format!("{l}\n"))
            .collect::<String>()
    };
    assert_eq!(text(4), MONT_ADX_4);
    assert_eq!(text(8), MONT_ADX_8);
}
