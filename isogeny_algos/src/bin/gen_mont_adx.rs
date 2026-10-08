//! Print the MULX/ADCX/ADOX Montgomery multiplication routines embedded in `fpm::adx`.
fn main() {
    print!("{}", isogeny_algos::adx_gen::source());
}
