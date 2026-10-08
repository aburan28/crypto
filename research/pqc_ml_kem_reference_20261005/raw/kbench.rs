fn main() {
    for (n, c) in mlkem_dev::ml_kem::kernel_cycles() {
        println!("{n}: {c}");
    }
}
