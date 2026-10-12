#[allow(dead_code)]
#[path = "worker.rs"]
mod legacy;
fn main() {
    legacy::diagnostic_main();
}
