//! Isolated build entry point; the active primary executable is never replaced.
#[path = "../../../src/bin/ecbench.rs"]
mod driver;
fn main() -> std::process::ExitCode { driver::main() }
