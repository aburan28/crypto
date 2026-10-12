#[path = "../producer/mod.rs"]
mod producer;

use std::process::ExitCode;

fn main() -> ExitCode {
    match producer::run(std::env::args().skip(1)) {
        Ok(output) => {
            print!("{output}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
