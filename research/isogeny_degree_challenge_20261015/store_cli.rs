//! Minimal maintainer build of the repository's unchanged artifact-store CLI.
#![allow(dead_code)]
#[path = "../../src/hash/sha256.rs"]
pub mod digest_source;
mod hash {
    pub use crate::digest_source as sha256;
}
#[path = "../../src/store/mod.rs"]
mod store;
use clap::{Parser, Subcommand};
#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}
#[derive(Subcommand)]
enum Commands {
    Store(store::cli::StoreArgs),
}
fn main() {
    let Cli {
        command: Commands::Store(args),
    } = Cli::parse();
    if let Err(e) = store::cli::run(&args) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
