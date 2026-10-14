//! Build the same native publisher without the large cryptographic library.
extern crate self as crypto_lib;
pub mod cryptanalysis {
    pub mod ecbench {
        pub mod canonical {
            use sha2::{Digest, Sha256};
            pub fn sha256_hex(bytes: &[u8]) -> String {
                format!("{:x}", Sha256::digest(bytes))
            }
        }
    }
}
#[path = "../../../../examples/ecbench_open_questions_publish.rs"]
mod publish;
fn main() -> Result<(), String> {
    publish::main()
}
