//! Native diagnostic of the index metadata regression, with the same map
//! key/value types and zero-sized BuildHasherDefault as the production map.
#![allow(dead_code)]
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, DefaultHasher};
use std::mem::size_of;
type Standard = HashMap<u64, usize>;
type Fast = HashMap<u64, usize, BuildHasherDefault<DefaultHasher>>;
enum Original {
    Standard(Standard),
    Fast(Fast),
}
enum Flat {
    Standard(Standard),
    Fast(Fast),
    Dense(Vec<u16>),
}
enum FastIndex {
    Hashed(Fast),
    Dense(Vec<u16>),
}
enum Nested {
    Standard(Standard),
    Fast(FastIndex),
}
fn main() {
    println!(
        "{{\"original_index_bytes\":{},\"flat_index_bytes\":{},\"nested_index_bytes\":{}}}",
        size_of::<Original>(),
        size_of::<Flat>(),
        size_of::<Nested>()
    );
}
