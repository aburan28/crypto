// Native metadata calculation for the existing in-memory walker.
use crypto_lib::cryptanalysis::isogeny_walk::walk::Node;
fn main() {
    let count: u128 = std::env::args()
        .nth(1)
        .expect("CLASS_COUNT")
        .parse()
        .unwrap();
    let bytes = std::mem::size_of::<Node>() as u128;
    println!("{{\"class_count_input\":\"{count}\",\"node_header_bytes\":{bytes},\"all_node_headers_bytes\":\"{}\",\"container_memory_limit_bytes\":\"17179869184\",\"excludes\":\"all node-owned heap allocations, edge and kernel storage, maps, records and factor-base work\"}}",count.checked_mul(bytes).unwrap());
}
