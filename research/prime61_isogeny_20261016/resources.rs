// Native metadata calculation for the existing in-memory walker.
use crypto_lib::cryptanalysis::isogeny_walk::walk::{Edge, Node};
fn main() {
    let count: u128 = std::env::args()
        .nth(1)
        .expect("CLASS_COUNT")
        .parse()
        .unwrap();
    let bytes = std::mem::size_of::<Node>() as u128;
    let edge_bytes = std::mem::size_of::<Edge>() as u128;
    let node_total = count.checked_mul(bytes).unwrap();
    let minimum_edges = count.checked_sub(1).unwrap();
    let combined = node_total
        .checked_add(minimum_edges.checked_mul(edge_bytes).unwrap())
        .unwrap();
    let limit: u128 = std::fs::read_to_string("/sys/fs/cgroup/memory.max")
        .expect("Linux cgroup-v2 memory limit")
        .trim()
        .parse()
        .expect("finite numeric container memory limit");
    println!("{{\"class_count_input\":\"{count}\",\"node_header_bytes\":{bytes},\"edge_header_bytes\":{edge_bytes},\"all_node_headers_bytes\":\"{node_total}\",\"minimum_spanning_edges\":\"{minimum_edges}\",\"minimum_node_and_edge_headers_bytes\":\"{combined}\",\"container_memory_limit_bytes\":\"{limit}\",\"excludes\":\"all node-owned heap allocations, kernel coefficients, hash tables, maps, records and factor-base work\"}}");
}
