//! Compare deterministic answer/counter records, excluding measured clocks and RSS.
//! Build: rustc --edition=2021 -O replay.rs -o replay.
use std::{env, fs};

fn main() {
    let mut args = env::args().skip(1);
    let first = fs::read_to_string(args.next().expect("usage: replay run.tsv repeat.tsv")).unwrap();
    let second = fs::read_to_string(args.next().expect("repeat.tsv missing")).unwrap();
    let mut left = first.lines();
    let mut right = second.lines();
    let head = left.next().unwrap();
    assert_eq!(head, right.next().unwrap(), "schema mismatch");
    let stop = head
        .split('\t')
        .position(|name| name == "wall_ns")
        .expect("time boundary missing");
    let left: Vec<_> = left.collect();
    let right: Vec<_> = right.collect();
    assert_eq!(left.len(), 768);
    assert_eq!(right.len(), 768);
    let mut count = 0usize;
    for (a, b) in left.iter().zip(&right) {
        let a: Vec<_> = a.split('\t').collect();
        let b: Vec<_> = b.split('\t').collect();
        assert_eq!(
            &a[..stop],
            &b[..stop],
            "deterministic row {} differs",
            count + 1
        );
        count += 1;
    }
    assert_eq!(count, 768);
    println!("replay=identical deterministic_rows={count} compared_fields_per_row={stop} clocks_and_rss_excluded=true");
}
