// Visualizes only the recorded summary.tsv; it performs no research calculation.
use std::fmt::Write as _;

fn main() {
    let input = std::fs::read_to_string("summary.tsv").expect("read summary.tsv");
    let rows: Vec<(u64, u64, u64)> = input
        .lines()
        .skip(2)
        .map(|line| {
            let values: Vec<u64> = line
                .split('\t')
                .map(|v| v.parse::<u64>().expect("integer cell"))
                .collect();
            assert_eq!(values.len(), 3);
            (values[0], values[1], values[2])
        })
        .collect();
    assert_eq!(rows.len(), 21);
    let mut svg = String::new();
    writeln!(
        svg,
        "<svg xmlns='http://www.w3.org/2000/svg' width='1100' height='490' viewBox='0 0 1100 490'>"
    )
    .unwrap();
    writeln!(svg, "<rect width='1100' height='490' fill='white'/>").unwrap();
    writeln!(svg, "<text x='48' y='34' font-size='21' font-family='sans-serif'>Bounded census: non-rational order-3 kernels</text>").unwrap();
    writeln!(svg, "<text x='48' y='56' font-size='13' font-family='sans-serif'>All p ≤ 199, p ≡ 1 mod 3; one row per nonsquare b in y² = x³ + b</text>").unwrap();
    for tick in [0, 25, 50, 75, 100] {
        let y = 390 - tick * 3;
        writeln!(svg, "<line x1='48' y1='{y}' x2='1080' y2='{y}' stroke='#dedede'/><text x='15' y='{}' font-size='14' font-family='sans-serif'>{tick}</text>", y + 4).unwrap();
    }
    for (index, (p, all, witnesses)) in rows.iter().enumerate() {
        let x = 53 + index * 49;
        let all_height = all * 3;
        let witness_height = witnesses * 3;
        writeln!(
            svg,
            "<rect x='{x}' y='{}' width='18' height='{all_height}' fill='#a9c6dd'/>",
            390 - all_height
        )
        .unwrap();
        writeln!(
            svg,
            "<rect x='{}' y='{}' width='18' height='{witness_height}' fill='#0a6f8c'/>",
            x + 19,
            390 - witness_height
        )
        .unwrap();
        writeln!(svg, "<text x='{}' y='409' font-size='14' font-family='sans-serif' transform='rotate(45 {} 409)'>{p}</text>", x + 6, x + 6).unwrap();
    }
    writeln!(svg, "<rect x='580' y='64' width='16' height='12' fill='#a9c6dd'/><text x='602' y='75' font-size='14' font-family='sans-serif'>nonsquare b scanned</text>").unwrap();
    writeln!(svg, "<rect x='775' y='64' width='16' height='12' fill='#0a6f8c'/><text x='797' y='75' font-size='14' font-family='sans-serif'>prime-order witness</text>").unwrap();
    writeln!(svg, "<text x='48' y='470' font-size='14' font-family='sans-serif'>Source: summary.tsv; exact enumeration, no timing or DLP measurement.</text></svg>").unwrap();
    std::fs::write("census.svg", svg).expect("write census.svg");
}
