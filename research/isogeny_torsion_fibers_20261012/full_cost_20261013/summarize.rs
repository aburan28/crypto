//! Native paired target summary; rustc --edition=2021 -O summarize.rs -o summarize.
use std::collections::BTreeMap;
use std::{env, fs};

#[derive(Clone, Copy, Default)]
struct Totals {
    wall: u64,
    cpu: u64,
    ops: u64,
    map: u64,
    setup: u64,
    collection: u64,
    recovery: u64,
    descent: u64,
    verify: u64,
    rho: u64,
    other: u64,
    lifts: u64,
    lookups: u64,
    rows: u64,
    steps: u64,
    rss: u64,
    n: u64,
}
impl Totals {
    fn add(&mut self, b: Self) {
        self.wall += b.wall;
        self.cpu += b.cpu;
        self.ops += b.ops;
        self.map += b.map;
        self.setup += b.setup;
        self.collection += b.collection;
        self.recovery += b.recovery;
        self.descent += b.descent;
        self.verify += b.verify;
        self.rho += b.rho;
        self.other += b.other;
        self.lifts += b.lifts;
        self.lookups += b.lookups;
        self.rows += b.rows;
        self.steps += b.steps;
        self.rss = self.rss.max(b.rss);
        self.n += b.n;
    }
}
struct SplitMix64(u64);
impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
}
fn parse(row: &[&str], columns: &BTreeMap<&str, usize>, key: &str) -> u64 {
    row[columns[key]].parse().unwrap()
}
fn interval(targets: &[[Totals; 4]; 8], a: usize, b: usize, seed: u64, cpu: bool) -> (f64, f64) {
    let mut rng = SplitMix64(seed);
    let mut ratios = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let mut aa = 0u64;
        let mut bb = 0u64;
        for _ in 0..8 {
            let idx = (rng.next() % 8) as usize;
            aa += if cpu {
                targets[idx][a].cpu
            } else {
                targets[idx][a].wall
            };
            bb += if cpu {
                targets[idx][b].cpu
            } else {
                targets[idx][b].wall
            };
        }
        ratios.push(aa as f64 / bb as f64);
    }
    ratios.sort_by(|a, b| a.total_cmp(b));
    (ratios[250], ratios[9749])
}
fn main() {
    let raw = fs::read_to_string(env::args().nth(1).expect("usage: summarize raw.tsv")).unwrap();
    let mut lines = raw.lines();
    let names: Vec<_> = lines.next().unwrap().split('\t').collect();
    let columns: BTreeMap<&str, usize> = names.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let mut data = BTreeMap::<(u64, u64), [[Totals; 4]; 8]>::new();
    let mut rows = 0usize;
    for line in lines {
        let row: Vec<_> = line.split('\t').collect();
        assert_eq!(row.len(), names.len());
        if parse(&row, &columns, "warmup") == 1 {
            continue;
        }
        assert_eq!(row[columns["status"]], "verified");
        let p = parse(&row, &columns, "p");
        let r = parse(&row, &columns, "r");
        let target = parse(&row, &columns, "target_index") as usize;
        assert!(target < 8);
        let arm = match row[columns["arm"]] {
            "source_relations" => 0,
            "transport_relations" => 1,
            "rho_negation" => 2,
            "rho_aa" => 3,
            _ => panic!("unknown arm"),
        };
        let x = Totals {
            wall: parse(&row, &columns, "wall_ns"),
            cpu: parse(&row, &columns, "cpu_ns"),
            ops: parse(&row, &columns, "total_gae"),
            map: parse(&row, &columns, "map_ns"),
            setup: parse(&row, &columns, "setup_ns"),
            collection: parse(&row, &columns, "collection_ns"),
            recovery: parse(&row, &columns, "recovery_ns"),
            descent: parse(&row, &columns, "descent_ns"),
            verify: parse(&row, &columns, "verification_ns"),
            rho: parse(&row, &columns, "rho_ns"),
            other: parse(&row, &columns, "other_ns"),
            lifts: parse(&row, &columns, "fb_lifts"),
            lookups: parse(&row, &columns, "relation_lookups")
                + parse(&row, &columns, "descent_lookups"),
            rows: parse(&row, &columns, "row_mul_subs"),
            steps: parse(&row, &columns, "rho_steps"),
            rss: parse(&row, &columns, "process_peak_rss_kib"),
            n: 1,
        };
        data.entry((p, r)).or_default()[target][arm].add(x);
        rows += 1;
    }
    assert_eq!(rows, 640);
    assert_eq!(data.len(), 4);
    println!("p\tr\tarm\tverified_rows\tmean_cpu_ms\tcpu_vs_rho\tcpu_ci_lo\tcpu_ci_hi\tcpu_vs_source\tcpu_source_ci_lo\tcpu_source_ci_hi\tmean_wall_ms\twall_vs_rho\twall_ci_lo\twall_ci_hi\tmean_group_ops\tS_group\tgroup_vs_floor\tgroup_vs_rho\tmean_map_ms\tmean_setup_ms\tmean_collection_ms\tmean_recovery_ms\tmean_descent_ms\tmean_verification_ms\tmean_rho_ms\tmean_other_ms\ttotal_lifts\ttotal_lookups\ttotal_row_ops\ttotal_rho_steps\tprocess_peak_rss_kib");
    let labels = [
        "source_relations",
        "transport_relations",
        "rho_negation",
        "rho_aa",
    ];
    for ((p, r), targets) in data {
        let mut all = [Totals::default(); 4];
        for target in &targets {
            for (arm, &item) in target.iter().enumerate() {
                all[arm].add(item);
            }
        }
        assert!(all.iter().all(|a| a.n == 40));
        assert_eq!(all[2].ops, all[3].ops);
        for (arm, t) in all.iter().enumerate() {
            let (cpu_lo, cpu_hi) = interval(
                &targets,
                arm,
                2,
                0x4953_4f46_554c_4c43 ^ p ^ (arm as u64),
                true,
            );
            let (lo, hi) = interval(
                &targets,
                arm,
                2,
                0x5741_4c4c_5449_4d45 ^ p ^ (arm as u64),
                false,
            );
            let (source_lo, source_hi) = interval(
                &targets,
                arm,
                0,
                0x534f_5552_4345_4349 ^ p ^ (arm as u64),
                true,
            );
            let ms = |n: u64| format!("{:.6}", n as f64 / 40.0 / 1e6);
            let fields = vec![
                p.to_string(),
                r.to_string(),
                labels[arm].into(),
                t.n.to_string(),
                ms(t.cpu),
                format!("{:.6}", t.cpu as f64 / all[2].cpu as f64),
                format!("{cpu_lo:.6}"),
                format!("{cpu_hi:.6}"),
                format!("{:.6}", t.cpu as f64 / all[0].cpu as f64),
                format!("{source_lo:.6}"),
                format!("{source_hi:.6}"),
                ms(t.wall),
                format!("{:.6}", t.wall as f64 / all[2].wall as f64),
                format!("{lo:.6}"),
                format!("{hi:.6}"),
                format!("{:.3}", t.ops as f64 / 40.0),
                format!("{:.6}", t.ops as f64 / 40.0 / (r as f64).sqrt()),
                format!(
                    "{:.6}",
                    t.ops as f64 / 40.0 / (r as f64).sqrt() / (std::f64::consts::PI / 4.0).sqrt()
                ),
                format!("{:.6}", t.ops as f64 / all[2].ops as f64),
                ms(t.map),
                ms(t.setup),
                ms(t.collection),
                ms(t.recovery),
                ms(t.descent),
                ms(t.verify),
                ms(t.rho),
                ms(t.other),
                t.lifts.to_string(),
                t.lookups.to_string(),
                t.rows.to_string(),
                t.steps.to_string(),
                t.rss.to_string(),
            ];
            println!("{}", fields.join("\t"));
        }
    }
}
