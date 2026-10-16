use super::*;
#[test]
fn workload_counts_and_crossover_matching() {
    let root = Path::new(ROOT);
    for (name, count) in [
        ("crossover", 240),
        ("yield", 400),
        ("rho", 400),
        ("orbit", 800),
        ("dose", 320),
    ] {
        let v = read(root.join(format!("prepared/{name}.json")));
        assert_eq!(v.as_array().unwrap().len(), count);
    }
    let v = read(root.join("prepared/crossover.json"));
    let mut orders = std::collections::BTreeMap::<Vec<String>, usize>::new();
    for block in v.as_array().unwrap().chunks_exact(3) {
        assert!(block
            .iter()
            .all(|x| x["pair"] == block[0]["pair"] && x["rep"] == block[0]["rep"]));
        let cm = block.iter().find(|x| x["arm"] == "crater-cm").unwrap();
        let rec = block
            .iter()
            .find(|x| x["arm"] == "crater-recurrence")
            .unwrap();
        for key in ["curve", "Qx", "Qy", "query_seed", "B", "query_limit"] {
            assert_eq!(cm[key], rec[key]);
        }
        *orders
            .entry(
                block
                    .iter()
                    .map(|x| x["arm"].as_str().unwrap().to_string())
                    .collect(),
            )
            .or_default() += 1;
    }
    assert_eq!(orders.len(), 6);
    assert!(orders.values().all(|&n| n == 13 || n == 14));
}
#[test]
fn targets_and_partition_canonicalization() {
    let j = &read(Path::new(ROOT).join("prepared/control.json"))[0];
    let i = inst(&j["curve"]);
    let q = target(&i.curve, i.r, "test/public");
    assert_eq!(q, target(&i.curve, i.r, "test/public"));
    assert_ne!(q, target(&i.curve, i.r, "test/other"));
    let canon = FixedNeg { r: i.r };
    for k in 1..30 {
        let p = i.curve.mul(&mut GroupOps::default(), q, k);
        let a = Salted {
            c: &i.curve,
            salt: 0,
            keys: Cell::new(0),
            negs: Cell::new(0),
        };
        let b = Salted {
            c: &i.curve,
            salt: 98765,
            keys: Cell::new(0),
            negs: Cell::new(0),
        };
        assert_eq!(
            canon.canon(&a, &mut GroupOps::default(), p),
            canon.canon(&b, &mut GroupOps::default(), p)
        );
    }
}
#[test]
fn independent_replay_rejects_corrupt_relation_and_query() {
    let source = Path::new(ROOT).join("runs/control-0000");
    let j = read(source.join("input.json"));
    let result = read(source.join("result.json"));
    let raw = &result["result"]["native"];
    assert_eq!(raw["verified"], true);
    let temp = std::env::temp_dir().join(format!("cm-followup-corruption-{}", std::process::id()));
    fs::create_dir_all(temp.join("factor-bases")).unwrap();
    fs::create_dir_all(temp.join("relations")).unwrap();
    fs::copy(
        source.join("relations/job-0.jsonl"),
        temp.join("relations/job-0.jsonl"),
    )
    .unwrap();
    fs::copy(
        source.join("factor-bases/job-0.csv"),
        temp.join("factor-bases/job-0.csv"),
    )
    .unwrap();
    let events = fs::read_to_string(source.join("events.jsonl")).unwrap();
    assert!(!events.is_empty());
    fs::write(temp.join("events.jsonl"), &events).unwrap();
    assert_eq!(replay(&j, &temp, raw)["ok"], true);
    for mode in ["sign", "query"] {
        let mut rows: Vec<Value> = events
            .lines()
            .map(|s| serde_json::from_str(s).unwrap())
            .collect();
        if mode == "sign" {
            rows[0]["si"] = json!(-rows[0]["si"].as_i64().unwrap());
        } else {
            rows[0]["query"] = raw["queries"].clone();
        }
        let text = rows
            .iter()
            .map(|v| serde_json::to_string(v).unwrap() + "\n")
            .collect::<String>();
        fs::write(temp.join("events.jsonl"), text).unwrap();
        assert!(std::panic::catch_unwind(|| replay(&j, &temp, raw)).is_err());
    }
    fs::remove_dir_all(temp).unwrap();
}

#[test]
fn dose_classes_hold_field_and_group_fixed() {
    let text =
        fs::read_to_string(Path::new(ROOT).join("results/dose-construction-v3.jsonl")).unwrap();
    let models: Vec<Value> = text
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    assert_eq!(models.len(), 32);
    for class in 1..=8 {
        let rows: Vec<_> = models.iter().filter(|r| r["class"] == class).collect();
        assert_eq!(rows.len(), 4);
        let strata: BTreeSet<_> = rows.iter().map(|r| u(r, "stratum")).collect();
        assert_eq!(strata, BTreeSet::from([1, 2, 3, 4]));
        for r in &rows {
            for k in ["p", "N", "t", "fpi", "D"] {
                assert_eq!(r[k], rows[0][k]);
            }
        }
        let f = u(rows[0], "fpi");
        let conductors: BTreeSet<_> = rows.iter().map(|r| u(r, "conductor")).collect();
        assert_eq!(conductors.len(), 4);
        assert!(conductors.contains(&1) && conductors.contains(&f));
        let primes: Vec<_> = conductors
            .iter()
            .filter(|&&c| c != 1 && c != f)
            .copied()
            .collect();
        assert_eq!(primes[0] * primes[1], f);
        for ell in primes {
            assert!(crypto_lib::cryptanalysis::residual_walk::is_prime_u64(ell));
        }
    }
}
