//! Select a complete SMT sibling set for the Linux F5 isolation reservation.

#[cfg(target_os = "linux")]
fn main() {
    use serde_json::json;
    use std::collections::BTreeSet;
    use std::fs;

    fn cpu_list(text: &str) -> BTreeSet<usize> {
        text.trim()
            .split(',')
            .flat_map(|part| {
                let (first, last) = part.split_once('-').unwrap_or((part, part));
                first.parse::<usize>().unwrap()..=last.parse::<usize>().unwrap()
            })
            .collect()
    }
    fn csv(cpus: &BTreeSet<usize>) -> String {
        cpus.iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",")
    }

    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 3, "usage: f5_support_split_cpus THREADS OUTPUT");
    let threads: usize = args[1].parse().unwrap();
    let mut set: libc::cpu_set_t = unsafe { std::mem::zeroed() };
    let code = unsafe { libc::sched_getaffinity(0, std::mem::size_of_val(&set), &mut set) };
    assert_eq!(code, 0, "sched_getaffinity failed");
    let allowed: BTreeSet<usize> = (0..libc::CPU_SETSIZE as usize)
        .filter(|&cpu| unsafe { libc::CPU_ISSET(cpu, &set) })
        .collect();
    let mut groups: Vec<BTreeSet<usize>> = Vec::new();
    for &cpu in &allowed {
        let siblings = cpu_list(
            &fs::read_to_string(format!(
                "/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"
            ))
            .unwrap(),
        );
        if siblings.is_subset(&allowed) && !groups.contains(&siblings) {
            groups.push(siblings);
        }
    }
    let mut choice = None;
    for start in 0..groups.len() {
        let mut reserved = BTreeSet::new();
        for group in groups.iter().skip(start) {
            reserved.extend(group);
            if reserved.len() >= threads {
                if reserved.len() < allowed.len() {
                    choice = Some(reserved);
                }
                break;
            }
        }
        if choice.is_some() {
            break;
        }
    }
    let result = if let Some(reserved) = &choice {
        let pin: BTreeSet<_> = reserved.iter().take(threads).copied().collect();
        println!("reserve={}", csv(reserved));
        println!("pin={}", csv(&pin));
        json!({"status":"ok", "allowed":allowed, "reserved":reserved, "pin":pin, "threads":threads})
    } else {
        json!({"status":"no_reservable_cpu_set", "allowed":allowed, "threads":threads})
    };
    fs::write(
        &args[2],
        format!("{}\n", serde_json::to_string_pretty(&result).unwrap()),
    )
    .unwrap();
    assert!(choice.is_some(), "no reservable CPU set");
}

#[cfg(not(target_os = "linux"))]
fn main() {
    panic!("F5 CPU reservation requires Linux");
}
