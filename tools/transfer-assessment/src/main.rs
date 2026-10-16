use serde_json::{json, Value};
use std::{env, fs, path::PathBuf, process};
use transfer_assessment::validate_assessment;

fn usage() -> ! {
    eprintln!("usage: transfer-assessment [--output PATH] ASSESSMENT.json...");
    process::exit(2);
}

fn main() {
    let mut inputs = Vec::new();
    let mut output = None;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--output" {
            output = args.next().map(PathBuf::from);
            if output.is_none() {
                usage();
            }
        } else if arg.starts_with('-') {
            usage();
        } else {
            inputs.push(PathBuf::from(arg));
        }
    }
    if inputs.is_empty() {
        usage();
    }

    let mut results = Vec::new();
    let mut all_valid = true;
    for path in &inputs {
        let contents = match fs::read_to_string(path) {
            Ok(contents) => contents,
            Err(error) => {
                eprintln!("{}: {error}", path.display());
                process::exit(1);
            }
        };
        let assessment: Value = match serde_json::from_str(&contents) {
            Ok(assessment) => assessment,
            Err(error) => {
                eprintln!("{}: invalid JSON: {error}", path.display());
                process::exit(2);
            }
        };
        let validation = validate_assessment(&assessment);
        all_valid &= validation.is_valid();
        results.push(validation.to_json(&assessment));
    }

    let summary = json!({
        "format": "transfer-assessment-validation-summary-v1",
        "validator": env!("CARGO_PKG_VERSION"),
        "assessments_checked": results.len(),
        "all_valid": all_valid,
        "results": results,
    });
    let rendered = serde_json::to_string_pretty(&summary).expect("summary serialization") + "\n";
    if let Some(path) = output {
        if let Err(error) = fs::write(&path, rendered) {
            eprintln!("{}: {error}", path.display());
            process::exit(1);
        }
    } else {
        print!("{rendered}");
    }
    if !all_valid {
        process::exit(2);
    }
}
