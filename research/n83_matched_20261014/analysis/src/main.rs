use std::{env, fs::OpenOptions, io::Write};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: n83-matched-analysis INPUT_JSON NEW_OUTPUT_JSON");
        std::process::exit(2);
    }
    let result = (|| -> Result<(), Box<dyn std::error::Error>> {
        let input = serde_json::from_slice(&std::fs::read(&args[1])?)?;
        let output = n83_matched_analysis::summarize(&input);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&args[2])?;
        serde_json::to_writer_pretty(&mut file, &output)?;
        file.write_all(b"\n")?;
        println!(
            "{}",
            serde_json::json!({"output":args[2], "schema":output["schema"], "all_required_rows_valid":output["all_required_rows_valid"]})
        );
        Ok(())
    })();
    if let Err(error) = result {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
