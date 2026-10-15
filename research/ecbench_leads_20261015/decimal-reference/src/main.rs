//! Independent historical reader; no float_roundtrip/arbitrary_precision.
use serde_json::json;

fn main() {
    let tokens=["3934.0137668849675","0.37222665730402227","1.5030185107919957e-7",
        "1.2345678901234567e-9","1e-23","1e-100","1e-300","2.2250738585072014e-308",
        "5e-324","-5e-324","0.0","-0.0","1.7976931348623157e308","1e308"];
    let rows:Vec<_>=tokens.iter().map(|token| {
        let old=serde_json::from_str::<f64>(token);
        json!({"token":token,"historical_bits":old.as_ref().ok().map(|v|format!("{:016x}",v.to_bits())),
            "historical_error":old.err().map(|e|e.to_string()),
            "rust_correct_bits":token.parse::<f64>().ok().map(|v|format!("{:016x}",v.to_bits()))})
    }).collect();
    println!("{}",serde_json::to_string_pretty(&json!({"schema":"ecbench-decimal-reference/v1",
        "serde_json":"1.0.145","features":["std"],"float_roundtrip":false,"arbitrary_precision":false,"rows":rows})).unwrap());
}
