//! Native read-only export; existing graph/replay certificates remain authoritative.
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, io::Write, path::Path};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
fn graph(root: &Path, name: &str) -> Result<Value> {
    let path = root.join(format!("{name}.sqlite"));
    let hash = sha256_hex(&fs::read(&path)?);
    let db = Connection::open_with_flags(&path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
    let plan: String = db.query_row("SELECT value FROM metadata WHERE key='plan'", [], |r| {
        r.get(0)
    })?;
    let plan: Value = serde_json::from_str(&plan)?;
    let verification: String = db.query_row(
        "SELECT value FROM metadata WHERE key='verification'",
        [],
        |r| r.get(0),
    )?;
    let verification: Value = serde_json::from_str(&verification)?;
    let mut nodes = Vec::new();
    let mut strata = BTreeMap::<u64, u64>::new();
    let mut statement = db.prepare("SELECT id,a,b,j,conductor,expanded FROM nodes ORDER BY id")?;
    for row in statement.query_map([], |r| {
        Ok((
            r.get::<_, u64>(0)?,
            r.get::<_, u64>(1)?,
            r.get::<_, u64>(2)?,
            r.get::<_, u64>(3)?,
            r.get::<_, u64>(4)?,
            r.get::<_, u64>(5)?,
        ))
    })? {
        let (id, a, b, j, f, expanded) = row?;
        if id != nodes.len() as u64 || expanded != 1 || f > 9_007_199_254_740_991 {
            return Err(
                "export requires contiguous ids, complete expansion and safe layout integers"
                    .into(),
            );
        }
        *strata.entry(f).or_default() += 1;
        // Large field values stay strings: JavaScript Number cannot hold them.
        nodes.push(json!([id, a.to_string(), b.to_string(), j.to_string(), f]));
    }
    if plan["expected_vertices"].as_u64() != Some(nodes.len() as u64) {
        return Err("incomplete vertex population".into());
    }
    let planned: BTreeMap<u64, u64> = plan["strata"]
        .as_array()
        .ok_or("strata")?
        .iter()
        .map(|s| {
            Ok((
                s["conductor"].as_u64().ok_or("conductor")?,
                s["expected_vertices"].as_u64().ok_or("population")?,
            ))
        })
        .collect::<Result<_>>()?;
    if strata != planned {
        return Err("stratum population mismatch".into());
    }
    let mut edges = Vec::new();
    let mut statement =
        db.prepare("SELECT source,ell,target,u,kernel FROM edges ORDER BY source,ell,kernel")?;
    for row in statement.query_map([], |r| {
        Ok((
            r.get::<_, u64>(0)?,
            r.get::<_, u64>(1)?,
            r.get::<_, u64>(2)?,
            r.get::<_, u64>(3)?,
            r.get::<_, String>(4)?,
        ))
    })? {
        let (source, ell, target, u, kernel) = row?;
        if source >= nodes.len() as u64
            || target >= nodes.len() as u64
            || !plan["degrees"]
                .as_array()
                .ok_or("degrees")?
                .contains(&json!(ell))
        {
            return Err("invalid export edge".into());
        }
        edges.push(json!([source, ell, target, u.to_string(), kernel]));
    }
    Ok(
        json!({"name":name,"source_slug":plan["source_identity"]["slug"],"p":plan["input"]["p"].as_u64().ok_or("p")?.to_string(),"trace":plan["trace"].to_string(),"order":plan["order"].as_u64().ok_or("order")?.to_string(),"D_K":plan["fundamental_discriminant"],"frobenius_conductor":plan["frobenius_conductor"],"degrees":plan["degrees"],"sha256":hash,"nodes":nodes,"edges":edges,"strata":strata,"recorded_verification":verification}),
    )
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err("usage: export GRAPH_DIR CHECKPOINT_INVENTORY TEMPLATE FRESH_HTML".into());
    }
    let inventory: Value = serde_json::from_slice(&fs::read(&args[2])?)?;
    if inventory["ok"] != true {
        return Err("checkpoint inventory did not pass".into());
    }
    let graphs = ["two", "three", "mixed"]
        .iter()
        .map(|n| graph(Path::new(&args[1]), n))
        .collect::<Result<Vec<_>>>()?;
    let counts: Vec<_> = graphs.iter().map(|g|json!({"graph":g["name"],"vertices":g["nodes"].as_array().unwrap().len(),"maps":g["edges"].as_array().unwrap().len(),"strata":g["strata"].as_object().unwrap().len()})).collect();
    let data = json!({"schema":"volcano-graph-view/v1","node_columns":["id","a","b","j","conductor"],"edge_columns":["source","ell","target","u","kernel_json"],"graphs":graphs,"checkpoint_inventory":inventory});
    let data = serde_json::to_string(&data)?
        .replace('<', "\\u003c")
        .replace('>', "\\u003e")
        .replace('&', "\\u0026");
    let template = fs::read_to_string(&args[3])?;
    if template.matches("__VOLCANO_DATA__").count() != 1 {
        return Err("expected exactly one data placeholder".into());
    }
    let html = template.replace("__VOLCANO_DATA__", &data);
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[4])?;
    out.write_all(html.as_bytes())?;
    out.sync_all()?;
    println!(
        "{}",
        json!({"output":args[4],"bytes":html.len(),"sha256":sha256_hex(html.as_bytes()),"graphs":counts,"scope":"Read-only lossless export and population checks; not a new mathematical map replay."})
    );
    Ok(())
}
