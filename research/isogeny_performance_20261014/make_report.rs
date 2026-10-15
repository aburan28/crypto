//! Native presentation builder: derive tables and graphs only from a PASS comparison.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::fmt::Write;
use std::path::Path;

fn string<'a>(j: &'a Json,k: &str)-> &'a str {j.get(k).unwrap().as_str().unwrap()}
fn main() {
    let args:Vec<_>=std::env::args().collect();
    if args.get(1).is_some_and(|a| a=="--manifest") {
        assert_eq!(args.len(),4,"--manifest DIRECTORY FRESH_OUTPUT");
        let root=Path::new(&args[2]); let output=Path::new(&args[3]);
        let mut files=vec![];
        fn visit(root: &Path, at: &Path, output: &Path, files: &mut Vec<Json>) {
            let mut entries:Vec<_>=std::fs::read_dir(at).unwrap().map(Result::unwrap).collect();
            entries.sort_by_key(|e|e.file_name());
            for entry in entries {
                let path=entry.path(); if path==output {continue;}
                let kind=entry.file_type().unwrap();
                if kind.is_dir() {visit(root,&path,output,files);} else if kind.is_file() {
                    let bytes=std::fs::read(&path).unwrap();
                    files.push(Json::obj(vec![("path",Json::str(path.strip_prefix(root).unwrap().to_string_lossy())),
                        ("bytes",Json::str(bytes.len())),("sha256",Json::str(sha256_hex(&bytes)))]));
                } else {panic!("unexpected non-file: {}",path.display());}
            }
        }
        visit(root,root,output,&mut files);
        let record=Json::obj(vec![("schema",Json::str("isogeny-perf.files/v1")),("files",Json::Arr(files))]);
        use std::io::Write;
        let mut f=std::fs::OpenOptions::new().write(true).create_new(true).open(output).unwrap();
        writeln!(f,"{}",record.dump()).unwrap(); return;
    }
    assert_eq!(args.len(),3,"make-report COMPARISON.json FRESH_DIRECTORY");
    let bytes=std::fs::read(&args[1]).unwrap();
    let data=Json::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(string(&data,"status"),"PASS");
    let Json::Arr(cases)=data.get("cases").unwrap() else {panic!("missing cases")};
    let rows:Vec<_>=cases.iter().filter(|r|string(r,"bits")=="224" && string(r,"ell")=="101").collect();
    assert_eq!(rows.len(),8);
    let root=Path::new(&args[2]);
    if root.exists() {
        assert_eq!(std::fs::read(root.join("COMPARISON.json")).unwrap(), bytes, "refresh must preserve the frozen comparison");
    } else { std::fs::create_dir(root).unwrap(); }
    let mut md=String::from("# BMSS constructor optimization — 2026-10-14\n\nMeasured generic Field calls, synthetic 224-bit fixtures, prime degree 101. Counts exclude fixture preparation and verification. Inversions are opaque; this is not an instruction or wall-time ratio.\n\n| Algorithm | Baseline mul calls | Candidate mul calls | Call quotient | Baseline inv calls | Candidate inv calls |\n|---|---:|---:|---:|---:|---:|\n");
    let mut cells=String::new();
    let mut svg=String::from("<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='550' viewBox='0 0 1000 550'><rect width='1000' height='550' fill='#fff'/><g font-family='Arial,sans-serif' fill='#17324d'><text x='28' y='35' font-size='23'>BMSS multiplication-call quotient: baseline / candidate</text><text x='28' y='64' font-size='15'>224-bit public fixture, degree 101; deterministic generic Field model, one run per arm</text><text x='28' y='86' font-size='15'>Exact output equality checked in all 96 cases; wall-time improvement not established</text>");
    for (i,row) in rows.iter().enumerate() {
        let name=string(row,"algorithm");let ops=row.get("operations").unwrap();
        let m=ops.get("mul").unwrap();let inv=ops.get("inv").unwrap();
        let ratio=m.get("baseline_over_candidate").unwrap().as_f64().unwrap();
        writeln!(md,"| {name} | {} | {} | {ratio:.3} | {} | {} |",string(m,"baseline"),string(m,"candidate"),string(inv,"baseline"),string(inv,"candidate")).unwrap();
        writeln!(cells,"[{}], [{}], [{}], [{ratio:.3}], [{}], [{}],", name.replace('_'," "),string(m,"baseline"),string(m,"candidate"),string(inv,"baseline"),string(inv,"candidate")).unwrap();
        let y=120+i*43;let width=ratio.log2()*90.0;
        let x=225.0+width.min(0.0);
        let color=if ratio<1.0 {"#b65d35"} else {"#238779"};
        writeln!(svg,"<text x='28' y='{}' font-size='17'>{}</text><line x1='225' x2='225' y1='{y}' y2='{}' stroke='#b8c9c5'/><rect x='{x:.2}' y='{y}' width='{:.2}' height='26' fill='{color}'/><text x='{:.2}' y='{}' font-size='16'>{ratio:.3}</text>",y+19,name.replace('_'," "),y+26,width.abs(),235.0+width.max(0.0),y+19).unwrap();
    }
    for i in 0..=6 {let x=225+i*90; writeln!(svg,"<text x='{x}' y='492' font-size='14'>{}</text>",1u64<<i).unwrap();}
    svg.push_str("<text x='225' y='520' font-size='15'>Logarithmic axis. 1 = unchanged multiplication calls. No priced total-cost claim.</text></g></svg>");
    let hash=sha256_hex(&bytes);
    writeln!(md,"\nSource: comparison SHA256 `{hash}`; baseline run SHA256 `{}`; candidate run SHA256 `{}`.\n\nThe linear-algebra rung moves row-independent series products into a column loop, without retaining an extra polynomial cache. The second rung batches fixed series/Newton denominators above their small-size thresholds. Separate rung receipts distinguish them. A rise in multiplications for batching is reported beside the decrease in inversions.\n\nNative macOS runs used an Apple M4 Pro (14 logical cores, 48 GiB); this shared host had unrelated CPU jobs. No native wall samples were collected or promoted. The four-architecture GitHub matrix is implemented; its receipts must be examined before any other CPU measurement is credited.\n\nRemaining requirements: qualified per-CPU A/A and interleaved A/B timings, additional constructors/finders/path families in INVENTORY.md, named P192/P224 whole-search/replay totals, and actual ReleaseMe asset publication. The earlier 2,000-map search requirement remains open.\n\nCanonical graphs reviewed: docs/index-calculus-scoreboard.html, docs/ic/progress-timeline.json, docs/performance-gains/data.json and comparisons.csv, docs/curves/ and figures/. They describe DLP/rho, unrelated polynomial specialization or curve edges; this constructor operation-count result changes none of those records. It appears in the study graph and benchmarks index.\n\n![Count comparison](COUNTS.svg)\n",string(&data,"baseline_run_sha256"),string(&data,"candidate_run_sha256")).unwrap();
    let typ=format!(r##"#set page(paper: "a4", margin: 16mm)
#set text(font: "Arial", size: 9pt)
#set heading(numbering: none)
= BMSS constructor optimization
2026-10-14 · Native exact-output evidence

All eight BMSS methods at 192, 224, 256 and 511 bits, degrees 7, 31 and 101: 96 exact baseline/candidate output matches. Kernel, codomain and rational map coefficients agree with the Kohel reference. Four public point transports are checked per constructor case.

== Arithmetic counts at 224 bits, degree 101
#table(columns: (1.8fr, 1fr, 1fr, 0.8fr, 0.8fr, 0.8fr), inset: 4pt,
  [Algorithm], [Baseline mul], [Candidate mul], [Quotient], [Baseline inv], [Candidate inv],
  {cells}
)
#image("COUNTS.svg", width: 100%)

The quotient counts generic Field multiplication calls. Inversions are opaque and counted separately. Setup, preparing sigma, independent reference construction and verification are excluded. Counts are deterministic for these fixed inputs; their model is not a CPU instruction count or a priced total cost.

== Changes and limits
Row-independent series products are computed once per column. The second rung batches fixed series and Newton denominators, with direct division retained at small sizes. The third uses forward elimination followed by back substitution. Raw receipts distinguish the rungs and retain any additional multiplications incurred by batching.

This shared Apple M4 Pro host has 14 logical cores and 48 GiB RAM. Wall time is unqualified; no wall samples are promoted. Qualified multi-CPU A/A and A/B comparisons, the other inventory families, named NIST whole-search/replay totals, and actual release publication remain open. The 2,000-additional-map search remains separate and unfinished.

Comparison SHA256:
#text(size: 7pt)[{hash}]
Baseline and candidate run digests and complete case inputs are in COMPARISON.json and the source-linked raw JSONL records. Read PROTOCOL.md, INVENTORY.md and ATTEMPTS.md for boundaries and censored attempts.
"##);
    for (file,content) in [("REPORT.md",md),("REPORT.typ",typ),("COUNTS.svg",svg),
        ("graph.typ",String::from("#set page(width: 1000pt, height: 550pt, margin: 0pt)\n#image(\"COUNTS.svg\", width: 100%)\n"))] {std::fs::write(root.join(file),content).unwrap();}
    let diagram=r##"<svg xmlns="http://www.w3.org/2000/svg" width="900" height="160" viewBox="0 0 900 160"><rect width="900" height="160" fill="white"/><g font-family="Arial,sans-serif" font-size="16" fill="#17324d"><text x="15" y="24">Same public inputs and constructor boundary for each algorithm revision</text><g fill="#e7f2ef" stroke="#238779"><rect x="15" y="45" width="185" height="75"/><rect x="240" y="45" width="185" height="75"/><rect x="465" y="45" width="185" height="75"/><rect x="690" y="45" width="195" height="75"/></g><text x="25" y="73">Frozen E, E', degree</text><text x="25" y="100">Sigma when required</text><text x="250" y="73">Baseline / candidate</text><text x="250" y="100">Separate field counts</text><text x="475" y="73">Exact map fingerprint</text><text x="475" y="100">Kohel + transports</text><text x="700" y="73">A/A, then interleaved</text><text x="700" y="100">A/B + Linux isolation</text><text x="205" y="87">→</text><text x="430" y="87">→</text><text x="655" y="87">→</text><text x="15" y="146">Timing qualification comes from the native isolation receipt; counts do not establish whole-search time.</text></g></svg>"##;
    std::fs::write(root.join("PIPELINE.svg"),diagram).unwrap();
    std::fs::write(root.join("pipeline.typ"),"#set page(width: 900pt, height: 160pt, margin: 0pt)\n#image(\"PIPELINE.svg\", width: 100%)\n").unwrap();
    std::fs::write(root.join("RUN_REPORT.typ"),format!(r##"#set page(paper: "a4", margin: 18mm)
#set text(font: "Arial", size: 10pt)
= Native BMSS run record
2026-10-14 · Source-linked baseline and candidate

Requested: baseline each native isogeny-finding method across multiple processor types, then iterate on speed. This run implements the first eight-method series; the larger inventory remains open.

#image("PIPELINE.svg", width: 100%)

The baseline was frozen before algorithm edits. Deterministic supersingular fixtures span 192, 224, 256 and 511 field bits, at prime degrees 7, 31 and 101. These fixtures have their full coefficients and hashes in the raw JSONL; they are not the named P192/P224 models. Counts and native constructor verification run separately. Inversion internals are opaque.

The committed native comparison accepts all 96 exact output fingerprints and matches each input digest. Independent Kohel kernel and rational-map agreement and four public point transports guard every case. Missing terminal receipts, changed inputs and missing algorithm rows fail closed.

This machine is an Apple M4 Pro, 14 logical cores, 48 GiB, macOS arm64, compiler rustc 1.93.1 (Homebrew), as recorded by Cargo's compiler fingerprint. Cargo was invoked through rustup 1.98 but resolved the Homebrew rustc on PATH. It was shared with unrelated native jobs. No constructor wall-time sample is credited. All preparation, verification, counting, compilation and interrupted attempts are outside a claimed constructor interval, and no full-process speedup is claimed.

Failed attempts include system-disk ENOSPC, an interrupted queued build, an invalid revision argument, a 192-bit fixture using too many limbs, and incomplete sparse-checkout test fixtures. Each remains recorded in ATTEMPTS.md and the raw logs. Corrected runs use fresh paths. Compiler outputs are on the second SSD; evidence remains in the source worktree.

Comparison SHA256: #text(size: 7pt)[{hash}]

See REPORT.pdf for arithmetic counts; PROTOCOL.md for requirements; INVENTORY.md for the remaining methods; and the JSONL/COMPARISON.json files for raw inputs, binary hashes and exact output matches. CPU matrix measurements are credited only when their actual receipts are available.
"##)).unwrap();
    std::fs::write(root.join("COMPARISON.json"),bytes).unwrap();
}
