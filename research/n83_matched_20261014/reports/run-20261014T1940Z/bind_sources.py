#!/usr/bin/env python3
"""Copy presentation inputs and attach source links; no result recomputation."""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent.parent
COMMIT = "fb26ccaef5c8f5e97aa81b3503d42e5501b030e4"
BASE = f"https://github.com/aburan28/crypto/blob/{COMMIT}/"
PATH_BASE = "research/n83_matched_20261014/"
RECORDS = "records/postprocess-20261015T0018Z/"


def href(label, relative, repository_path=False):
    path = relative if repository_path else PATH_BASE + relative
    return f'<link href="{BASE}{path}" color="#087F83">{label}</link>'


def main():
    data = ROOT / "data"
    data.mkdir(exist_ok=True)
    copied = {}
    for name in ("statistics.json", "statistics_input.json", "presentation_supplement.json",
                 "independent_replay_v3_receipt.json", "collector_extraction_receipt.json",
                 "independent_replay_v3_meter.jsonl"):
        original = STUDY / RECORDS / name
        target = data / name
        shutil.copy2(original, target)
        copied[name] = {"source":str(original.relative_to(STUDY)),
                        "copy":"data/"+name,
                        "sha256":hashlib.sha256(target.read_bytes()).hexdigest()}
    (ROOT/"input_copy_receipt.json").write_text(json.dumps({
        "kind":"report_presentation_input_copies", "native_evidence_commit":COMMIT,
        "research_results_recomputed":False, "inputs":copied}, indent=2)+"\n")
    notes = {
        "0": [("Protocol","PROTOCOL.md"),("Native statistics",RECORDS+"statistics.json"),("Independent replay",RECORDS+"independent_replay_v3_receipt.json"),("Collector seal",RECORDS+"collector_extraction_receipt.json")],
        "1": [("Exact curve registry","docs/curves/registry.json",True),("Frozen protocol","provenance/protocol.json"),("Build/source pins","provenance/build_receipt.json")],
        "2": [("Fixture generator","verification/fixtures/generator_receipt.json"),("Independent fixtures","verification/fixture_replay_receipt.json"),("Frozen execution matrix","provenance/frozen_campaign_manifest.json")],
        "3": [("IC timer source","examples/koblitz_orbit_dlp_fast.rs",True),("Rho timer source","examples/n83_rho_matched_eval.rs",True),("Audit meter",RECORDS+"independent_replay_v3_meter.jsonl")],
        "4": [("All cold rows and native estimator",RECORDS+"statistics.json"),("Raw cold process K01-ic","records/campaign-20261014T1940Z/runs/K01-ic/observation.json")],
        "5": [("All shared cohorts and native estimator",RECORDS+"statistics.json"),("Raw cohort C01-ic","records/campaign-20261014T1940Z/runs/C01-ic/observation.json")],
        "6": [("All native query distributions/intervals",RECORDS+"statistics.json"),("Native statistical source","analysis/src/lib.rs")],
        "7": [("All native target rows",RECORDS+"statistics.json"),("Raw query input",RECORDS+"statistics_input.json")],
        "8": [("All native target rows",RECORDS+"statistics.json"),("Raw query input",RECORDS+"statistics_input.json")],
        "9": [("Exact resource/work supplement",RECORDS+"presentation_supplement.json"),("Native process records",RECORDS+"statistics_input.json")],
        "10": [("Native additive model and finite/nonfinite mass",RECORDS+"statistics.json"),("Frozen native model source","analysis/src/lib.rs")],
        "11": [("Native reservation","records/main-20261014T1940Z/main_isolation.jsonl"),("Native A/A statistics",RECORDS+"statistics.json"),("Build receipt","provenance/build_receipt.json")],
        "12": [("Final independent replay",RECORDS+"independent_replay_v3_receipt.json"),("V3 numeric-copy correction","verification/revisions/20261015-v3/README.md"),("Separate collect-only extraction","postprocessing/collector-only/README.md"),("Pre-timing checks","checks/README.md")],
        "13": [("Executed workload statistics",RECORDS+"statistics.json"),("Accounting protocol","PROTOCOL.md"),("Preserved numeric-copy diagnosis","verification/revisions/20261015-v3/diagnosis.json")]
    }
    meta = json.loads((ROOT/"presentation_metadata.json").read_text())
    meta["section_source_notes"] = {key:"Sources: "+", ".join(href(*item) for item in items)+"." for key,items in notes.items()}
    meta["phase_source_note"] = "Sources: " + href("Native signed phase/residual supplement",RECORDS+"presentation_supplement.json") + ", " + href("Native phase exporter source","postprocessing/src/lib.rs") + ", " + href("Raw producer fields",RECORDS+"statistics_input.json") + "."
    if "Native evidence is pinned at " not in meta["reproduction_details"]:
        meta["reproduction_details"] += " Native evidence is pinned at "+COMMIT+". " + href("Authoritative evidence",RECORDS+"collector_extraction_receipt.json") + "."
    (ROOT/"presentation_metadata.json").write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n")

    links = json.loads((ROOT/"source_links.json").read_text())
    links["evidence_status"] = "NATIVE_RECORDS_VERIFIED_LAYOUT_ONLY"
    links["native_evidence_commit"] = COMMIT
    for key,name in (("M11","statistics.json"),("M12","independent_replay_v3_receipt.json"),("M13","collector_extraction_receipt.json"),("M14","independent_replay_v3_meter.jsonl")):
        relative = RECORDS+name
        links["links"][key]["path"] = "../../"+relative
        links["links"][key]["url"] = BASE+PATH_BASE+relative
    (ROOT/"source_links.json").write_text(json.dumps(links,indent=2)+"\n")
    contract = {
        "contract_status":"COMPLETED_NATIVE_INPUT_BINDING_FOR_LAYOUT_NOT_NEW_EVIDENCE",
        "purpose":"Display copied native results only; the presentation does not rerun statistics or independent verification.",
        "run_id":"run-20261014T1940Z", "native_evidence_commit":COMMIT,
        "measurement_source_freeze_commit":links["measurement_source_freeze_commit"],
        "protocol_sha256":links["protocol_sha256"],
        "frozen_manifest_sha256":links["frozen_manifest_sha256"],
        "inputs":copied,
        "ratio_direction":"IC/rho; below1 favors IC",
        "primary_boundaries":["Actual8pairedsingletonprocesses", "Actual3paired32targetprocesses"],
        "query_diagnostic":"32targetmeans/3sharedpreparations; unequalvalidation/input/outputclocks; secondary",
        "model_status":"N>32unexecuted; nonfinitebootstrapmassretained; rhojumpschargedpertarget",
        "historical_K600_creation_cost":None,
        "memory_limit_hit_detection":None,
        "calibrated_S":None,"ecbench_L2":None,"independent_host":None,
        "posthoc_corrections":"VersionedV3checkerandseparatecollect-onlyextractor; failurespreserved; frozenmeasuredbytesunchanged",
        "research_results_recomputed":False
    }
    (ROOT/"native_report_data_contract.json").write_text(json.dumps(contract,indent=2)+"\n")


if __name__ == "__main__":
    main()
