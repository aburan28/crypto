#!/usr/bin/env python3
"""Correct the metadata-only test-count parser; never runs native tests."""
import datetime, hashlib, json, pathlib, re
S=pathlib.Path(__file__).resolve().parent
def hashbytes(b):return hashlib.sha256(b).hexdigest()
def rec(p):
    b=pathlib.Path(p).read_bytes();return {"path":str(p),"bytes":len(b),"sha256":hashbytes(b)}
before=json.loads((S/"provenance-before.json").read_text())
old=json.loads((S/"control-check.json").read_text())
assert old["errors"]==["reference:exact parent test count","candidate:exact parent test count"]
manifest=json.loads((S/"raw-manifest.json").read_text())
for r in manifest["files"]:
    b=(S/"raw"/r["path"]).read_bytes()
    assert len(b)==r["bytes"] and hashbytes(b)==r["sha256"],r["path"]
test="replay_trace::tests::flushed_begin_survives_bounded_child_interruption"
arms={};errors=[]
for arm,exe in before["native_unit_executables"].items():
    d=S/"raw"/arm
    text=(d/"parent.stdout.log").read_text()
    selected=re.findall(r"^test (.+?) \.\.\. (.+)$",text,re.M)
    assert selected==[(test,"ok")],selected
    assert text.count("running 1 test")==1
    summaries=re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out",text)
    assert len(summaries)==1,summaries
    passed,failed,ignored,measured,filtered=map(int,summaries[0])
    assert (passed,failed,ignored,measured)==(1,0,0,0)
    assert (d/"exit-code.txt").read_text().strip()=="0"
    records=[json.loads(x) for x in (d/"retained/interrupted.jsonl").read_text().splitlines()]
    assert [(x["stage"],x["event"]) for x in records]==[("total","begin"),("interruption_fixture_wait","begin")]
    control=json.loads((d/"retained/control.json").read_text())
    assert control["cryptographic_work"] is False
    assert control["flushed_begin_observed_before_kill"] is True
    assert control["unfinished_stage_marked_complete"] is False
    assert control["child_exit_code"] is None
    observed=rec(exe["path"]);assert observed==exe
    assert (d/"executable-before.sha256").read_bytes()==(d/"executable-after.sha256").read_bytes()
    assert (d/"executable-before.sha256").read_text().split()[0]==exe["sha256"]
    time_text=(d/"os-time.txt").read_text()
    times=re.search(r"([0-9.]+)\s+real\s+([0-9.]+)\s+user\s+([0-9.]+)\s+sys",time_text)
    rss=re.search(r"(\d+)\s+maximum resident set size",time_text)
    assert times and rss
    arms[arm]={"status":"PASS","parent_test":{"passed":passed,"failed":failed,"ignored":ignored,"measured":measured,"filtered":filtered,"selected":[test]},"parent_exit":0,"trace_records":records,"control":control,"child_signal":"not measured","unit_executable_after":observed,"os_receipt":{"real_seconds":float(times[1]),"user_seconds":float(times[2]),"system_seconds":float(times[3]),"max_rss_bytes":int(rss[1]),"interpretation":"OS operational receipt for pure logging; rounded/shared host; no speed or cryptographic memory-admission claim"}}
for r in before["scientific_executables"].values():assert rec(r["path"])==r
for r in before["frozen_files"].values():assert rec(r["path"])==r
assert rec(before["gtimeout"]["path"])==before["gtimeout"]
out={"schema":"isogeny-trace-retake-control-check/v2","utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"status":"PASS","arms":arms,"errors":[],"raw_files":len(manifest["files"]),"scientific_runs":0,"source_gate_complete":False,"historical_missing_raw_restored":False,"metadata_correction":{"original_check":rec(S/"control-check.json"),"original_inspector":rec(S/"metadata.py"),"reason":"The first metadata checker incorrectly required128filtered from aggregate multi-executable suite totals. This native unit executable selected exactly1namedtest and filtered85others in both arms. The frozen protocol requires exactlyoneparent test; it does not require128filtered.","native_reruns":0,"originals_preserved":True}}
with (S/"control-check-v2.json").open("x") as f:json.dump(out,f,indent=2);f.write("\n")
print(json.dumps({"status":"PASS","raw_files":len(manifest["files"]),"arms":{a:v["parent_test"] for a,v in arms.items()},"new_native_invocations":0},indent=2))
