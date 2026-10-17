#!/usr/bin/env python3
"""Metadata-only full archive/readback and payload verification."""
import datetime,hashlib,json,pathlib,tarfile
s=pathlib.Path(__file__).resolve().parent
a=json.loads((s/"archive-info.json").read_text())
p=s/"native-controls.tar.gz";b=p.read_bytes()
assert len(b)==a["bytes"] and hashlib.sha256(b).hexdigest()==a["sha256"]
cu=s/"custody";assert (cu/"transfer.exit-code.txt").read_text().strip()=="0"
idx=json.loads((cu/"store-index.json").read_text())
logical="evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/interruption-retake-v1"
assert list(idx["entries"])==[logical]
e=idx["entries"][logical];assert e["sha256"]==a["sha256"] and e["bytes"]==a["bytes"]
push=json.loads((cu/"push.json").read_text())
assert push["uploaded"] is True and push["replaced"] is False and push["cached"] is False
assert push["sha256"]==a["sha256"] and push["bytes"]==a["bytes"]
path=pathlib.Path((cu/"pull.stdout.log").read_text().strip());rb=path.read_bytes()
assert rb==b and hashlib.sha256(rb).hexdigest()==a["sha256"]
assert (cu/"verify.stdout.log").read_text().strip()==logical+": cache ok, remote ok"
assert (cu/"store-client-before.sha256").read_bytes()==(cu/"store-client-after.sha256").read_bytes()
manifest=json.loads((s/"payload-manifest.json").read_text())
with tarfile.open(p) as t:
    ms=t.getmembers();assert len(ms)==a["regular_members"]
    for x in manifest["files"]:
        ob=(s/x["path"]).read_bytes();ab=t.extractfile(x["path"]).read()
        assert ob==ab and len(ab)==x["bytes"] and hashlib.sha256(ab).hexdigest()==x["sha256"],x["path"]
    assert t.extractfile("payload-manifest.json").read()==(s/"payload-manifest.json").read_bytes()
out={"schema":"isogeny-trace-custody-verification/v1","utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"status":"VERIFIED_REMOTE_AND_FULL_READBACK","logical_names":1,"distinct_objects":1,"logical":logical,"remote":"s3://crypto-autoresearcher/store","archive":{"sha256":a["sha256"],"bytes":a["bytes"],"regular_members":a["regular_members"]},"readback_path":str(path),"native_remote_verification":"cache ok, remote ok","independent_sha256":"PASS","full_byte_equality":"PASS","payload_members_match_archive_and_worktree":"PASS","native_store_client_sha256":(cu/"store-client-before.sha256").read_text().split()[0],"store_client_build_source":"not independently established","unit_executable_payloads":0,"scientific_runs":0,"source_gate_complete":False,"historical_raw_restored":False}
with (cu/"verification.json").open("x") as f:json.dump(out,f,indent=2);f.write("\n")
print(json.dumps(out,indent=2))
