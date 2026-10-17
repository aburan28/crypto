#!/usr/bin/env python3
"""Archive retained metadata/report/raw bytes only. No native execution."""
import datetime,hashlib,json,pathlib,tarfile
s=pathlib.Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
payload=s/"native-controls.tar.gz"
assert not payload.exists()
files=[p for p in sorted(s.rglob("*")) if p.is_file() and "custody" not in p.relative_to(s).parts and p.name not in ("native-controls.tar.gz","payload-manifest.json","archive-info.json","package.stdout.log","package.stderr.log")]
assert all(p.stat().st_size<10_000_000 for p in files),"unexpected large payload"
entries=[{"path":str(p.relative_to(s)),"bytes":p.stat().st_size,"sha256":digest(p)} for p in files]
manifest={"schema":"isogeny-trace-payload-manifest/v1","task_id":"TASK-20261016-7b4402","files":entries,"manifest_self_hash":"Recorded externally in archive-info.json to avoid self-reference.","unit_executable_payloads":0,"scientific_runs":0}
m=s/"payload-manifest.json";m.write_text(json.dumps(manifest,indent=2)+"\n")
with tarfile.open(payload,"w:gz") as t:
    for p in files+[m]:t.add(p,arcname=str(p.relative_to(s)),recursive=False)
with tarfile.open(payload) as t:
    members=t.getmembers()
    assert len(members)==len(files)+1
    assert all(x.isfile() for x in members)
    for e in entries:
        b=t.extractfile(e["path"]).read()
        assert len(b)==e["bytes"] and hashlib.sha256(b).hexdigest()==e["sha256"],e["path"]
    assert t.extractfile("payload-manifest.json").read()==m.read_bytes()
info={"schema":"isogeny-trace-archive-info/v1","utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"path":str(payload),"bytes":payload.stat().st_size,"sha256":digest(payload),"regular_members":len(members),"payload_files_excluding_manifest":len(entries),"manifest_sha256":digest(m),"member_verification":"PASS","unit_executable_payloads":0}
(s/"archive-info.json").write_text(json.dumps(info,indent=2)+"\n")
payload.chmod(0o444)
print(json.dumps(info,indent=2))
