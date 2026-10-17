#!/usr/bin/env python3
"""Metadata/hash/archive inspection only. Never launches a test or crypto tool."""
import datetime, hashlib, json, pathlib, platform, re, shutil, subprocess, sys, tarfile
CRYPTO=pathlib.Path("/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015")
HA=pathlib.Path("/Volumes/SSD990-2/crypto-autoresearcher-isogeny-frontier-20261015")
REL="research/isogeny_world_frontier_20261016/replay-optimization-v1"
STUDY=CRYPTO/REL/"interruption-retake-015-v1"
NREL="coordination/goals/GOAL-ISOU-f30a07/next-native-source-20261016/interruption-retake-015-v1"
N=HA/NREL
OUT=pathlib.Path("/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/interruption-retake-015-v1")
SOURCE="015fa7bb504ef63953e0e04e464bb266f4f207bf"
TEST="replay_trace::tests::flushed_begin_survives_bounded_child_interruption"
FROZEN={"dispatch_queue.json":"fbd5ecd235c382527e028b334495c6621ba84bf378a510b4a76008a5e7b1825a","native-command.sh":"98bc11563029aa04712bbd92217529c980a217d0513deb646248c7a02cb18028","PROTOCOL.md":"995a25364158569a136683ab5b1b5f29d8011048941706d3a0a913b003783ffd","bindings.json":"aaf7649bcead2d0253ac11228558cbc210b059d7c5ecc2067483e8670d3eb26d"}
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(b): return hashlib.sha256(b).hexdigest()
def file_record(p):
    b=pathlib.Path(p).read_bytes()
    return {"path":str(p),"bytes":len(b),"sha256":digest(b)}
def git(repo,*args): return subprocess.run(["git","-C",str(repo),*args],check=True,stdout=subprocess.PIPE).stdout
def write(path,obj):
    path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x") as f: json.dump(obj,f,indent=2);f.write("\n")
def scientific():
    base=OUT.parent/"clean-retake-v1"
    return {arm+"/"+name:file_record(base/arm/"install"/name) for arm in ("reference","candidate") for name in ("isogeny","fixed-modulus-bench")}
def prepare():
    assert not OUT.exists(),"fresh native output already exists"
    for name,expected in FROZEN.items():
        b=(N/name).read_bytes();assert digest(b)==expected,name
        assert b==git(HA,"show","1517801dfce137ba966bced3b16dfb50a82ddc61:"+NREL+"/"+name),name
    bindings=json.loads((N/"bindings.json").read_text())
    assert file_record(bindings["gtimeout"]["path"])["sha256"]==bindings["gtimeout"]["sha256"]
    for arm,v in bindings["arms"].items():
        observed=file_record(v["path"]);assert observed["sha256"]==v["sha256"] and observed["bytes"]==v["bytes"],arm
    inp=STUDY/"inputs";inp.mkdir()
    for name in FROZEN: shutil.copyfile(N/name,inp/name)
    trace=git(CRYPTO,"show",SOURCE+":tools/isogeny-cli/src/replay_trace.rs")
    assert trace==(CRYPTO/"tools/isogeny-cli/src/replay_trace.rs").read_bytes()
    (inp/"replay_trace.rs").write_bytes(trace)
    archive=CRYPTO/REL/"clean-retake-015fa7bb504e/native-controls.tar.gz"
    members=["source-commit.txt","source-tree.txt","rustc.txt","host.txt","os.txt","flags.txt","explicit-build-environment.json"]
    for arm in ("reference","candidate"):
        members += [arm+"/"+x for x in ("build-command.txt","feature-arguments.txt","build.exit","build.stderr","tests.exit","tests.stderr","tests.stdout","git-status-before.txt","git-status-after.txt")]
    records=[]
    with tarfile.open(archive) as t:
        for name in members:
            member=t.getmember("./"+name);assert member.isfile()
            b=t.extractfile(member).read()
            p=inp/"original-clean-build"/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
            records.append({"member":member.name,"bytes":len(b),"sha256":digest(b)})
    assert (inp/"original-clean-build/source-commit.txt").read_text().strip()==SOURCE
    claim=N/"claims/TASK-20261016-7b4402.1.claim.json"
    write(STUDY/"provenance-before.json",{
      "schema":"isogeny-trace-retake-provenance/v1","task_id":"TASK-20261016-7b4402","utc":now(),
      "source_commit":SOURCE,"ha_head":git(HA,"rev-parse","HEAD").decode().strip(),
      "crypto_publication_base":git(CRYPTO,"rev-parse","HEAD").decode().strip(),
      "crypto_tracked_status":git(CRYPTO,"status","--porcelain","--untracked-files=no").decode(),
      "protocol_content_commit":"1517801dfce137ba966bced3b16dfb50a82ddc61",
      "frozen_files":{x:file_record(N/x) for x in FROZEN},
      "native_unit_executables":{a:file_record(v["path"]) for a,v in bindings["arms"].items()},
      "gtimeout":file_record(bindings["gtimeout"]["path"]),"scientific_executables":scientific(),
      "source_trace":file_record(inp/"replay_trace.rs"),"claim":file_record(claim),
      "original_archive":{**file_record(archive),"crypto_commit":"c6fb2a332ad42f0be2b30f1c0b0b6e1ec04221f8","git_blob":git(CRYPTO,"rev-parse","c6fb2a332ad42f0be2b30f1c0b0b6e1ec04221f8:"+str(archive.relative_to(CRYPTO))).decode().strip(),"selected_members":records},
      "host":{"uname":list(platform.uname()),"platform":platform.platform()},
      "environment_whitelist":{"LC_ALL":"C","ISOGENY_TRACE_CONTROL_CHILD":"unset by frozen launcher","ISOGENY_REPLAY_TRACE":"unset by frozen launcher","ISOGENY_TRACE_UNIT_EVIDENCE":"fresh per-arm absolute OUT/ARM/retained"},
      "outer_argv":["/bin/sh",str(N/"native-command.sh"),str(OUT)],
      "per_arm_outer_argv":{a:["/opt/homebrew/bin/gtimeout","--kill-after=5s","20s","/usr/bin/time","-l","-o",str(OUT/a/"os-time.txt"),"/usr/bin/env","-u","ISOGENY_TRACE_CONTROL_CHILD","-u","ISOGENY_REPLAY_TRACE","ISOGENY_TRACE_UNIT_EVIDENCE="+str(OUT/a/"retained"),v["path"],"--exact",TEST,"--nocapture","--test-threads=1"] for a,v in bindings["arms"].items()},
      "historical_invocation_hashes_retained":False,"new_invocation_only":True,
      "inference":{"runtime":"native Codex","requested_policy":"executor-implementation","requested_effort":"medium","actual_model":"unknown","actual_effort":"unknown","fallback_used":False},
      "scientific_runs":0,"source_gate_complete":False
    })
    print("PREPARED: frozen metadata, two unit hashes, source trace, and old build inputs bound; no native invocation.")
def verify():
    before=json.loads((STUDY/"provenance-before.json").read_text())
    bindings=json.loads((N/"bindings.json").read_text())
    errors=[];arms={}
    for name,expected in FROZEN.items():
        if digest((N/name).read_bytes())!=expected: errors.append("frozen file changed: "+name)
    for arm,v in bindings["arms"].items():
        d=OUT/arm;issues=[]
        if not d.exists(): arms[arm]={"status":"NOT_STARTED"};errors.append(arm+":not started");continue
        try:
            trace=[json.loads(x) for x in (d/"retained/interrupted.jsonl").read_text().splitlines()]
            control=json.loads((d/"retained/control.json").read_text())
            if len(trace)!=2 or [(x.get("stage"),x.get("event")) for x in trace]!=[("total","begin"),("interruption_fixture_wait","begin")]:issues.append("trace event contract")
            if control.get("cryptographic_work") is not False or control.get("flushed_begin_observed_before_kill") is not True or control.get("unfinished_stage_marked_complete") is not False:issues.append("control contract")
            text=(d/"parent.stdout.log").read_text()
            if TEST+" ... ok" not in text or not re.search(r"test result: ok\. 1 passed; 0 failed; 0 ignored; 0 measured; 128 filtered out",text):issues.append("exact parent test count")
            if (d/"exit-code.txt").read_text().strip()!="0":issues.append("parent exit")
            observed=file_record(v["path"])
            if observed["sha256"]!=v["sha256"] or observed["bytes"]!=v["bytes"]:issues.append("unit executable identity")
            if (d/"executable-before.sha256").read_bytes()!=(d/"executable-after.sha256").read_bytes():issues.append("before/after executable hash")
            time_text=(d/"os-time.txt").read_text()
            match=re.search(r"([0-9.]+)\s+real\s+([0-9.]+)\s+user\s+([0-9.]+)\s+sys",time_text)
            rss=re.search(r"(\d+)\s+maximum resident set size",time_text)
            arms[arm]={"status":"PASS" if not issues else "INCOMPLETE","issues":issues,"parent_exit":int((d/"exit-code.txt").read_text()),"trace_records":trace,"control":control,"child_signal":"not measured","unit_executable_after":observed,"os_receipt":{"real_seconds":float(match[1]) if match else None,"user_seconds":float(match[2]) if match else None,"system_seconds":float(match[3]) if match else None,"max_rss_bytes":int(rss[1]) if rss else None,"interpretation":"OS operational receipt for pure logging, rounded and shared host; no speed or cryptographic memory-admission claim"}}
        except Exception as e: issues.append(type(e).__name__+": "+str(e));arms[arm]={"status":"INCOMPLETE","issues":issues}
        errors.extend([arm+":"+x for x in issues])
    after=scientific()
    if before["scientific_executables"]!=after:errors.append("scientific executable changed")
    if file_record(bindings["gtimeout"]["path"])["sha256"]!=bindings["gtimeout"]["sha256"]:errors.append("gtimeout changed")
    raw=STUDY/"raw";raw.mkdir()
    for p in OUT.rglob("*"):
        if p.is_file():
            target=raw/p.relative_to(OUT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    raw_files=[{**file_record(p),"path":str(p.relative_to(raw))} for p in sorted(raw.rglob("*")) if p.is_file()]
    write(STUDY/"raw-manifest.json",{"schema":"isogeny-trace-raw-manifest/v1","files":raw_files})
    result={"schema":"isogeny-trace-retake-control-check/v1","utc":now(),"status":"PASS" if not errors else "INCOMPLETE","arms":arms,"errors":errors,"scientific_executables_after":after,"raw_files":len(raw_files),"scientific_runs":0,"source_gate_complete":False,"historical_missing_raw_restored":False}
    write(STUDY/"control-check.json",result)
    print(json.dumps({"status":result["status"],"raw_files":len(raw_files),"arm_status":{k:v["status"] for k,v in arms.items()},"errors":errors},indent=2))
{"prepare":prepare,"verify":verify}[sys.argv[1]]()
