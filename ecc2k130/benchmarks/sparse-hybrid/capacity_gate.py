"""Bind a fresh compile only when its source, native code and resources match the reviewed fixture."""
from pathlib import Path
import hashlib,json,re

HERE=Path(__file__).resolve().parent
NAMES={'denseKernel','hybridKernel','dense4','hybrid4','dense8','hybrid8'}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def code(text):
    result={}
    for section in text.split('Function : ')[1:]:
        name=section.splitlines()[0].strip()
        if name not in NAMES:continue
        assert name not in result,'duplicate kernel section'
        words=[];pending=None
        for line in section.splitlines():
            m=re.match(r'\s*/\*([0-9a-f]+)\*/.*?;\s*/\* (0x[0-9a-f]{16}) \*/',line)
            if m:
                assert pending is None
                pending=[int(m[1],16),m[2]]
            else:
                m=re.match(r'\s*/\* (0x[0-9a-f]{16}) \*/',line)
                if m and pending is not None:words.append(pending+[m[1]]);pending=None
        assert words and pending is None
        assert [x[0] for x in words]==list(range(0,16*len(words),16))
        result[name]=words
    assert set(result)==NAMES,'missing reviewed kernel'
    return result

def resources(text):
    out={}
    for name in NAMES:
        rows=re.findall(r'Function '+name+r':\s*REG:(\d+) STACK:(\d+) SHARED:(\d+) LOCAL:(\d+)',text)
        assert len(rows)==1,'missing/duplicate resources'
        out[name]=tuple(map(int,rows[0]))
    return out

def bind_review(compiled_path,output_path):
    original_path=HERE/'evidence/capacity-compile-result.json'
    original=json.loads(original_path.read_text())
    review_path=HERE/'evidence/capacity-code-review.json'
    review=json.loads(review_path.read_text())
    actual=json.loads(Path(compiled_path).read_text())
    assert review['valid'] and sha(original_path)==review['compileRawSha256']
    assert actual['valid'] and actual['compileAttempts']==1 and actual['image']==review['image']
    assert actual['sourceHashes']==actual['expectedSourceHashes']==review['sourceHashes']
    for name,digest in actual['sourceHashes'].items():assert sha(HERE/name)==digest,'source changed since review'
    commands={r['label']:r for r in actual['commands']}
    old_commands={r['label']:r for r in original['commands']}
    assert set(commands)==set(old_commands)
    for name,row in commands.items():
        assert row['returncode']==0 and row['timedOut'] is False
        assert row['command']==old_commands[name]['command'],'compiler/inspection command changed'
    assert all('V13.3.73' in commands[n]['output'] for n in ('nvcc version','ptxas version','cuobjdump version'))
    assert code(commands['SASS']['output'])==code(old_commands['SASS']['output']),'native code differs; review before timing'
    assert resources(commands['resources']['output'])==resources(old_commands['resources']['output']),'resources differ'
    bound=dict(valid=True,kind='Exact reproduction of reviewed source/native code/resources',
               originalCodeReviewSha256=sha(review_path),originalCompileSha256=sha(original_path),
               compileRawSha256=sha(compiled_path),binarySha256=actual['binarySha256'],
               sourceHashes=actual['sourceHashes'],kernels=review['kernels'],
               limits=['This reuses the original semantic/code review through exact source and native encoding equality; it is not a fresh performance result.'])
    Path(output_path).write_text(json.dumps(bound,indent=2)+'\n')
    return bound
