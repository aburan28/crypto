"""Exact structural 2:4 sparse mapping of the prior radix128 raw128 product."""

if not __debug__:
    raise RuntimeError("This benchmark requires Python assertions; remove -O, -OO and PYTHONOPTIMIZE.")

from pathlib import Path
import hashlib,importlib.util,json,random

import reference as old
prior=Path(__file__).with_name('reference.py')

# k64..112 occur in every row; pair each with an explicit zero. The two
# moving edges are complementary: exactly one of49+t and113+t is in row r.
pairs=[(k,None) for k in range(64,113)]+[(49+t,113+t) for t in range(15)]
columns=[k for pair in pairs for k in pair]
assert all(not (0 <= 112+r-k < 64) for r in range(16) for k in range(32))
assert len(pairs)==64 and len(columns)==128
assert sorted(k for k in columns if k is not None)==list(range(49,128))

def present(r,k):return k is not None and 0<=112+r-k<64

metadata=[]
for r in range(16):
    assert all(sum(present(r,k) for k in pair)==1 for pair in pairs)
    selected=[];row=[]
    for group in range(32):
        ks=columns[group*4:group*4+4]
        idx=[i for i,k in enumerate(ks) if present(r,k)]
        assert len(idx)==2 and idx[0]<idx[1]
        code=idx[0]+4*idx[1]
        assert code in (4,8,9,12,13,14)
        selected += [ks[i] for i in idx]
        row.append(code)
    assert sorted(112+r-k for k in selected)==list(range(64))
    metadata.append(row)
    for j in range(8):
        actual=sorted((112+r-k,k+16*j-112) for k in selected if 0<=k+16*j-112<64)
        expected=[(a,16*j+r-a) for a in range(64) if 0<=16*j+r-a<64]
        assert actual==expected

hybrid_pairs=[(96+t,32+t) for t in range(17)]+[(49+t,113+t) for t in range(15)]
hybrid_columns=[k for pair in hybrid_pairs for k in pair]
assert sorted(hybrid_columns+list(range(64,96)))==list(range(32,128))
hybrid_metadata=[]
for r in range(16):
    row=[];selected=list(range(64,96))
    for group in range(16):
        ks=hybrid_columns[group*4:group*4+4]
        idx=[i for i,k in enumerate(ks) if present(r,k)]
        assert len(idx)==2
        code=idx[0]+4*idx[1]
        assert code in (4,8,9,12,13,14)
        row.append(code);selected += [ks[i] for i in idx]
    assert sorted(112+r-k for k in selected)==list(range(64))
    for j in range(8):
        actual=sorted((112+r-k,k+16*j-112) for k in selected if 0<=k+16*j-112<64)
        assert actual==[(a,16*j+r-a) for a in range(64) if 0<=16*j+r-a<64]
    hybrid_metadata.append(row)

def sparse_dot(a,b):
    aa,bb=old.encode(a),old.encode(b);out=[0]*128
    for r in range(16):
        for j in range(8):
            value=0
            for tile in range(2):
                for group in range(tile*16,(tile+1)*16):
                    nibble=metadata[r][group]
                    for index in (nibble&3,nibble>>2):
                        k=columns[group*4+index]
                        assert k is not None
                        bi=k+16*j-112
                        if 0<=bi<64:value+=aa[112+r-k]*bb[bi]
            out[16*j+r]=value
    return out

def hybrid_dot(a,b):
    aa,bb=old.encode(a),old.encode(b);out=[0]*128
    for r in range(16):
        kept=list(range(64,96))
        for group,nibble in enumerate(hybrid_metadata[r]):
            kept += [hybrid_columns[4*group+(nibble&3)],hybrid_columns[4*group+(nibble>>2)]]
        for j in range(8):
            out[16*j+r]=sum(aa[112+r-k]*bb[k+16*j-112] for k in kept if 0<=k+16*j-112<64)
    return out

rng=random.Random(13124)
edges=[0,1,old.MASK,old.MASK^1,old.MASK//3,2*(old.MASK//3),1<<127,old.MASK^(1<<127)]
tests=[(a,b) for a in edges for b in edges]+[(rng.getrandbits(128),rng.getrandbits(128)) for _ in range(64)]
for a,b in tests:
    c=sparse_dot(a,b)
    assert c==old.raw_convolution(a,b)
    assert hybrid_dot(a,b)==c
    assert old.reconstruct(c,a,b)==old.serial(a,b)
    assert max(c)<=64*129*129

result=dict(valid=True,kind='Exact sparse matrix mapping, no CUDA or performance result',
            sourceMappingSha256=hashlib.sha256(prior.read_bytes()).hexdigest(),
            logicalColumns=columns,metadataNibblesByRow=metadata,
            rows=16,outputsPerRow=8,coefficientInputTerms=64,
            completeSymbolicCoefficientMaps=128,edgeAndDensePairs=len(tests),
            sourceDenseTiles=3,sourceDenseShape='m16n8k32',
            candidateSparseTiles=2,candidateSparseShape='m16n8k64',
            designatedStoredAValuesPerFour=2,metadataDataDependent=False,
            operandPayloadBytes=dict(denseA=1536,denseB=768,sparseA=1024,sparseB=1024,sparseMetadata=256),
            hybrid=dict(denseTiles=1,denseShape='m16n8k32',denseK=list(range(64,96)),
                        sparseTiles=1,sparseShape='m16n8k64',sparseK=hybrid_columns,
                        metadataNibblesByRow=hybrid_metadata,addedPaddingColumns=0,
                        operandPayloadBytes=dict(A=1024,B=768,metadata=128,total=1920),
                        preferredForNextLoweringScreen=True,
                        reason='Keep the original all-active middle K32 slice dense. Pair the17 always-active upper-edge columns with17 original all-zero lower-edge columns, and pair15 complementary edges. This folds two dense slices into one sparse slice without additional padding.'),
            sparseInstructionsPerSecondAt15B=15_000_000_000*(165/32)*2,
            preserved='Every integer convolution coefficient, radix128 carry bound and all-ones correction are unchanged.',
            explanation='49 always-present columns pair with49 zeros;15 complementary edge pairs make64 pairs. Every row has exactly one designated element per pair, hence2per4 after grouping pairs. One global column permutation applies to both A and B.',
            documentation='https://docs.nvidia.com/cuda/parallel-thread-execution/index.html',
            limits=['Matrix-format proof only; CUDA fragment packing, metadata register layout and native instruction lowering are not tested.',
                    'Two sparse instructions versus three dense instructions is not a measured speedup or an equivalence of issue costs.',
                    'Preparation, reconstruction, high-three-bit correction, field reduction, inverse and full walk remain unimplemented for this mapping.',
                    'Payload counts are matrix storage counts, not measured memory traffic. No field/walk throughput or15B claim.'],
            cudaCompiles=0,gpuRuns=0)
destination=Path(__file__).with_name('build')/'mapping-result.json'
destination.parent.mkdir(exist_ok=True)
destination.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('logicalColumns','metadataNibblesByRow')},indent=2))
