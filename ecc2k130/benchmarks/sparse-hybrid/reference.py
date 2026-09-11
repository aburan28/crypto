"""Integer convolution and independent carryless reference; raw128 only."""
MASK=(1<<128)-1

def encode(a):
    assert 0<=a<=MASK
    return [((a>>(2*i))&1)+128*((a>>(2*i+1))&1) for i in range(64)]

def serial(a,b):
    out=0
    while b:
        bit=b&-b;out^=a<<(bit.bit_length()-1);b^=bit
    return out

def raw_convolution(a,b):
    aa,bb=encode(a),encode(b);c=[0]*128
    # Sparse enumeration makes the complete basis panel inexpensive.
    for i,x in enumerate(aa):
        if x:
            for j,y in enumerate(bb):
                if y:c[i+j]+=x*y
    return c

def reconstruct(c,a,b,correct=True):
    assert len(c)==128 and c[127]==0
    out=0
    for l,value in enumerate(c):
        out^=(value&1)<<(2*l)
        out^=((value>>7)&1)<<(2*l+1)
        out^=((value>>14)&1)<<(2*l+2)
    if correct and a==MASK and b==MASK:out^=1<<128
    return out
