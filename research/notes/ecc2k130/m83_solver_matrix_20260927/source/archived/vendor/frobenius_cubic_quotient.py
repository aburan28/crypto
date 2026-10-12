#!/usr/bin/env python3
"""Invertible rotation-orbit quotient on normal-basis coordinate words.

For prime p with ord_p(2)=p-1, set K=F2[T]/Phi_p(T).
A p-bit coordinate word c maps linearly to (parity(c), z=c(T) mod Phi_p).
Rotation multiplies z by T. Thus (parity(c), z**p) separates rotation orbits.
The map has Boolean degree <= popcount(p), hence <=3 for p=131.

This auxiliary K has degree p-1. It is NOT the curve's field F_(2^p).
The quotient does not preserve that field's multiplication or elliptic addition.
"""
from dataclasses import dataclass
import math


class CyclotomicField:
    def __init__(self,p):
        if p<3 or any(p%d==0 for d in range(2,math.isqrt(p)+1)):
            raise ValueError('p must be an odd prime')
        order=next(k for k in range(1,p) if pow(2,k,p)==1)
        if order!=p-1:raise ValueError('Phi_p is not irreducible over F2')
        self.p,self.n=p,p-1
        self.modulus=(1<<p)-1
        self.size=1<<self.n
        self.group_order=self.size-1
        self.image_order=self.group_order//p
        if math.gcd(p,self.image_order)!=1:
            raise ValueError('this explicit section requires gcd(p,(2^(p-1)-1)/p)=1')
        self.root_exponent=pow(p,-1,self.image_order)
        basis=[self.mul(1<<j,1<<j) for j in range(self.n)]
        self.square_tables=[]
        for offset in range(0,self.n,8):
            table=[]
            for byte in range(256):
                out=0
                for j in range(8):
                    if offset+j<self.n and ((byte>>j)&1):out^=basis[offset+j]
                table.append(out)
            self.square_tables.append(table)
        assert self.pow(2,p)==1 and 2!=1

    def mul(self,a,b):
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a>>self.n:a^=self.modulus
        return out

    def sq(self,a):
        out=0
        for table in self.square_tables:
            out^=table[a&255];a>>=8
        return out

    def pow(self,a,exponent):
        if exponent==0:return 1
        out=a
        for bit in bin(exponent)[3:]:
            out=self.sq(out)
            if bit=='1':out=self.mul(out,a)
        return out


@dataclass(frozen=True)
class OrbitKey:
    parity:int
    value:int


class FrobeniusQuotient:
    def __init__(self,p=131):
        self.field=CyclotomicField(p)
        self.p=p
        self.mask=(1<<p)-1
        self.index_masks=[sum(1<<i for i in range(p) if (i>>j)&1)
                          for j in range((p-1).bit_length())]
        self.weight_inverses=[0]+[pow(w,-1,p) for w in range(1,p)]

    def rotate(self,word,phase):
        phase%=self.p
        return ((word<<phase)|(word>>(self.p-phase)))&self.mask

    def project(self,word):
        if not 0<=word<=self.mask:raise ValueError('invalid coordinate word')
        z=word^self.mask if word>>(self.p-1) else word
        return word.bit_count()&1,z

    def unproject(self,parity,z):
        if parity not in (0,1) or not 0<=z<self.field.size:
            raise ValueError('invalid split coordinates')
        return z^self.mask if (z.bit_count()&1)!=parity else z

    def encode(self,word):
        parity,z=self.project(word)
        return OrbitKey(parity,self.field.pow(z,self.p))

    def decode(self,key,validate=True):
        """Deterministic representative, using a fixed field exponent, no search.

        With validate=False the caller must supply a key in the encoder image.
        Validation costs another exponentiation; invalid values raise ValueError.
        """
        if key.parity not in (0,1) or not 0<=key.value<self.field.size:
            raise ValueError('invalid orbit key')
        z=self.field.pow(key.value,self.field.root_exponent) if key.value else 0
        if validate and self.field.pow(z,self.p)!=key.value:
            raise ValueError('key is outside quotient image')
        return self.unproject(key.parity,z)

    def canonical_rotation(self,word):
        return min(self.rotate(word,a) for a in range(self.p))

    def centroid_canonical(self,word):
        """Unique rotation with support's first moment zero modulo prime p.

        Returns (canonical_word, phase). Constants have phase zero.
        This is a combinatorial canonicalizer, not a low-degree Boolean map.
        """
        if not 0<=word<=self.mask:raise ValueError('invalid coordinate word')
        weight=word.bit_count()
        if weight in (0,self.p):return word,0
        moment=sum((word&mask).bit_count()<<j for j,mask in enumerate(self.index_masks))
        phase=(-moment*self.weight_inverses[weight])%self.p
        return self.rotate(word,phase),phase

    def canonical_point_words(self,xword,yword):
        """Canonicalize a nonconstant-x curve point in normal coordinates.

        Return ((canonical_x,canonical_y), phase, sign), representing
        sign * Frobenius^phase(P). Negation is (x,y+x) on our binary curve.
        The caller handles infinity and points with x in F2 separately.
        """
        if xword in (0,self.mask):raise ValueError('handle constant-x points separately')
        x,phase=self.centroid_canonical(xword)
        y=self.rotate(yword,phase);neg_y=y^x
        sign=-1 if neg_y<y else 1
        return (x,min(y,neg_y)),phase,sign

    def phase(self,source,destination):
        for a in range(self.p):
            if self.rotate(source,a)==destination:return a
        raise ValueError('words are in different rotation orbits')


class CyclotomicFactorSpace:
    """Union of rotations of an s-dimensional trace-zero subspace.

    Require s properly divides p-1. Basis words come from the multiplicative
    cosets of <2^s> in Z/pZ, with bit zero adjusted for even parity.
    Nonzero parameters are unique: s payload bits and one rotation phase.
    Membership is epsilon=0 and (z^p)^(2^s)=z^p in the auxiliary field.
    For p=131 this gives rank-(130-s) cubic coordinate conditions.
    """
    def __init__(self,p=131,s=26):
        if not 0<s<p-1 or (p-1)%s:raise ValueError('s must properly divide p-1')
        self.q=FrobeniusQuotient(p);self.p,self.s=p,s
        self.multiplier=pow(2,s,p)
        unseen=set(range(1,p));self.cosets=[];self.basis=[]
        while unseen:
            first=min(unseen);orbit=[];v=first
            while v not in orbit:orbit.append(v);v=v*self.multiplier%p
            unseen.difference_update(orbit);self.cosets.append(orbit)
            word=sum(1<<i for i in orbit)
            if len(orbit)%2:word^=1
            self.basis.append(word)
        assert len(self.basis)==s
        for word in self.basis:
            parity,z=self.q.project(word);assert parity==0
            zz=z
            for _ in range(s):zz=self.q.field.sq(zz)
            assert zz==z

    def encode(self,payload,phase=0):
        if not 0<=payload<(1<<self.s):raise ValueError('invalid payload')
        word=0
        for j,basis in enumerate(self.basis):
            if (payload>>j)&1:word^=basis
        return self.q.rotate(word,phase)

    def decode(self,word):
        if word==0:return 0,0
        canonical,phase_to_center=self.q.centroid_canonical(word)
        payload=sum(((canonical>>orbit[0])&1)<<j for j,orbit in enumerate(self.cosets))
        if self.encode(payload)!=canonical:raise ValueError('word not in factor space')
        return payload,(-phase_to_center)%self.p

    def contains_cubic(self,word):
        key=self.q.encode(word)
        if key.parity:return False
        powered=key.value
        for _ in range(self.s):powered=self.q.field.sq(powered)
        return powered==key.value

    def cubic_constraint_masks(self):
        """Independent linear forms on q=z^p; their composition is cubic at 131.

        Each returned mask imposes parity(mask & q)==0. These p-1-s forms
        are selected from the coordinate rows of Frobenius_K^s - identity.
        A separate parity(word)==0 condition is required.
        """
        f=self.q.field;columns=[]
        for j in range(f.n):
            v=1<<j
            for _ in range(self.s):v=f.sq(v)
            columns.append(v^(1<<j))
        rows=[sum(((v>>i)&1)<<j for j,v in enumerate(columns)) for i in range(f.n)]
        pivots={};selected=[]
        for row in rows:
            v=row
            while v:
                i=v.bit_length()-1
                if i in pivots:v^=pivots[i]
                else:pivots[i]=v;selected.append(row);break
        assert len(selected)==f.n-self.s
        return selected
