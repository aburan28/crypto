  201940:	push   %rbp
  201941:	push   %r15
  201943:	push   %r14
  201945:	push   %r13
  201947:	push   %r12
  201949:	push   %rbx
  20194a:	sub    $0x68,%rsp
  20194e:	mov    %r9,0x18(%rsp)
  201953:	mov    %r8,(%rsp)
  201957:	mov    %rcx,0x20(%rsp)
  20195c:	mov    %rsi,0x8(%rsp)
  201961:	mov    %rdi,0x10(%rsp)
  201966:	mov    (%rdi),%rax
  201969:	mov    %rax,0x28(%rsp)
  20196e:	mov    (%rax),%rax
  201971:	cmp    $0x2,%rax
  201975:	jb     201b54 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x214>
  20197b:	mov    %rdx,%r12
  20197e:	mov    0x10(%rsp),%rcx
  201983:	mov    0x8(%rcx),%r14
  201987:	mov    0x10(%rcx),%rdx
  20198b:	mov    %rdx,0x40(%rsp)
  201990:	mov    0x18(%rcx),%rcx
  201994:	mov    %rcx,0x38(%rsp)
  201999:	mov    0x100(%r12),%rbp
  2019a1:	mov    %r14,0x48(%rsp)
  2019a6:	cs nopw 0x0(%rax,%rax,1)
  2019b0:	mov    %rax,0x30(%rsp)
  2019b5:	lea    -0x1(%rax),%r14
  2019b9:	bsr    %r14,%rcx
  2019bd:	not    %ecx
  2019bf:	mov    %r14,%r15
  2019c2:	shl    %cl,%r15
  2019c5:	dec    %r15
  2019c8:	jmp    2019f2 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xb2>
  2019ca:	nopw   0x0(%rax,%rax,1)
  2019d0:	lea    0x2(%rbp),%rcx
  2019d4:	mov    %rcx,0x100(%r12)
  2019dc:	mov    (%r12,%rbp,4),%rax
  2019e0:	mov    %rcx,%rbp
  2019e3:	mul    %r14
  2019e6:	mov    %rdx,%r13
  2019e9:	cmp    %rax,%r15
  2019ec:	jae    201a8b <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x14b>
  2019f2:	cmp    $0x3f,%rbp
  2019f6:	jb     2019d0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x90>
  2019f8:	je     201a27 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xe7>
  2019fa:	lea    0x110(%r12),%rdi
  201a02:	mov    $0x6,%esi
  201a07:	mov    %r12,%rdx
  201a0a:	call   *0x1c7fd8(%rip)        # 3c99e8 <_DYNAMIC+0xca8>
  201a10:	movq   $0x2,0x100(%r12)
  201a1c:	mov    (%r12),%rax
  201a20:	mov    $0x2,%ebp
  201a25:	jmp    2019e3 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xa3>
  201a27:	mov    0xfc(%r12),%ebx
  201a2f:	lea    0x110(%r12),%rdi
  201a37:	mov    $0x6,%esi
  201a3c:	mov    %r12,%rdx
  201a3f:	call   *0x1c7fa3(%rip)        # 3c99e8 <_DYNAMIC+0xca8>
  201a45:	movq   $0x1,0x100(%r12)
  201a51:	mov    (%r12),%eax
  201a55:	shl    $0x20,%rax
  201a59:	or     %rbx,%rax
  201a5c:	mov    $0x1,%ebp
  201a61:	jmp    2019e3 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xa3>
  201a63:	data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  201a70:	lea    0x2(%rbp),%rcx
  201a74:	mov    %rcx,0x100(%r12)
  201a7c:	mov    (%r12,%rbp,4),%rax
  201a80:	mov    %rcx,%rbp
  201a83:	mul    %r14
  201a86:	cmp    %rax,%r15
  201a89:	jae    201b00 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x1c0>
  201a8b:	cmp    $0x3f,%rbp
  201a8f:	jb     201a70 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x130>
  201a91:	je     201ac0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x180>
  201a93:	lea    0x110(%r12),%rdi
  201a9b:	mov    $0x6,%esi
  201aa0:	mov    %r12,%rdx
  201aa3:	call   *0x1c7f3f(%rip)        # 3c99e8 <_DYNAMIC+0xca8>
  201aa9:	movq   $0x2,0x100(%r12)
  201ab5:	mov    (%r12),%rax
  201ab9:	mov    $0x2,%ebp
  201abe:	jmp    201a83 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x143>
  201ac0:	mov    0xfc(%r12),%ebx
  201ac8:	lea    0x110(%r12),%rdi
  201ad0:	mov    $0x6,%esi
  201ad5:	mov    %r12,%rdx
  201ad8:	call   *0x1c7f0a(%rip)        # 3c99e8 <_DYNAMIC+0xca8>
  201ade:	movq   $0x1,0x100(%r12)
  201aea:	mov    (%r12),%eax
  201aee:	shl    $0x20,%rax
  201af2:	or     %rbx,%rax
  201af5:	mov    $0x1,%ebp
  201afa:	jmp    201a83 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x143>
  201afc:	nopl   0x0(%rax)
  201b00:	mov    %rdx,%rbx
  201b03:	inc    %r13
  201b06:	inc    %rbx
  201b09:	mov    0x40(%rsp),%rax
  201b0e:	mov    (%rax),%rdx
  201b11:	mov    0x8(%rax),%rcx
  201b15:	mov    0x38(%rsp),%rax
  201b1a:	mov    (%rax),%r8
  201b1d:	mov    0x8(%rax),%r9
  201b21:	lea    0x50(%rsp),%rdi
  201b26:	mov    0x48(%rsp),%r14
  201b2b:	mov    %r14,%rsi
  201b2e:	push   %rbx
  201b2f:	push   %r13
  201b31:	call   *0x1c7a91(%rip)        # 3c95c8 <_DYNAMIC+0x888>
  201b37:	add    $0x10,%rsp
  201b3b:	cmpl   $0x1,0x50(%rsp)
  201b40:	je     201b6d <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x22d>
  201b42:	mov    0x28(%rsp),%rax
  201b47:	mov    (%rax),%rax
  201b4a:	cmp    $0x1,%rax
  201b4e:	ja     2019b0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x70>
  201b54:	lea    -0x1e58b0(%rip),%rdi        # 1c2ab <anon.3c11ed13effa8077175ba73b06d1f882.407.llvm.7896473603554930564>
  201b5b:	lea    0x1c0886(%rip),%rdx        # 3c23e8 <anon.3c11ed13effa8077175ba73b06d1f882.409.llvm.7896473603554930564>
  201b62:	mov    $0x19,%esi
  201b67:	call   *0x1c748b(%rip)        # 3c8ff8 <_DYNAMIC+0x2b8>
  201b6d:	mov    0x60(%rsp),%rax
  201b72:	mov    0x10(%rsp),%rcx
  201b77:	mov    0x20(%rcx),%rcx
  201b7b:	cmpb   $0x0,0x28(%rcx)
  201b7f:	je     201bb2 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x272>
  201b81:	cmp    0x40(%r14),%rax
  201b85:	mov    0x8(%rsp),%rdi
  201b8a:	mov    (%rsp),%rsi
  201b8e:	jbe    201ca0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x360>
  201b94:	mov    0x18(%r14),%rdx
  201b98:	sub    %rax,%rdx
  201b9b:	mov    $0x1,%cl
  201b9d:	mov    %rdx,%rax
  201ba0:	mov    0xa0(%rsp),%r8
  201ba8:	cmp    %rsi,%rdi
  201bab:	jb     201bce <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x28e>
  201bad:	jmp    201cb3 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x373>
  201bb2:	xor    %ecx,%ecx
  201bb4:	mov    0x8(%rsp),%rdi
  201bb9:	mov    (%rsp),%rsi
  201bbd:	mov    0xa0(%rsp),%r8
  201bc5:	cmp    %rsi,%rdi
  201bc8:	jae    201cb3 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x373>
  201bce:	mov    0x58(%rsp),%rdx
  201bd3:	mov    0x20(%rsp),%rsi
  201bd8:	mov    %rdx,(%rsi,%rdi,8)
  201bdc:	cmp    %r8,%rdi
  201bdf:	jae    201cc0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x380>
  201be5:	mov    0x18(%rsp),%rdx
  201bea:	mov    %rax,(%rdx,%rdi,8)
  201bee:	mov    0xb0(%rsp),%rsi
  201bf6:	test   %cl,%cl
  201bf8:	je     201c22 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x2e2>
  201bfa:	cmp    %rsi,%rdi
  201bfd:	jae    201cd0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x390>
  201c03:	mov    0x30(%rsp),%rcx
  201c08:	mov    %rcx,%rax
  201c0b:	sub    %r13,%rax
  201c0e:	mov    0xa8(%rsp),%rdx
  201c16:	mov    %rax,(%rdx,%rdi,8)
  201c1a:	sub    %rbx,%rcx
  201c1d:	mov    %rcx,%rbx
  201c20:	jmp    201c37 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x2f7>
  201c22:	cmp    %rsi,%rdi
  201c25:	jae    201cd0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x390>
  201c2b:	mov    0xa8(%rsp),%rax
  201c33:	mov    %r13,(%rax,%rdi,8)
  201c37:	mov    0xc0(%rsp),%rsi
  201c3f:	cmp    %rsi,%rdi
  201c42:	jae    201cdd <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x39d>
  201c48:	mov    0xb8(%rsp),%rax
  201c50:	mov    %rbx,(%rax,%rdi,8)
  201c54:	mov    0xd0(%rsp),%rsi
  201c5c:	cmp    %rsi,%rdi
  201c5f:	jae    201cea <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x3aa>
  201c65:	mov    0xc8(%rsp),%rax
  201c6d:	movq   $0xffffffffffffffff,(%rax,%rdi,8)
  201c75:	mov    0xe0(%rsp),%rsi
  201c7d:	cmp    %rsi,%rdi
  201c80:	jae    201cf7 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x3b7>
  201c82:	mov    0xd8(%rsp),%rax
  201c8a:	movl   $0x0,(%rax,%rdi,4)
  201c91:	add    $0x68,%rsp
  201c95:	pop    %rbx
  201c96:	pop    %r12
  201c98:	pop    %r13
  201c9a:	pop    %r14
  201c9c:	pop    %r15
  201c9e:	pop    %rbp
  201c9f:	ret
  201ca0:	xor    %ecx,%ecx
  201ca2:	mov    0xa0(%rsp),%r8
  201caa:	cmp    %rsi,%rdi
  201cad:	jb     201bce <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x28e>
  201cb3:	lea    0x1bfcbe(%rip),%rdx        # 3c1978 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x68>
  201cba:	call   *0x1c7570(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  201cc0:	lea    0x1bfcc9(%rip),%rdx        # 3c1990 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x80>
  201cc7:	mov    %r8,%rsi
  201cca:	call   *0x1c7560(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  201cd0:	lea    0x1bfcd1(%rip),%rdx        # 3c19a8 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x98>
  201cd7:	call   *0x1c7553(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  201cdd:	lea    0x1bfcdc(%rip),%rdx        # 3c19c0 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0xb0>
  201ce4:	call   *0x1c7546(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  201cea:	lea    0x1bfce7(%rip),%rdx        # 3c19d8 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0xc8>
  201cf1:	call   *0x1c7539(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  201cf7:	lea    0x1bfcf2(%rip),%rdx        # 3c19f0 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0xe0>
  201cfe:	call   *0x1c752c(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
