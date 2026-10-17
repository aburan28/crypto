  1f4d10:	push   %rbp
  1f4d11:	push   %r15
  1f4d13:	push   %r14
  1f4d15:	push   %r13
  1f4d17:	push   %r12
  1f4d19:	push   %rbx
  1f4d1a:	sub    $0x68,%rsp
  1f4d1e:	mov    %r9,0x18(%rsp)
  1f4d23:	mov    %r8,(%rsp)
  1f4d27:	mov    %rcx,0x20(%rsp)
  1f4d2c:	mov    %rsi,0x8(%rsp)
  1f4d31:	mov    %rdi,0x10(%rsp)
  1f4d36:	mov    (%rdi),%rax
  1f4d39:	mov    %rax,0x28(%rsp)
  1f4d3e:	mov    (%rax),%rax
  1f4d41:	cmp    $0x2,%rax
  1f4d45:	jb     1f4f24 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x214>
  1f4d4b:	mov    %rdx,%r12
  1f4d4e:	mov    0x10(%rsp),%rcx
  1f4d53:	mov    0x8(%rcx),%r14
  1f4d57:	mov    0x10(%rcx),%rdx
  1f4d5b:	mov    %rdx,0x40(%rsp)
  1f4d60:	mov    0x18(%rcx),%rcx
  1f4d64:	mov    %rcx,0x38(%rsp)
  1f4d69:	mov    0x100(%r12),%rbp
  1f4d71:	mov    %r14,0x48(%rsp)
  1f4d76:	cs nopw 0x0(%rax,%rax,1)
  1f4d80:	mov    %rax,0x30(%rsp)
  1f4d85:	lea    -0x1(%rax),%r14
  1f4d89:	bsr    %r14,%rcx
  1f4d8d:	not    %ecx
  1f4d8f:	mov    %r14,%r15
  1f4d92:	shl    %cl,%r15
  1f4d95:	dec    %r15
  1f4d98:	jmp    1f4dc2 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xb2>
  1f4d9a:	nopw   0x0(%rax,%rax,1)
  1f4da0:	lea    0x2(%rbp),%rcx
  1f4da4:	mov    %rcx,0x100(%r12)
  1f4dac:	mov    (%r12,%rbp,4),%rax
  1f4db0:	mov    %rcx,%rbp
  1f4db3:	mul    %r14
  1f4db6:	mov    %rdx,%r13
  1f4db9:	cmp    %rax,%r15
  1f4dbc:	jae    1f4e5b <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x14b>
  1f4dc2:	cmp    $0x3f,%rbp
  1f4dc6:	jb     1f4da0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x90>
  1f4dc8:	je     1f4df7 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xe7>
  1f4dca:	lea    0x110(%r12),%rdi
  1f4dd2:	mov    $0x6,%esi
  1f4dd7:	mov    %r12,%rdx
  1f4dda:	call   *0x1d4da0(%rip)        # 3c9b80 <_DYNAMIC+0x950>
  1f4de0:	movq   $0x2,0x100(%r12)
  1f4dec:	mov    (%r12),%rax
  1f4df0:	mov    $0x2,%ebp
  1f4df5:	jmp    1f4db3 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xa3>
  1f4df7:	mov    0xfc(%r12),%ebx
  1f4dff:	lea    0x110(%r12),%rdi
  1f4e07:	mov    $0x6,%esi
  1f4e0c:	mov    %r12,%rdx
  1f4e0f:	call   *0x1d4d6b(%rip)        # 3c9b80 <_DYNAMIC+0x950>
  1f4e15:	movq   $0x1,0x100(%r12)
  1f4e21:	mov    (%r12),%eax
  1f4e25:	shl    $0x20,%rax
  1f4e29:	or     %rbx,%rax
  1f4e2c:	mov    $0x1,%ebp
  1f4e31:	jmp    1f4db3 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0xa3>
  1f4e33:	data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  1f4e40:	lea    0x2(%rbp),%rcx
  1f4e44:	mov    %rcx,0x100(%r12)
  1f4e4c:	mov    (%r12,%rbp,4),%rax
  1f4e50:	mov    %rcx,%rbp
  1f4e53:	mul    %r14
  1f4e56:	cmp    %rax,%r15
  1f4e59:	jae    1f4ed0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x1c0>
  1f4e5b:	cmp    $0x3f,%rbp
  1f4e5f:	jb     1f4e40 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x130>
  1f4e61:	je     1f4e90 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x180>
  1f4e63:	lea    0x110(%r12),%rdi
  1f4e6b:	mov    $0x6,%esi
  1f4e70:	mov    %r12,%rdx
  1f4e73:	call   *0x1d4d07(%rip)        # 3c9b80 <_DYNAMIC+0x950>
  1f4e79:	movq   $0x2,0x100(%r12)
  1f4e85:	mov    (%r12),%rax
  1f4e89:	mov    $0x2,%ebp
  1f4e8e:	jmp    1f4e53 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x143>
  1f4e90:	mov    0xfc(%r12),%ebx
  1f4e98:	lea    0x110(%r12),%rdi
  1f4ea0:	mov    $0x6,%esi
  1f4ea5:	mov    %r12,%rdx
  1f4ea8:	call   *0x1d4cd2(%rip)        # 3c9b80 <_DYNAMIC+0x950>
  1f4eae:	movq   $0x1,0x100(%r12)
  1f4eba:	mov    (%r12),%eax
  1f4ebe:	shl    $0x20,%rax
  1f4ec2:	or     %rbx,%rax
  1f4ec5:	mov    $0x1,%ebp
  1f4eca:	jmp    1f4e53 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x143>
  1f4ecc:	nopl   0x0(%rax)
  1f4ed0:	mov    %rdx,%rbx
  1f4ed3:	inc    %r13
  1f4ed6:	inc    %rbx
  1f4ed9:	mov    0x40(%rsp),%rax
  1f4ede:	mov    (%rax),%rdx
  1f4ee1:	mov    0x8(%rax),%rcx
  1f4ee5:	mov    0x38(%rsp),%rax
  1f4eea:	mov    (%rax),%r8
  1f4eed:	mov    0x8(%rax),%r9
  1f4ef1:	lea    0x50(%rsp),%rdi
  1f4ef6:	mov    0x48(%rsp),%r14
  1f4efb:	mov    %r14,%rsi
  1f4efe:	push   %rbx
  1f4eff:	push   %r13
  1f4f01:	call   *0x1d4bb1(%rip)        # 3c9ab8 <_DYNAMIC+0x888>
  1f4f07:	add    $0x10,%rsp
  1f4f0b:	cmpl   $0x1,0x50(%rsp)
  1f4f10:	je     1f4f3d <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x22d>
  1f4f12:	mov    0x28(%rsp),%rax
  1f4f17:	mov    (%rax),%rax
  1f4f1a:	cmp    $0x1,%rax
  1f4f1e:	ja     1f4d80 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x70>
  1f4f24:	lea    -0x1d8950(%rip),%rdi        # 1c5db <anon.842e8684e504a58b9e8cc66d76e73277.373.llvm.18170108012268971648>
  1f4f2b:	lea    0x1cee5e(%rip),%rdx        # 3c3d90 <anon.842e8684e504a58b9e8cc66d76e73277.375.llvm.18170108012268971648>
  1f4f32:	mov    $0x19,%esi
  1f4f37:	call   *0x1d45ab(%rip)        # 3c94e8 <_DYNAMIC+0x2b8>
  1f4f3d:	mov    0x60(%rsp),%rax
  1f4f42:	mov    0x10(%rsp),%rcx
  1f4f47:	mov    0x20(%rcx),%rcx
  1f4f4b:	cmpb   $0x0,0x28(%rcx)
  1f4f4f:	je     1f4f82 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x272>
  1f4f51:	cmp    0x40(%r14),%rax
  1f4f55:	mov    0x8(%rsp),%rdi
  1f4f5a:	mov    (%rsp),%rsi
  1f4f5e:	jbe    1f5070 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x360>
  1f4f64:	mov    0x18(%r14),%rdx
  1f4f68:	sub    %rax,%rdx
  1f4f6b:	mov    $0x1,%cl
  1f4f6d:	mov    %rdx,%rax
  1f4f70:	mov    0xa0(%rsp),%r8
  1f4f78:	cmp    %rsi,%rdi
  1f4f7b:	jb     1f4f9e <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x28e>
  1f4f7d:	jmp    1f5083 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x373>
  1f4f82:	xor    %ecx,%ecx
  1f4f84:	mov    0x8(%rsp),%rdi
  1f4f89:	mov    (%rsp),%rsi
  1f4f8d:	mov    0xa0(%rsp),%r8
  1f4f95:	cmp    %rsi,%rdi
  1f4f98:	jae    1f5083 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x373>
  1f4f9e:	mov    0x58(%rsp),%rdx
  1f4fa3:	mov    0x20(%rsp),%rsi
  1f4fa8:	mov    %rdx,(%rsi,%rdi,8)
  1f4fac:	cmp    %r8,%rdi
  1f4faf:	jae    1f5090 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x380>
  1f4fb5:	mov    0x18(%rsp),%rdx
  1f4fba:	mov    %rax,(%rdx,%rdi,8)
  1f4fbe:	mov    0xb0(%rsp),%rsi
  1f4fc6:	test   %cl,%cl
  1f4fc8:	je     1f4ff2 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x2e2>
  1f4fca:	cmp    %rsi,%rdi
  1f4fcd:	jae    1f50a0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x390>
  1f4fd3:	mov    0x30(%rsp),%rcx
  1f4fd8:	mov    %rcx,%rax
  1f4fdb:	sub    %r13,%rax
  1f4fde:	mov    0xa8(%rsp),%rdx
  1f4fe6:	mov    %rax,(%rdx,%rdi,8)
  1f4fea:	sub    %rbx,%rcx
  1f4fed:	mov    %rcx,%rbx
  1f4ff0:	jmp    1f5007 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x2f7>
  1f4ff2:	cmp    %rsi,%rdi
  1f4ff5:	jae    1f50a0 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x390>
  1f4ffb:	mov    0xa8(%rsp),%rax
  1f5003:	mov    %r13,(%rax,%rdi,8)
  1f5007:	mov    0xc0(%rsp),%rsi
  1f500f:	cmp    %rsi,%rdi
  1f5012:	jae    1f50ad <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x39d>
  1f5018:	mov    0xb8(%rsp),%rax
  1f5020:	mov    %rbx,(%rax,%rdi,8)
  1f5024:	mov    0xd0(%rsp),%rsi
  1f502c:	cmp    %rsi,%rdi
  1f502f:	jae    1f50ba <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x3aa>
  1f5035:	mov    0xc8(%rsp),%rax
  1f503d:	movq   $0xffffffffffffffff,(%rax,%rdi,8)
  1f5045:	mov    0xe0(%rsp),%rsi
  1f504d:	cmp    %rsi,%rdi
  1f5050:	jae    1f50c7 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x3b7>
  1f5052:	mov    0xd8(%rsp),%rax
  1f505a:	movl   $0x0,(%rax,%rdi,4)
  1f5061:	add    $0x68,%rsp
  1f5065:	pop    %rbx
  1f5066:	pop    %r12
  1f5068:	pop    %r13
  1f506a:	pop    %r14
  1f506c:	pop    %r15
  1f506e:	pop    %rbp
  1f506f:	ret
  1f5070:	xor    %ecx,%ecx
  1f5072:	mov    0xa0(%rsp),%r8
  1f507a:	cmp    %rsi,%rdi
  1f507d:	jb     1f4f9e <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}+0x28e>
  1f5083:	lea    0x1ccd46(%rip),%rdx        # 3c1dd0 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x68>
  1f508a:	call   *0x1d4690(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  1f5090:	lea    0x1ccd51(%rip),%rdx        # 3c1de8 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x80>
  1f5097:	mov    %r8,%rsi
  1f509a:	call   *0x1d4680(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  1f50a0:	lea    0x1ccd59(%rip),%rdx        # 3c1e00 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x98>
  1f50a7:	call   *0x1d4673(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  1f50ad:	lea    0x1ccd64(%rip),%rdx        # 3c1e18 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0xb0>
  1f50b4:	call   *0x1d4666(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  1f50ba:	lea    0x1ccd6f(%rip),%rdx        # 3c1e30 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0xc8>
  1f50c1:	call   *0x1d4659(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  1f50c7:	lea    0x1ccd7a(%rip),%rdx        # 3c1e48 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0xe0>
  1f50ce:	call   *0x1d464c(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
