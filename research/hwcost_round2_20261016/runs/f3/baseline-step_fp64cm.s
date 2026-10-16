0000000000243af0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm>:
  243af0:	push   %rbp
  243af1:	push   %r15
  243af3:	push   %r14
  243af5:	push   %r13
  243af7:	push   %r12
  243af9:	push   %rbx
  243afa:	sub    $0x58,%rsp
  243afe:	mov    %rcx,0x50(%rsp)
  243b03:	mov    %rdi,%r13
  243b06:	mov    0x10(%rsi),%r15
  243b0a:	test   %r15,%r15
  243b0d:	mov    %rsi,0x8(%rsp)
  243b12:	mov    %rdx,0x48(%rsp)
  243b17:	je     243c0a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x11a>
  243b1d:	mov    0x40(%rsi),%rdi
  243b21:	mov    0x8(%rsi),%rax
  243b25:	mov    0x38(%rsi),%r9
  243b29:	mov    0x120(%r13),%r10
  243b30:	xor    %ecx,%ecx
  243b32:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  243b40:	cmp    %rcx,%rdi
  243b43:	je     243edc <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3ec>
  243b49:	mov    (%rax,%rcx,8),%r11
  243b4d:	mov    %r11d,%r8d
  243b50:	and    $0x3e,%r8d
  243b54:	mov    (%rdx,%r8,4),%r8
  243b58:	sub    %r11,%r8
  243b5b:	mov    $0x0,%r11d
  243b61:	cmovb  %r10,%r11
  243b65:	add    %r11,%r8
  243b68:	mov    %r8,(%r9,%rcx,8)
  243b6c:	inc    %rcx
  243b6f:	cmp    %rcx,%r15
  243b72:	jne    243b40 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x50>
  243b74:	mov    0x100(%r13),%rcx
  243b7b:	mov    0x128(%r13),%r11
  243b82:	mov    0x50(%rsi),%rbx
  243b86:	mov    0x58(%rsi),%r12
  243b8a:	mov    0x60(%rsi),%r14
  243b8e:	xor    %ebp,%ebp
  243b90:	jmp    243bd5 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xe5>
  243b92:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  243ba0:	mul    %rcx
  243ba3:	mov    %rax,%r8
  243ba6:	mov    %r11,%rax
  243ba9:	imul   %r8,%rax
  243bad:	mov    %rdx,%rcx
  243bb0:	mul    %r10
  243bb3:	cmp    $0x1,%r8
  243bb7:	sbb    $0xffffffffffffffff,%rcx
  243bbb:	add    %rdx,%rcx
  243bbe:	cmp    %r10,%rcx
  243bc1:	mov    $0x0,%eax
  243bc6:	cmovae %r10,%rax
  243bca:	sub    %rax,%rcx
  243bcd:	inc    %rbp
  243bd0:	cmp    %rbp,%r15
  243bd3:	je     243c11 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x121>
  243bd5:	cmp    %rbp,%r12
  243bd8:	je     243eec <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3fc>
  243bde:	mov    %rcx,(%rbx,%rbp,8)
  243be2:	cmp    %rbp,%rdi
  243be5:	je     243eff <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x40f>
  243beb:	mov    (%r9,%rbp,8),%rax
  243bef:	test   %rax,%rax
  243bf2:	jne    243ba0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xb0>
  243bf4:	inc    %r14
  243bf7:	mov    0x8(%rsp),%rax
  243bfc:	mov    %r14,0x60(%rax)
  243c00:	inc    %rbp
  243c03:	cmp    %rbp,%r15
  243c06:	jne    243bd5 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xe5>
  243c08:	jmp    243c11 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x121>
  243c0a:	mov    0x100(%r13),%rcx
  243c11:	mov    %r13,%rdi
  243c14:	mov    %rcx,%rsi
  243c17:	call   *0x164d2b(%rip)        # 3a8948 <_DYNAMIC+0xf90>
  243c1d:	test   $0x1,%al
  243c1f:	je     243ec4 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3d4>
  243c25:	mov    0x8(%rsp),%rax
  243c2a:	mov    0x38(%rax),%rcx
  243c2e:	mov    %rcx,0x18(%rsp)
  243c33:	mov    0x40(%rax),%rcx
  243c37:	mov    0x120(%r13),%rbp
  243c3e:	mov    0x128(%r13),%rax
  243c45:	mov    %rax,0x30(%rsp)
  243c4a:	lea    -0x1(%r15),%rax
  243c4e:	mov    %r15,%rdi
  243c51:	mov    %rcx,0x10(%rsp)
  243c56:	cmp    %rcx,%rax
  243c59:	jae    243d22 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x232>
  243c5f:	mov    0x8(%rsp),%rax
  243c64:	mov    0x58(%rax),%rcx
  243c68:	mov    0x50(%rax),%rsi
  243c6c:	xor    %r11d,%r11d
  243c6f:	mov    %r15,%rdi
  243c72:	mov    0x30(%rsp),%r14
  243c77:	nopw   0x0(%rax,%rax,1)
  243c80:	test   %rdi,%rdi
  243c83:	je     243d2b <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x23b>
  243c89:	mov    0x18(%rsp),%rax
  243c8e:	mov    -0x8(%rax,%rdi,8),%r9
  243c93:	dec    %rdi
  243c96:	test   %r9,%r9
  243c99:	je     243c80 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x190>
  243c9b:	cmp    %rcx,%rdi
  243c9e:	jae    243f4e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x45e>
  243ca4:	mov    %rdx,%rax
  243ca7:	mov    %rdx,%rbx
  243caa:	mulq   (%rsi,%rdi,8)
  243cae:	mov    %rax,%r10
  243cb1:	mov    %rdx,%r8
  243cb4:	mov    %r14,%rax
  243cb7:	imul   %r10,%rax
  243cbb:	mul    %rbp
  243cbe:	cmp    $0x1,%r10
  243cc2:	sbb    $0xffffffffffffffff,%r8
  243cc6:	add    %rdx,%r8
  243cc9:	cmp    %rbp,%r8
  243ccc:	mov    %rbp,%rax
  243ccf:	cmovb  %r11,%rax
  243cd3:	sub    %rax,%r8
  243cd6:	mov    %r9,%rax
  243cd9:	mul    %rbx
  243cdc:	mov    %rax,%r9
  243cdf:	mov    %r14,%rax
  243ce2:	imul   %r9,%rax
  243ce6:	mov    %rdx,%r10
  243ce9:	mul    %rbp
  243cec:	mov    %rdx,%rax
  243cef:	mov    %r10,%rdx
  243cf2:	cmp    $0x1,%r9
  243cf6:	sbb    $0xffffffffffffffff,%rdx
  243cfa:	add    %rax,%rdx
  243cfd:	cmp    %rbp,%rdx
  243d00:	mov    %rbp,%rax
  243d03:	cmovb  %r11,%rax
  243d07:	mov    0x18(%rsp),%r9
  243d0c:	mov    %r8,(%r9,%rdi,8)
  243d10:	sub    %rax,%rdx
  243d13:	lea    -0x1(%rdi),%rax
  243d17:	cmp    0x10(%rsp),%rax
  243d1c:	jb     243c80 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x190>
  243d22:	test   %rdi,%rdi
  243d25:	jne    243f39 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x449>
  243d2b:	test   %r15,%r15
  243d2e:	je     243ecd <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3dd>
  243d34:	mov    0x8(%rsp),%rax
  243d39:	mov    0x8(%rax),%rcx
  243d3d:	mov    %rcx,0x40(%rsp)
  243d42:	mov    0x20(%rax),%rcx
  243d46:	mov    %rcx,0x38(%rsp)
  243d4b:	mov    0x28(%rax),%rcx
  243d4f:	mov    %rcx,0x28(%rsp)
  243d54:	mov    0x68(%rax),%rax
  243d58:	mov    %rax,0x20(%rsp)
  243d5d:	xor    %esi,%esi
  243d5f:	mov    0x40(%rsp),%rbx
  243d64:	mov    0x48(%rsp),%r12
  243d69:	mov    0x38(%rsp),%r13
  243d6e:	jmp    243d7c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x28c>
  243d70:	inc    %rsi
  243d73:	cmp    %rsi,%r15
  243d76:	je     243ecd <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3dd>
  243d7c:	cmp    %rsi,0x28(%rsp)
  243d81:	je     243f0f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x41f>
  243d87:	cmp    %rsi,0x10(%rsp)
  243d8c:	je     243f24 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x434>
  243d92:	mov    (%rbx,%rsi,8),%rdi
  243d96:	mov    %edi,%r11d
  243d99:	and    $0x3e,%r11d
  243d9d:	mov    0x100(%r12,%r11,4),%rax
  243da5:	mov    0x0(%r13,%rsi,8),%r14
  243daa:	sub    %r14,%rax
  243dad:	mov    $0x0,%ecx
  243db2:	cmovb  %rbp,%rcx
  243db6:	add    %rcx,%rax
  243db9:	mov    0x18(%rsp),%rcx
  243dbe:	mulq   (%rcx,%rsi,8)
  243dc2:	mov    %rax,%r8
  243dc5:	mov    %rdx,%rcx
  243dc8:	mov    0x30(%rsp),%r10
  243dcd:	mov    %r10,%rax
  243dd0:	imul   %r8,%rax
  243dd4:	mul    %rbp
  243dd7:	cmp    $0x1,%r8
  243ddb:	sbb    $0xffffffffffffffff,%rcx
  243ddf:	add    %rdx,%rcx
  243de2:	cmp    %rbp,%rcx
  243de5:	mov    %rbp,%rax
  243de8:	mov    $0x0,%edx
  243ded:	cmovb  %rdx,%rax
  243df1:	sub    %rax,%rcx
  243df4:	mov    %rcx,%rax
  243df7:	mul    %rcx
  243dfa:	mov    %rax,%r8
  243dfd:	mov    %rdx,%r9
  243e00:	mov    %r10,%rax
  243e03:	imul   %r8,%rax
  243e07:	mul    %rbp
  243e0a:	cmp    $0x1,%r8
  243e0e:	sbb    $0xffffffffffffffff,%r9
  243e12:	add    %rdx,%r9
  243e15:	cmp    %rbp,%r9
  243e18:	mov    %rbp,%rax
  243e1b:	mov    $0x0,%edx
  243e20:	cmovb  %rdx,%rax
  243e24:	sub    %rax,%r9
  243e27:	sub    %rdi,%r9
  243e2a:	mov    $0x0,%eax
  243e2f:	cmovb  %rbp,%rax
  243e33:	add    %r9,%rax
  243e36:	sub    (%r12,%r11,4),%rax
  243e3a:	mov    $0x0,%edx
  243e3f:	cmovb  %rbp,%rdx
  243e43:	lea    (%rdx,%rax,1),%r9
  243e47:	sub    %r9,%rdi
  243e4a:	mov    $0x0,%eax
  243e4f:	cmovb  %rbp,%rax
  243e53:	add    %rdi,%rax
  243e56:	mul    %rcx
  243e59:	mov    %rax,%rcx
  243e5c:	mov    %rdx,%r8
  243e5f:	imul   %rax,%r10
  243e63:	mov    %r10,%rax
  243e66:	mul    %rbp
  243e69:	cmp    $0x1,%rcx
  243e6d:	sbb    $0xffffffffffffffff,%r8
  243e71:	add    %rdx,%r8
  243e74:	cmp    %rbp,%r8
  243e77:	mov    %rbp,%rax
  243e7a:	mov    $0x0,%ecx
  243e7f:	cmovb  %rcx,%rax
  243e83:	sub    %rax,%r8
  243e86:	sub    %r14,%r8
  243e89:	mov    $0x0,%eax
  243e8e:	cmovb  %rbp,%rax
  243e92:	add    %r8,%rax
  243e95:	mov    %r9,(%rbx,%rsi,8)
  243e99:	mov    %rax,0x0(%r13,%rsi,8)
  243e9e:	test   %r9,0x50(%rsp)
  243ea3:	jne    243d70 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x280>
  243ea9:	mov    0x20(%rsp),%rcx
  243eae:	inc    %rcx
  243eb1:	mov    0x8(%rsp),%rax
  243eb6:	mov    %rcx,0x20(%rsp)
  243ebb:	mov    %rcx,0x68(%rax)
  243ebf:	jmp    243d70 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x280>
  243ec4:	mov    0x8(%rsp),%rax
  243ec9:	incq   0x60(%rax)
  243ecd:	add    $0x58,%rsp
  243ed1:	pop    %rbx
  243ed2:	pop    %r12
  243ed4:	pop    %r13
  243ed6:	pop    %r14
  243ed8:	pop    %r15
  243eda:	pop    %rbp
  243edb:	ret
  243edc:	lea    0x15d715(%rip),%rdx        # 3a15f8 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x168>
  243ee3:	mov    %rdi,%rsi
  243ee6:	call   *0x163fb4(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243eec:	lea    0x15d6d5(%rip),%rdx        # 3a15c8 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x138>
  243ef3:	mov    %r12,%rdi
  243ef6:	mov    %r12,%rsi
  243ef9:	call   *0x163fa1(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243eff:	lea    0x15d6da(%rip),%rdx        # 3a15e0 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x150>
  243f06:	mov    %rdi,%rsi
  243f09:	call   *0x163f91(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243f0f:	lea    0x15d63a(%rip),%rdx        # 3a1550 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0xc0>
  243f16:	mov    0x28(%rsp),%rdi
  243f1b:	mov    %rdi,%rsi
  243f1e:	call   *0x163f7c(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243f24:	lea    0x15d63d(%rip),%rdx        # 3a1568 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0xd8>
  243f2b:	mov    0x10(%rsp),%rdi
  243f30:	mov    %rdi,%rsi
  243f33:	call   *0x163f67(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243f39:	dec    %rdi
  243f3c:	lea    0x15d655(%rip),%rdx        # 3a1598 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x108>
  243f43:	mov    0x10(%rsp),%rsi
  243f48:	call   *0x163f52(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243f4e:	lea    0x15d65b(%rip),%rdx        # 3a15b0 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x120>
  243f55:	mov    %rcx,%rsi
  243f58:	call   *0x163f42(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243f5e:	int3
  243f5f:	int3

