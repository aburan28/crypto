00000000002c3760 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm>:
  2c3760:	push   %rbp
  2c3761:	push   %r15
  2c3763:	push   %r14
  2c3765:	push   %r13
  2c3767:	push   %r12
  2c3769:	push   %rbx
  2c376a:	sub    $0x58,%rsp
  2c376e:	mov    %rcx,0x50(%rsp)
  2c3773:	mov    %rdi,%r13
  2c3776:	mov    0x10(%rsi),%r15
  2c377a:	test   %r15,%r15
  2c377d:	mov    %rsi,0x8(%rsp)
  2c3782:	mov    %rdx,0x48(%rsp)
  2c3787:	je     2c387a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x11a>
  2c378d:	mov    0x40(%rsi),%rdi
  2c3791:	mov    0x8(%rsi),%rax
  2c3795:	mov    0x38(%rsi),%r9
  2c3799:	mov    0x120(%r13),%r10
  2c37a0:	xor    %ecx,%ecx
  2c37a2:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  2c37b0:	cmp    %rcx,%rdi
  2c37b3:	je     2c3b4c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3ec>
  2c37b9:	mov    (%rax,%rcx,8),%r11
  2c37bd:	mov    %r11d,%r8d
  2c37c0:	and    $0x3e,%r8d
  2c37c4:	mov    (%rdx,%r8,4),%r8
  2c37c8:	sub    %r11,%r8
  2c37cb:	mov    $0x0,%r11d
  2c37d1:	cmovb  %r10,%r11
  2c37d5:	add    %r11,%r8
  2c37d8:	mov    %r8,(%r9,%rcx,8)
  2c37dc:	inc    %rcx
  2c37df:	cmp    %rcx,%r15
  2c37e2:	jne    2c37b0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x50>
  2c37e4:	mov    0x100(%r13),%rcx
  2c37eb:	mov    0x128(%r13),%r11
  2c37f2:	mov    0x50(%rsi),%rbx
  2c37f6:	mov    0x58(%rsi),%r12
  2c37fa:	mov    0x60(%rsi),%r14
  2c37fe:	xor    %ebp,%ebp
  2c3800:	jmp    2c3845 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xe5>
  2c3802:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  2c3810:	mul    %rcx
  2c3813:	mov    %rax,%r8
  2c3816:	mov    %r11,%rax
  2c3819:	imul   %r8,%rax
  2c381d:	mov    %rdx,%rcx
  2c3820:	mul    %r10
  2c3823:	cmp    $0x1,%r8
  2c3827:	sbb    $0xffffffffffffffff,%rcx
  2c382b:	add    %rdx,%rcx
  2c382e:	cmp    %r10,%rcx
  2c3831:	mov    $0x0,%eax
  2c3836:	cmovae %r10,%rax
  2c383a:	sub    %rax,%rcx
  2c383d:	inc    %rbp
  2c3840:	cmp    %rbp,%r15
  2c3843:	je     2c3881 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x121>
  2c3845:	cmp    %rbp,%r12
  2c3848:	je     2c3b5c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3fc>
  2c384e:	mov    %rcx,(%rbx,%rbp,8)
  2c3852:	cmp    %rbp,%rdi
  2c3855:	je     2c3b6f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x40f>
  2c385b:	mov    (%r9,%rbp,8),%rax
  2c385f:	test   %rax,%rax
  2c3862:	jne    2c3810 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xb0>
  2c3864:	inc    %r14
  2c3867:	mov    0x8(%rsp),%rax
  2c386c:	mov    %r14,0x60(%rax)
  2c3870:	inc    %rbp
  2c3873:	cmp    %rbp,%r15
  2c3876:	jne    2c3845 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xe5>
  2c3878:	jmp    2c3881 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x121>
  2c387a:	mov    0x100(%r13),%rcx
  2c3881:	mov    %r13,%rdi
  2c3884:	mov    %rcx,%rsi
  2c3887:	call   *0x1068e3(%rip)        # 3ca170 <_DYNAMIC+0xf40>
  2c388d:	test   $0x1,%al
  2c388f:	je     2c3b34 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3d4>
  2c3895:	mov    0x8(%rsp),%rax
  2c389a:	mov    0x38(%rax),%rcx
  2c389e:	mov    %rcx,0x18(%rsp)
  2c38a3:	mov    0x40(%rax),%rcx
  2c38a7:	mov    0x120(%r13),%rbp
  2c38ae:	mov    0x128(%r13),%rax
  2c38b5:	mov    %rax,0x30(%rsp)
  2c38ba:	lea    -0x1(%r15),%rax
  2c38be:	mov    %r15,%rdi
  2c38c1:	mov    %rcx,0x10(%rsp)
  2c38c6:	cmp    %rcx,%rax
  2c38c9:	jae    2c3992 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x232>
  2c38cf:	mov    0x8(%rsp),%rax
  2c38d4:	mov    0x58(%rax),%rcx
  2c38d8:	mov    0x50(%rax),%rsi
  2c38dc:	xor    %r11d,%r11d
  2c38df:	mov    %r15,%rdi
  2c38e2:	mov    0x30(%rsp),%r14
  2c38e7:	nopw   0x0(%rax,%rax,1)
  2c38f0:	test   %rdi,%rdi
  2c38f3:	je     2c399b <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x23b>
  2c38f9:	mov    0x18(%rsp),%rax
  2c38fe:	mov    -0x8(%rax,%rdi,8),%r9
  2c3903:	dec    %rdi
  2c3906:	test   %r9,%r9
  2c3909:	je     2c38f0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x190>
  2c390b:	cmp    %rcx,%rdi
  2c390e:	jae    2c3bbe <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x45e>
  2c3914:	mov    %rdx,%rax
  2c3917:	mov    %rdx,%rbx
  2c391a:	mulq   (%rsi,%rdi,8)
  2c391e:	mov    %rax,%r10
  2c3921:	mov    %rdx,%r8
  2c3924:	mov    %r14,%rax
  2c3927:	imul   %r10,%rax
  2c392b:	mul    %rbp
  2c392e:	cmp    $0x1,%r10
  2c3932:	sbb    $0xffffffffffffffff,%r8
  2c3936:	add    %rdx,%r8
  2c3939:	cmp    %rbp,%r8
  2c393c:	mov    %rbp,%rax
  2c393f:	cmovb  %r11,%rax
  2c3943:	sub    %rax,%r8
  2c3946:	mov    %r9,%rax
  2c3949:	mul    %rbx
  2c394c:	mov    %rax,%r9
  2c394f:	mov    %r14,%rax
  2c3952:	imul   %r9,%rax
  2c3956:	mov    %rdx,%r10
  2c3959:	mul    %rbp
  2c395c:	mov    %rdx,%rax
  2c395f:	mov    %r10,%rdx
  2c3962:	cmp    $0x1,%r9
  2c3966:	sbb    $0xffffffffffffffff,%rdx
  2c396a:	add    %rax,%rdx
  2c396d:	cmp    %rbp,%rdx
  2c3970:	mov    %rbp,%rax
  2c3973:	cmovb  %r11,%rax
  2c3977:	mov    0x18(%rsp),%r9
  2c397c:	mov    %r8,(%r9,%rdi,8)
  2c3980:	sub    %rax,%rdx
  2c3983:	lea    -0x1(%rdi),%rax
  2c3987:	cmp    0x10(%rsp),%rax
  2c398c:	jb     2c38f0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x190>
  2c3992:	test   %rdi,%rdi
  2c3995:	jne    2c3ba9 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x449>
  2c399b:	test   %r15,%r15
  2c399e:	je     2c3b3d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3dd>
  2c39a4:	mov    0x8(%rsp),%rax
  2c39a9:	mov    0x8(%rax),%rcx
  2c39ad:	mov    %rcx,0x40(%rsp)
  2c39b2:	mov    0x20(%rax),%rcx
  2c39b6:	mov    %rcx,0x38(%rsp)
  2c39bb:	mov    0x28(%rax),%rcx
  2c39bf:	mov    %rcx,0x28(%rsp)
  2c39c4:	mov    0x68(%rax),%rax
  2c39c8:	mov    %rax,0x20(%rsp)
  2c39cd:	xor    %esi,%esi
  2c39cf:	mov    0x40(%rsp),%rbx
  2c39d4:	mov    0x48(%rsp),%r12
  2c39d9:	mov    0x38(%rsp),%r13
  2c39de:	jmp    2c39ec <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x28c>
  2c39e0:	inc    %rsi
  2c39e3:	cmp    %rsi,%r15
  2c39e6:	je     2c3b3d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3dd>
  2c39ec:	cmp    %rsi,0x28(%rsp)
  2c39f1:	je     2c3b7f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x41f>
  2c39f7:	cmp    %rsi,0x10(%rsp)
  2c39fc:	je     2c3b94 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x434>
  2c3a02:	mov    (%rbx,%rsi,8),%rdi
  2c3a06:	mov    %edi,%r11d
  2c3a09:	and    $0x3e,%r11d
  2c3a0d:	mov    0x100(%r12,%r11,4),%rax
  2c3a15:	mov    0x0(%r13,%rsi,8),%r14
  2c3a1a:	sub    %r14,%rax
  2c3a1d:	mov    $0x0,%ecx
  2c3a22:	cmovb  %rbp,%rcx
  2c3a26:	add    %rcx,%rax
  2c3a29:	mov    0x18(%rsp),%rcx
  2c3a2e:	mulq   (%rcx,%rsi,8)
  2c3a32:	mov    %rax,%r8
  2c3a35:	mov    %rdx,%rcx
  2c3a38:	mov    0x30(%rsp),%r10
  2c3a3d:	mov    %r10,%rax
  2c3a40:	imul   %r8,%rax
  2c3a44:	mul    %rbp
  2c3a47:	cmp    $0x1,%r8
  2c3a4b:	sbb    $0xffffffffffffffff,%rcx
  2c3a4f:	add    %rdx,%rcx
  2c3a52:	cmp    %rbp,%rcx
  2c3a55:	mov    %rbp,%rax
  2c3a58:	mov    $0x0,%edx
  2c3a5d:	cmovb  %rdx,%rax
  2c3a61:	sub    %rax,%rcx
  2c3a64:	mov    %rcx,%rax
  2c3a67:	mul    %rcx
  2c3a6a:	mov    %rax,%r8
  2c3a6d:	mov    %rdx,%r9
  2c3a70:	mov    %r10,%rax
  2c3a73:	imul   %r8,%rax
  2c3a77:	mul    %rbp
  2c3a7a:	cmp    $0x1,%r8
  2c3a7e:	sbb    $0xffffffffffffffff,%r9
  2c3a82:	add    %rdx,%r9
  2c3a85:	cmp    %rbp,%r9
  2c3a88:	mov    %rbp,%rax
  2c3a8b:	mov    $0x0,%edx
  2c3a90:	cmovb  %rdx,%rax
  2c3a94:	sub    %rax,%r9
  2c3a97:	sub    %rdi,%r9
  2c3a9a:	mov    $0x0,%eax
  2c3a9f:	cmovb  %rbp,%rax
  2c3aa3:	add    %r9,%rax
  2c3aa6:	sub    (%r12,%r11,4),%rax
  2c3aaa:	mov    $0x0,%edx
  2c3aaf:	cmovb  %rbp,%rdx
  2c3ab3:	lea    (%rdx,%rax,1),%r9
  2c3ab7:	sub    %r9,%rdi
  2c3aba:	mov    $0x0,%eax
  2c3abf:	cmovb  %rbp,%rax
  2c3ac3:	add    %rdi,%rax
  2c3ac6:	mul    %rcx
  2c3ac9:	mov    %rax,%rcx
  2c3acc:	mov    %rdx,%r8
  2c3acf:	imul   %rax,%r10
  2c3ad3:	mov    %r10,%rax
  2c3ad6:	mul    %rbp
  2c3ad9:	cmp    $0x1,%rcx
  2c3add:	sbb    $0xffffffffffffffff,%r8
  2c3ae1:	add    %rdx,%r8
  2c3ae4:	cmp    %rbp,%r8
  2c3ae7:	mov    %rbp,%rax
  2c3aea:	mov    $0x0,%ecx
  2c3aef:	cmovb  %rcx,%rax
  2c3af3:	sub    %rax,%r8
  2c3af6:	sub    %r14,%r8
  2c3af9:	mov    $0x0,%eax
  2c3afe:	cmovb  %rbp,%rax
  2c3b02:	add    %r8,%rax
  2c3b05:	mov    %r9,(%rbx,%rsi,8)
  2c3b09:	mov    %rax,0x0(%r13,%rsi,8)
  2c3b0e:	test   %r9,0x50(%rsp)
  2c3b13:	jne    2c39e0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x280>
  2c3b19:	mov    0x20(%rsp),%rcx
  2c3b1e:	inc    %rcx
  2c3b21:	mov    0x8(%rsp),%rax
  2c3b26:	mov    %rcx,0x20(%rsp)
  2c3b2b:	mov    %rcx,0x68(%rax)
  2c3b2f:	jmp    2c39e0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x280>
  2c3b34:	mov    0x8(%rsp),%rax
  2c3b39:	incq   0x60(%rax)
  2c3b3d:	add    $0x58,%rsp
  2c3b41:	pop    %rbx
  2c3b42:	pop    %r12
  2c3b44:	pop    %r13
  2c3b46:	pop    %r14
  2c3b48:	pop    %r15
  2c3b4a:	pop    %rbp
  2c3b4b:	ret
  2c3b4c:	lea    0x1001f5(%rip),%rdx        # 3c3d48 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0xd8>
  2c3b53:	mov    %rdi,%rsi
  2c3b56:	call   *0x105bc4(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3b5c:	lea    0x1001b5(%rip),%rdx        # 3c3d18 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0xa8>
  2c3b63:	mov    %r12,%rdi
  2c3b66:	mov    %r12,%rsi
  2c3b69:	call   *0x105bb1(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3b6f:	lea    0x1001ba(%rip),%rdx        # 3c3d30 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0xc0>
  2c3b76:	mov    %rdi,%rsi
  2c3b79:	call   *0x105ba1(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3b7f:	lea    0x100132(%rip),%rdx        # 3c3cb8 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x48>
  2c3b86:	mov    0x28(%rsp),%rdi
  2c3b8b:	mov    %rdi,%rsi
  2c3b8e:	call   *0x105b8c(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3b94:	lea    0x100135(%rip),%rdx        # 3c3cd0 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x60>
  2c3b9b:	mov    0x10(%rsp),%rdi
  2c3ba0:	mov    %rdi,%rsi
  2c3ba3:	call   *0x105b77(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3ba9:	dec    %rdi
  2c3bac:	lea    0x100135(%rip),%rdx        # 3c3ce8 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x78>
  2c3bb3:	mov    0x10(%rsp),%rsi
  2c3bb8:	call   *0x105b62(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3bbe:	lea    0x10013b(%rip),%rdx        # 3c3d00 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x90>
  2c3bc5:	mov    %rcx,%rsi
  2c3bc8:	call   *0x105b52(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c3bce:	int3
  2c3bcf:	int3

