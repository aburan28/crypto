0000000000244630 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64>:
  244630:	push   %rbp
  244631:	push   %r15
  244633:	push   %r14
  244635:	push   %r13
  244637:	push   %r12
  244639:	push   %rbx
  24463a:	sub    $0x38,%rsp
  24463e:	mov    %rcx,0x30(%rsp)
  244643:	mov    %rdx,%r15
  244646:	mov    %rdi,%rbx
  244649:	mov    0x10(%rsi),%r13
  24464d:	test   %r13,%r13
  244650:	mov    %rsi,0x10(%rsp)
  244655:	mov    %rdi,0x8(%rsp)
  24465a:	je     24474a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x11a>
  244660:	mov    0x40(%rsi),%rdi
  244664:	mov    0x8(%rsi),%rax
  244668:	mov    0x38(%rsi),%r9
  24466c:	mov    0xf0(%rbx),%r10
  244673:	xor    %ecx,%ecx
  244675:	data16 cs nopw 0x0(%rax,%rax,1)
  244680:	cmp    %rcx,%rdi
  244683:	je     244a2c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3fc>
  244689:	mov    (%rax,%rcx,8),%rdx
  24468d:	mov    %edx,%r8d
  244690:	and    $0x3e,%r8d
  244694:	mov    (%r15,%r8,4),%r8
  244698:	sub    %rdx,%r8
  24469b:	mov    $0x0,%edx
  2446a0:	cmovb  %r10,%rdx
  2446a4:	add    %r8,%rdx
  2446a7:	mov    %rdx,(%r9,%rcx,8)
  2446ab:	inc    %rcx
  2446ae:	cmp    %rcx,%r13
  2446b1:	jne    244680 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x50>
  2446b3:	mov    0xf8(%rbx),%r11
  2446ba:	mov    0x100(%rbx),%rcx
  2446c1:	mov    0x50(%rsi),%rbx
  2446c5:	mov    0x58(%rsi),%r12
  2446c9:	mov    0x60(%rsi),%r14
  2446cd:	xor    %ebp,%ebp
  2446cf:	jmp    244715 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xe5>
  2446d1:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  2446e0:	mul    %rcx
  2446e3:	mov    %rax,%r8
  2446e6:	mov    %r11,%rax
  2446e9:	imul   %r8,%rax
  2446ed:	mov    %rdx,%rcx
  2446f0:	mul    %r10
  2446f3:	cmp    $0x1,%r8
  2446f7:	sbb    $0xffffffffffffffff,%rcx
  2446fb:	add    %rdx,%rcx
  2446fe:	cmp    %r10,%rcx
  244701:	mov    $0x0,%eax
  244706:	cmovae %r10,%rax
  24470a:	sub    %rax,%rcx
  24470d:	inc    %rbp
  244710:	cmp    %rbp,%r13
  244713:	je     244751 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x121>
  244715:	cmp    %rbp,%r12
  244718:	je     244a3c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x40c>
  24471e:	mov    %rcx,(%rbx,%rbp,8)
  244722:	cmp    %rbp,%rdi
  244725:	je     244a4f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x41f>
  24472b:	mov    (%r9,%rbp,8),%rax
  24472f:	test   %rax,%rax
  244732:	jne    2446e0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xb0>
  244734:	inc    %r14
  244737:	mov    0x10(%rsp),%rax
  24473c:	mov    %r14,0x60(%rax)
  244740:	inc    %rbp
  244743:	cmp    %rbp,%r13
  244746:	jne    244715 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xe5>
  244748:	jmp    244751 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x121>
  24474a:	mov    0x100(%rbx),%rcx
  244751:	mov    0x8(%rsp),%rbx
  244756:	mov    %rbx,%rdi
  244759:	mov    %rcx,%rsi
  24475c:	call   *0x1641e6(%rip)        # 3a8948 <_DYNAMIC+0xf90>
  244762:	test   $0x1,%al
  244764:	je     244a14 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3e4>
  24476a:	mov    0x10(%rsp),%rsi
  24476f:	mov    0x38(%rsi),%rax
  244773:	mov    %rax,0x18(%rsp)
  244778:	mov    0x40(%rsi),%rcx
  24477c:	mov    0xf0(%rbx),%rbp
  244783:	mov    0xf8(%rbx),%r12
  24478a:	lea    -0x1(%r13),%rax
  24478e:	mov    %r13,%rdi
  244791:	mov    %rcx,0x8(%rsp)
  244796:	cmp    %rcx,%rax
  244799:	jae    244852 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x222>
  24479f:	mov    0x58(%rsi),%rcx
  2447a3:	mov    0x50(%rsi),%rsi
  2447a7:	xor    %r11d,%r11d
  2447aa:	mov    %r13,%rdi
  2447ad:	nopl   (%rax)
  2447b0:	test   %rdi,%rdi
  2447b3:	je     24485b <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x22b>
  2447b9:	mov    0x18(%rsp),%rax
  2447be:	mov    -0x8(%rax,%rdi,8),%r9
  2447c3:	dec    %rdi
  2447c6:	test   %r9,%r9
  2447c9:	je     2447b0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x180>
  2447cb:	cmp    %rcx,%rdi
  2447ce:	jae    244a9e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x46e>
  2447d4:	mov    %rdx,%rax
  2447d7:	mov    %rdx,%rbx
  2447da:	mulq   (%rsi,%rdi,8)
  2447de:	mov    %rax,%r10
  2447e1:	mov    %rdx,%r8
  2447e4:	mov    %r12,%rax
  2447e7:	imul   %r10,%rax
  2447eb:	mul    %rbp
  2447ee:	cmp    $0x1,%r10
  2447f2:	sbb    $0xffffffffffffffff,%r8
  2447f6:	add    %rdx,%r8
  2447f9:	cmp    %rbp,%r8
  2447fc:	mov    %rbp,%rax
  2447ff:	cmovb  %r11,%rax
  244803:	sub    %rax,%r8
  244806:	mov    %r9,%rax
  244809:	mul    %rbx
  24480c:	mov    %rax,%r9
  24480f:	mov    %r12,%rax
  244812:	imul   %r9,%rax
  244816:	mov    %rdx,%r10
  244819:	mul    %rbp
  24481c:	mov    %rdx,%rax
  24481f:	mov    %r10,%rdx
  244822:	cmp    $0x1,%r9
  244826:	sbb    $0xffffffffffffffff,%rdx
  24482a:	add    %rax,%rdx
  24482d:	cmp    %rbp,%rdx
  244830:	mov    %rbp,%rax
  244833:	cmovb  %r11,%rax
  244837:	mov    0x18(%rsp),%r9
  24483c:	mov    %r8,(%r9,%rdi,8)
  244840:	sub    %rax,%rdx
  244843:	lea    -0x1(%rdi),%rax
  244847:	cmp    0x8(%rsp),%rax
  24484c:	jb     2447b0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x180>
  244852:	test   %rdi,%rdi
  244855:	jne    244a89 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x459>
  24485b:	test   %r13,%r13
  24485e:	je     244a1d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ed>
  244864:	mov    0x10(%rsp),%rax
  244869:	mov    0x8(%rax),%r10
  24486d:	mov    0x20(%rax),%rbx
  244871:	mov    0x28(%rax),%rcx
  244875:	mov    %rcx,0x28(%rsp)
  24487a:	mov    0x68(%rax),%rax
  24487e:	mov    %rax,0x20(%rsp)
  244883:	xor    %r14d,%r14d
  244886:	jmp    2448b2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x282>
  244888:	nopl   0x0(%rax,%rax,1)
  244890:	add    %rcx,%rax
  244893:	mov    %r9,(%r10,%r14,8)
  244897:	mov    %rax,(%rbx,%r14,8)
  24489b:	test   %r9,0x30(%rsp)
  2448a0:	je     2449f0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3c0>
  2448a6:	inc    %r14
  2448a9:	cmp    %r14,%r13
  2448ac:	je     244a1d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ed>
  2448b2:	cmp    %r14,0x28(%rsp)
  2448b7:	je     244a5f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x42f>
  2448bd:	cmp    %r14,0x8(%rsp)
  2448c2:	je     244a74 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x444>
  2448c8:	mov    (%r10,%r14,8),%rdi
  2448cc:	mov    %edi,%r11d
  2448cf:	and    $0x3e,%r11d
  2448d3:	mov    0x100(%r15,%r11,4),%rax
  2448db:	mov    (%rbx,%r14,8),%rsi
  2448df:	mov    %rbp,%rcx
  2448e2:	sub    %rsi,%rax
  2448e5:	jb     2448e9 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x2b9>
  2448e7:	xor    %ecx,%ecx
  2448e9:	add    %rcx,%rax
  2448ec:	mov    0x18(%rsp),%rcx
  2448f1:	mulq   (%rcx,%r14,8)
  2448f5:	mov    %rax,%r8
  2448f8:	mov    %rdx,%rcx
  2448fb:	mov    %r12,%rax
  2448fe:	imul   %r8,%rax
  244902:	mul    %rbp
  244905:	cmp    $0x1,%r8
  244909:	sbb    $0xffffffffffffffff,%rcx
  24490d:	add    %rdx,%rcx
  244910:	mov    $0x0,%eax
  244915:	cmp    %rbp,%rcx
  244918:	jb     24491d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x2ed>
  24491a:	mov    %rbp,%rax
  24491d:	sub    %rax,%rcx
  244920:	mov    %rcx,%rax
  244923:	mul    %rcx
  244926:	mov    %rax,%r9
  244929:	mov    %rdx,%r8
  24492c:	mov    %r12,%rax
  24492f:	imul   %r9,%rax
  244933:	mul    %rbp
  244936:	cmp    $0x1,%r9
  24493a:	sbb    $0xffffffffffffffff,%r8
  24493e:	add    %rdx,%r8
  244941:	mov    $0x0,%eax
  244946:	cmp    %rbp,%r8
  244949:	jae    2449c0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x390>
  24494b:	sub    %rax,%r8
  24494e:	mov    %rbp,%rax
  244951:	sub    %rdi,%r8
  244954:	jae    2449ce <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x39e>
  244956:	add    %r8,%rax
  244959:	mov    %rbp,%rdx
  24495c:	sub    (%r15,%r11,4),%rax
  244960:	jae    2449dc <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ac>
  244962:	lea    (%rdx,%rax,1),%r9
  244966:	mov    %rbp,%rax
  244969:	sub    %r9,%rdi
  24496c:	jb     244970 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x340>
  24496e:	xor    %eax,%eax
  244970:	add    %rdi,%rax
  244973:	mul    %rcx
  244976:	mov    %rax,%r8
  244979:	mov    %rdx,%rcx
  24497c:	mov    %r12,%rax
  24497f:	imul   %r8,%rax
  244983:	mul    %rbp
  244986:	cmp    $0x1,%r8
  24498a:	sbb    $0xffffffffffffffff,%rcx
  24498e:	add    %rdx,%rcx
  244991:	mov    $0x0,%eax
  244996:	cmp    %rbp,%rcx
  244999:	jb     24499e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x36e>
  24499b:	mov    %rbp,%rax
  24499e:	sub    %rax,%rcx
  2449a1:	mov    %rbp,%rax
  2449a4:	sub    %rsi,%rcx
  2449a7:	jb     244890 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x260>
  2449ad:	xor    %eax,%eax
  2449af:	jmp    244890 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x260>
  2449b4:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  2449c0:	mov    %rbp,%rax
  2449c3:	sub    %rax,%r8
  2449c6:	mov    %rbp,%rax
  2449c9:	sub    %rdi,%r8
  2449cc:	jb     244956 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x326>
  2449ce:	xor    %eax,%eax
  2449d0:	add    %r8,%rax
  2449d3:	mov    %rbp,%rdx
  2449d6:	sub    (%r15,%r11,4),%rax
  2449da:	jb     244962 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x332>
  2449dc:	xor    %edx,%edx
  2449de:	lea    (%rdx,%rax,1),%r9
  2449e2:	mov    %rbp,%rax
  2449e5:	sub    %r9,%rdi
  2449e8:	jae    24496e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x33e>
  2449ea:	jmp    244970 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x340>
  2449ec:	nopl   0x0(%rax)
  2449f0:	mov    0x20(%rsp),%rcx
  2449f5:	inc    %rcx
  2449f8:	mov    0x10(%rsp),%rax
  2449fd:	mov    %rcx,0x20(%rsp)
  244a02:	mov    %rcx,0x68(%rax)
  244a06:	inc    %r14
  244a09:	cmp    %r14,%r13
  244a0c:	jne    2448b2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x282>
  244a12:	jmp    244a1d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ed>
  244a14:	mov    0x10(%rsp),%rax
  244a19:	incq   0x60(%rax)
  244a1d:	add    $0x38,%rsp
  244a21:	pop    %rbx
  244a22:	pop    %r12
  244a24:	pop    %r13
  244a26:	pop    %r14
  244a28:	pop    %r15
  244a2a:	pop    %rbp
  244a2b:	ret
  244a2c:	lea    0x15cbc5(%rip),%rdx        # 3a15f8 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x168>
  244a33:	mov    %rdi,%rsi
  244a36:	call   *0x163464(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244a3c:	lea    0x15cb85(%rip),%rdx        # 3a15c8 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x138>
  244a43:	mov    %r12,%rdi
  244a46:	mov    %r12,%rsi
  244a49:	call   *0x163451(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244a4f:	lea    0x15cb8a(%rip),%rdx        # 3a15e0 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x150>
  244a56:	mov    %rdi,%rsi
  244a59:	call   *0x163441(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244a5f:	lea    0x15caea(%rip),%rdx        # 3a1550 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0xc0>
  244a66:	mov    0x28(%rsp),%rdi
  244a6b:	mov    %rdi,%rsi
  244a6e:	call   *0x16342c(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244a74:	lea    0x15caed(%rip),%rdx        # 3a1568 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0xd8>
  244a7b:	mov    0x8(%rsp),%rdi
  244a80:	mov    %rdi,%rsi
  244a83:	call   *0x163417(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244a89:	dec    %rdi
  244a8c:	lea    0x15cb05(%rip),%rdx        # 3a1598 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x108>
  244a93:	mov    0x8(%rsp),%rsi
  244a98:	call   *0x163402(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244a9e:	lea    0x15cb0b(%rip),%rdx        # 3a15b0 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x120>
  244aa5:	mov    %rcx,%rsi
  244aa8:	call   *0x1633f2(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  244aae:	int3
  244aaf:	int3

