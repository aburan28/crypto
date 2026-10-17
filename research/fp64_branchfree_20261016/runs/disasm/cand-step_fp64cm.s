000000000025e7d0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm>:
  25e7d0:	push   %rbp
  25e7d1:	push   %r15
  25e7d3:	push   %r14
  25e7d5:	push   %r13
  25e7d7:	push   %r12
  25e7d9:	push   %rbx
  25e7da:	sub    $0x58,%rsp
  25e7de:	mov    %rcx,0x50(%rsp)
  25e7e3:	mov    %rdi,%r13
  25e7e6:	mov    0x10(%rsi),%r15
  25e7ea:	test   %r15,%r15
  25e7ed:	mov    %rsi,0x8(%rsp)
  25e7f2:	mov    %rdx,0x48(%rsp)
  25e7f7:	je     25e8ea <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x11a>
  25e7fd:	mov    0x40(%rsi),%rdi
  25e801:	mov    0x8(%rsi),%rax
  25e805:	mov    0x38(%rsi),%r9
  25e809:	mov    0x120(%r13),%r10
  25e810:	xor    %ecx,%ecx
  25e812:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  25e820:	cmp    %rcx,%rdi
  25e823:	je     25ebbc <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3ec>
  25e829:	mov    (%rax,%rcx,8),%r11
  25e82d:	mov    %r11d,%r8d
  25e830:	and    $0x3e,%r8d
  25e834:	mov    (%rdx,%r8,4),%r8
  25e838:	sub    %r11,%r8
  25e83b:	mov    $0x0,%r11d
  25e841:	cmovb  %r10,%r11
  25e845:	add    %r11,%r8
  25e848:	mov    %r8,(%r9,%rcx,8)
  25e84c:	inc    %rcx
  25e84f:	cmp    %rcx,%r15
  25e852:	jne    25e820 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x50>
  25e854:	mov    0x100(%r13),%rcx
  25e85b:	mov    0x128(%r13),%r11
  25e862:	mov    0x50(%rsi),%rbx
  25e866:	mov    0x58(%rsi),%r12
  25e86a:	mov    0x60(%rsi),%r14
  25e86e:	xor    %ebp,%ebp
  25e870:	jmp    25e8b5 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xe5>
  25e872:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  25e880:	mul    %rcx
  25e883:	mov    %rax,%r8
  25e886:	mov    %r11,%rax
  25e889:	imul   %r8,%rax
  25e88d:	mov    %rdx,%rcx
  25e890:	mul    %r10
  25e893:	cmp    $0x1,%r8
  25e897:	sbb    $0xffffffffffffffff,%rcx
  25e89b:	add    %rdx,%rcx
  25e89e:	cmp    %r10,%rcx
  25e8a1:	mov    $0x0,%eax
  25e8a6:	cmovae %r10,%rax
  25e8aa:	sub    %rax,%rcx
  25e8ad:	inc    %rbp
  25e8b0:	cmp    %rbp,%r15
  25e8b3:	je     25e8f1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x121>
  25e8b5:	cmp    %rbp,%r12
  25e8b8:	je     25ebcc <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3fc>
  25e8be:	mov    %rcx,(%rbx,%rbp,8)
  25e8c2:	cmp    %rbp,%rdi
  25e8c5:	je     25ebdf <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x40f>
  25e8cb:	mov    (%r9,%rbp,8),%rax
  25e8cf:	test   %rax,%rax
  25e8d2:	jne    25e880 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xb0>
  25e8d4:	inc    %r14
  25e8d7:	mov    0x8(%rsp),%rax
  25e8dc:	mov    %r14,0x60(%rax)
  25e8e0:	inc    %rbp
  25e8e3:	cmp    %rbp,%r15
  25e8e6:	jne    25e8b5 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0xe5>
  25e8e8:	jmp    25e8f1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x121>
  25e8ea:	mov    0x100(%r13),%rcx
  25e8f1:	mov    %r13,%rdi
  25e8f4:	mov    %rcx,%rsi
  25e8f7:	call   *0x16b2f3(%rip)        # 3c9bf0 <_DYNAMIC+0xeb0>
  25e8fd:	test   $0x1,%al
  25e8ff:	je     25eba4 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3d4>
  25e905:	mov    0x8(%rsp),%rax
  25e90a:	mov    0x38(%rax),%rcx
  25e90e:	mov    %rcx,0x18(%rsp)
  25e913:	mov    0x40(%rax),%rcx
  25e917:	mov    0x120(%r13),%rbp
  25e91e:	mov    0x128(%r13),%rax
  25e925:	mov    %rax,0x30(%rsp)
  25e92a:	lea    -0x1(%r15),%rax
  25e92e:	mov    %r15,%rdi
  25e931:	mov    %rcx,0x10(%rsp)
  25e936:	cmp    %rcx,%rax
  25e939:	jae    25ea02 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x232>
  25e93f:	mov    0x8(%rsp),%rax
  25e944:	mov    0x58(%rax),%rcx
  25e948:	mov    0x50(%rax),%rsi
  25e94c:	xor    %r11d,%r11d
  25e94f:	mov    %r15,%rdi
  25e952:	mov    0x30(%rsp),%r14
  25e957:	nopw   0x0(%rax,%rax,1)
  25e960:	test   %rdi,%rdi
  25e963:	je     25ea0b <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x23b>
  25e969:	mov    0x18(%rsp),%rax
  25e96e:	mov    -0x8(%rax,%rdi,8),%r9
  25e973:	dec    %rdi
  25e976:	test   %r9,%r9
  25e979:	je     25e960 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x190>
  25e97b:	cmp    %rcx,%rdi
  25e97e:	jae    25ec2e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x45e>
  25e984:	mov    %rdx,%rax
  25e987:	mov    %rdx,%rbx
  25e98a:	mulq   (%rsi,%rdi,8)
  25e98e:	mov    %rax,%r10
  25e991:	mov    %rdx,%r8
  25e994:	mov    %r14,%rax
  25e997:	imul   %r10,%rax
  25e99b:	mul    %rbp
  25e99e:	cmp    $0x1,%r10
  25e9a2:	sbb    $0xffffffffffffffff,%r8
  25e9a6:	add    %rdx,%r8
  25e9a9:	cmp    %rbp,%r8
  25e9ac:	mov    %rbp,%rax
  25e9af:	cmovb  %r11,%rax
  25e9b3:	sub    %rax,%r8
  25e9b6:	mov    %r9,%rax
  25e9b9:	mul    %rbx
  25e9bc:	mov    %rax,%r9
  25e9bf:	mov    %r14,%rax
  25e9c2:	imul   %r9,%rax
  25e9c6:	mov    %rdx,%r10
  25e9c9:	mul    %rbp
  25e9cc:	mov    %rdx,%rax
  25e9cf:	mov    %r10,%rdx
  25e9d2:	cmp    $0x1,%r9
  25e9d6:	sbb    $0xffffffffffffffff,%rdx
  25e9da:	add    %rax,%rdx
  25e9dd:	cmp    %rbp,%rdx
  25e9e0:	mov    %rbp,%rax
  25e9e3:	cmovb  %r11,%rax
  25e9e7:	mov    0x18(%rsp),%r9
  25e9ec:	mov    %r8,(%r9,%rdi,8)
  25e9f0:	sub    %rax,%rdx
  25e9f3:	lea    -0x1(%rdi),%rax
  25e9f7:	cmp    0x10(%rsp),%rax
  25e9fc:	jb     25e960 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x190>
  25ea02:	test   %rdi,%rdi
  25ea05:	jne    25ec19 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x449>
  25ea0b:	test   %r15,%r15
  25ea0e:	je     25ebad <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3dd>
  25ea14:	mov    0x8(%rsp),%rax
  25ea19:	mov    0x8(%rax),%rcx
  25ea1d:	mov    %rcx,0x40(%rsp)
  25ea22:	mov    0x20(%rax),%rcx
  25ea26:	mov    %rcx,0x38(%rsp)
  25ea2b:	mov    0x28(%rax),%rcx
  25ea2f:	mov    %rcx,0x28(%rsp)
  25ea34:	mov    0x68(%rax),%rax
  25ea38:	mov    %rax,0x20(%rsp)
  25ea3d:	xor    %esi,%esi
  25ea3f:	mov    0x40(%rsp),%rbx
  25ea44:	mov    0x48(%rsp),%r12
  25ea49:	mov    0x38(%rsp),%r13
  25ea4e:	jmp    25ea5c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x28c>
  25ea50:	inc    %rsi
  25ea53:	cmp    %rsi,%r15
  25ea56:	je     25ebad <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x3dd>
  25ea5c:	cmp    %rsi,0x28(%rsp)
  25ea61:	je     25ebef <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x41f>
  25ea67:	cmp    %rsi,0x10(%rsp)
  25ea6c:	je     25ec04 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x434>
  25ea72:	mov    (%rbx,%rsi,8),%rdi
  25ea76:	mov    %edi,%r11d
  25ea79:	and    $0x3e,%r11d
  25ea7d:	mov    0x100(%r12,%r11,4),%rax
  25ea85:	mov    0x0(%r13,%rsi,8),%r14
  25ea8a:	sub    %r14,%rax
  25ea8d:	mov    $0x0,%ecx
  25ea92:	cmovb  %rbp,%rcx
  25ea96:	add    %rcx,%rax
  25ea99:	mov    0x18(%rsp),%rcx
  25ea9e:	mulq   (%rcx,%rsi,8)
  25eaa2:	mov    %rax,%r8
  25eaa5:	mov    %rdx,%rcx
  25eaa8:	mov    0x30(%rsp),%r10
  25eaad:	mov    %r10,%rax
  25eab0:	imul   %r8,%rax
  25eab4:	mul    %rbp
  25eab7:	cmp    $0x1,%r8
  25eabb:	sbb    $0xffffffffffffffff,%rcx
  25eabf:	add    %rdx,%rcx
  25eac2:	cmp    %rbp,%rcx
  25eac5:	mov    %rbp,%rax
  25eac8:	mov    $0x0,%edx
  25eacd:	cmovb  %rdx,%rax
  25ead1:	sub    %rax,%rcx
  25ead4:	mov    %rcx,%rax
  25ead7:	mul    %rcx
  25eada:	mov    %rax,%r8
  25eadd:	mov    %rdx,%r9
  25eae0:	mov    %r10,%rax
  25eae3:	imul   %r8,%rax
  25eae7:	mul    %rbp
  25eaea:	cmp    $0x1,%r8
  25eaee:	sbb    $0xffffffffffffffff,%r9
  25eaf2:	add    %rdx,%r9
  25eaf5:	cmp    %rbp,%r9
  25eaf8:	mov    %rbp,%rax
  25eafb:	mov    $0x0,%edx
  25eb00:	cmovb  %rdx,%rax
  25eb04:	sub    %rax,%r9
  25eb07:	sub    %rdi,%r9
  25eb0a:	mov    $0x0,%eax
  25eb0f:	cmovb  %rbp,%rax
  25eb13:	add    %r9,%rax
  25eb16:	sub    (%r12,%r11,4),%rax
  25eb1a:	mov    $0x0,%edx
  25eb1f:	cmovb  %rbp,%rdx
  25eb23:	lea    (%rdx,%rax,1),%r9
  25eb27:	sub    %r9,%rdi
  25eb2a:	mov    $0x0,%eax
  25eb2f:	cmovb  %rbp,%rax
  25eb33:	add    %rdi,%rax
  25eb36:	mul    %rcx
  25eb39:	mov    %rax,%rcx
  25eb3c:	mov    %rdx,%r8
  25eb3f:	imul   %rax,%r10
  25eb43:	mov    %r10,%rax
  25eb46:	mul    %rbp
  25eb49:	cmp    $0x1,%rcx
  25eb4d:	sbb    $0xffffffffffffffff,%r8
  25eb51:	add    %rdx,%r8
  25eb54:	cmp    %rbp,%r8
  25eb57:	mov    %rbp,%rax
  25eb5a:	mov    $0x0,%ecx
  25eb5f:	cmovb  %rcx,%rax
  25eb63:	sub    %rax,%r8
  25eb66:	sub    %r14,%r8
  25eb69:	mov    $0x0,%eax
  25eb6e:	cmovb  %rbp,%rax
  25eb72:	add    %r8,%rax
  25eb75:	mov    %r9,(%rbx,%rsi,8)
  25eb79:	mov    %rax,0x0(%r13,%rsi,8)
  25eb7e:	test   %r9,0x50(%rsp)
  25eb83:	jne    25ea50 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x280>
  25eb89:	mov    0x20(%rsp),%rcx
  25eb8e:	inc    %rcx
  25eb91:	mov    0x8(%rsp),%rax
  25eb96:	mov    %rcx,0x20(%rsp)
  25eb9b:	mov    %rcx,0x68(%rax)
  25eb9f:	jmp    25ea50 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64cm+0x280>
  25eba4:	mov    0x8(%rsp),%rax
  25eba9:	incq   0x60(%rax)
  25ebad:	add    $0x58,%rsp
  25ebb1:	pop    %rbx
  25ebb2:	pop    %r12
  25ebb4:	pop    %r13
  25ebb6:	pop    %r14
  25ebb8:	pop    %r15
  25ebba:	pop    %rbp
  25ebbb:	ret
  25ebbc:	lea    0x1637dd(%rip),%rdx        # 3c23a0 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0xd8>
  25ebc3:	mov    %rdi,%rsi
  25ebc6:	call   *0x16a664(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ebcc:	lea    0x16379d(%rip),%rdx        # 3c2370 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0xa8>
  25ebd3:	mov    %r12,%rdi
  25ebd6:	mov    %r12,%rsi
  25ebd9:	call   *0x16a651(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ebdf:	lea    0x1637a2(%rip),%rdx        # 3c2388 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0xc0>
  25ebe6:	mov    %rdi,%rsi
  25ebe9:	call   *0x16a641(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ebef:	lea    0x16371a(%rip),%rdx        # 3c2310 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x48>
  25ebf6:	mov    0x28(%rsp),%rdi
  25ebfb:	mov    %rdi,%rsi
  25ebfe:	call   *0x16a62c(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ec04:	lea    0x16371d(%rip),%rdx        # 3c2328 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x60>
  25ec0b:	mov    0x10(%rsp),%rdi
  25ec10:	mov    %rdi,%rsi
  25ec13:	call   *0x16a617(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ec19:	dec    %rdi
  25ec1c:	lea    0x16371d(%rip),%rdx        # 3c2340 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x78>
  25ec23:	mov    0x10(%rsp),%rsi
  25ec28:	call   *0x16a602(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ec2e:	lea    0x163723(%rip),%rdx        # 3c2358 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x90>
  25ec35:	mov    %rcx,%rsi
  25ec38:	call   *0x16a5f2(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25ec3e:	int3
  25ec3f:	int3

