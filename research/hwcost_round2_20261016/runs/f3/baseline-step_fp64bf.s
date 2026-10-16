0000000000243670 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf>:
  243670:	push   %rbp
  243671:	push   %r15
  243673:	push   %r14
  243675:	push   %r13
  243677:	push   %r12
  243679:	push   %rbx
  24367a:	sub    $0x38,%rsp
  24367e:	mov    %rcx,0x30(%rsp)
  243683:	mov    %rdx,%r15
  243686:	mov    %rdi,%r11
  243689:	mov    0x10(%rsi),%r13
  24368d:	test   %r13,%r13
  243690:	mov    %rsi,0x10(%rsp)
  243695:	mov    %rdi,0x8(%rsp)
  24369a:	je     24378a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x11a>
  2436a0:	mov    0x40(%rsi),%rdi
  2436a4:	mov    0x8(%rsi),%rax
  2436a8:	mov    0x38(%rsi),%r9
  2436ac:	mov    0x120(%r11),%r10
  2436b3:	xor    %ecx,%ecx
  2436b5:	data16 cs nopw 0x0(%rax,%rax,1)
  2436c0:	cmp    %rcx,%rdi
  2436c3:	je     243a6c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3fc>
  2436c9:	mov    (%rax,%rcx,8),%rdx
  2436cd:	mov    %edx,%r8d
  2436d0:	and    $0x3e,%r8d
  2436d4:	mov    (%r15,%r8,4),%r8
  2436d8:	sub    %rdx,%r8
  2436db:	mov    $0x0,%edx
  2436e0:	cmovb  %r10,%rdx
  2436e4:	add    %r8,%rdx
  2436e7:	mov    %rdx,(%r9,%rcx,8)
  2436eb:	inc    %rcx
  2436ee:	cmp    %rcx,%r13
  2436f1:	jne    2436c0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x50>
  2436f3:	mov    0x100(%r11),%rcx
  2436fa:	mov    0x128(%r11),%r11
  243701:	mov    0x50(%rsi),%rbx
  243705:	mov    0x58(%rsi),%r12
  243709:	mov    0x60(%rsi),%r14
  24370d:	xor    %ebp,%ebp
  24370f:	jmp    243755 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0xe5>
  243711:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  243720:	mul    %rcx
  243723:	mov    %rax,%r8
  243726:	mov    %r11,%rax
  243729:	imul   %r8,%rax
  24372d:	mov    %rdx,%rcx
  243730:	mul    %r10
  243733:	cmp    $0x1,%r8
  243737:	sbb    $0xffffffffffffffff,%rcx
  24373b:	add    %rdx,%rcx
  24373e:	cmp    %r10,%rcx
  243741:	mov    $0x0,%eax
  243746:	cmovae %r10,%rax
  24374a:	sub    %rax,%rcx
  24374d:	inc    %rbp
  243750:	cmp    %rbp,%r13
  243753:	je     243791 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x121>
  243755:	cmp    %rbp,%r12
  243758:	je     243a7c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x40c>
  24375e:	mov    %rcx,(%rbx,%rbp,8)
  243762:	cmp    %rbp,%rdi
  243765:	je     243a8f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x41f>
  24376b:	mov    (%r9,%rbp,8),%rax
  24376f:	test   %rax,%rax
  243772:	jne    243720 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0xb0>
  243774:	inc    %r14
  243777:	mov    0x10(%rsp),%rax
  24377c:	mov    %r14,0x60(%rax)
  243780:	inc    %rbp
  243783:	cmp    %rbp,%r13
  243786:	jne    243755 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0xe5>
  243788:	jmp    243791 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x121>
  24378a:	mov    0x100(%r11),%rcx
  243791:	mov    0x8(%rsp),%rbx
  243796:	mov    %rbx,%rdi
  243799:	mov    %rcx,%rsi
  24379c:	call   *0x1651a6(%rip)        # 3a8948 <_DYNAMIC+0xf90>
  2437a2:	test   $0x1,%al
  2437a4:	je     243a54 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3e4>
  2437aa:	mov    0x10(%rsp),%rsi
  2437af:	mov    0x38(%rsi),%rax
  2437b3:	mov    %rax,0x18(%rsp)
  2437b8:	mov    0x40(%rsi),%rcx
  2437bc:	mov    0x120(%rbx),%rbp
  2437c3:	mov    0x128(%rbx),%r12
  2437ca:	lea    -0x1(%r13),%rax
  2437ce:	mov    %r13,%rdi
  2437d1:	mov    %rcx,0x8(%rsp)
  2437d6:	cmp    %rcx,%rax
  2437d9:	jae    243892 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x222>
  2437df:	mov    0x58(%rsi),%rcx
  2437e3:	mov    0x50(%rsi),%rsi
  2437e7:	xor    %r11d,%r11d
  2437ea:	mov    %r13,%rdi
  2437ed:	nopl   (%rax)
  2437f0:	test   %rdi,%rdi
  2437f3:	je     24389b <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x22b>
  2437f9:	mov    0x18(%rsp),%rax
  2437fe:	mov    -0x8(%rax,%rdi,8),%r9
  243803:	dec    %rdi
  243806:	test   %r9,%r9
  243809:	je     2437f0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x180>
  24380b:	cmp    %rcx,%rdi
  24380e:	jae    243ade <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x46e>
  243814:	mov    %rdx,%rax
  243817:	mov    %rdx,%rbx
  24381a:	mulq   (%rsi,%rdi,8)
  24381e:	mov    %rax,%r10
  243821:	mov    %rdx,%r8
  243824:	mov    %r12,%rax
  243827:	imul   %r10,%rax
  24382b:	mul    %rbp
  24382e:	cmp    $0x1,%r10
  243832:	sbb    $0xffffffffffffffff,%r8
  243836:	add    %rdx,%r8
  243839:	cmp    %rbp,%r8
  24383c:	mov    %rbp,%rax
  24383f:	cmovb  %r11,%rax
  243843:	sub    %rax,%r8
  243846:	mov    %r9,%rax
  243849:	mul    %rbx
  24384c:	mov    %rax,%r9
  24384f:	mov    %r12,%rax
  243852:	imul   %r9,%rax
  243856:	mov    %rdx,%r10
  243859:	mul    %rbp
  24385c:	mov    %rdx,%rax
  24385f:	mov    %r10,%rdx
  243862:	cmp    $0x1,%r9
  243866:	sbb    $0xffffffffffffffff,%rdx
  24386a:	add    %rax,%rdx
  24386d:	cmp    %rbp,%rdx
  243870:	mov    %rbp,%rax
  243873:	cmovb  %r11,%rax
  243877:	mov    0x18(%rsp),%r9
  24387c:	mov    %r8,(%r9,%rdi,8)
  243880:	sub    %rax,%rdx
  243883:	lea    -0x1(%rdi),%rax
  243887:	cmp    0x8(%rsp),%rax
  24388c:	jb     2437f0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x180>
  243892:	test   %rdi,%rdi
  243895:	jne    243ac9 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x459>
  24389b:	test   %r13,%r13
  24389e:	je     243a5d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3ed>
  2438a4:	mov    0x10(%rsp),%rax
  2438a9:	mov    0x8(%rax),%r10
  2438ad:	mov    0x20(%rax),%rbx
  2438b1:	mov    0x28(%rax),%rcx
  2438b5:	mov    %rcx,0x28(%rsp)
  2438ba:	mov    0x68(%rax),%rax
  2438be:	mov    %rax,0x20(%rsp)
  2438c3:	xor    %r14d,%r14d
  2438c6:	jmp    2438f2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x282>
  2438c8:	nopl   0x0(%rax,%rax,1)
  2438d0:	add    %rcx,%rax
  2438d3:	mov    %r9,(%r10,%r14,8)
  2438d7:	mov    %rax,(%rbx,%r14,8)
  2438db:	test   %r9,0x30(%rsp)
  2438e0:	je     243a30 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3c0>
  2438e6:	inc    %r14
  2438e9:	cmp    %r14,%r13
  2438ec:	je     243a5d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3ed>
  2438f2:	cmp    %r14,0x28(%rsp)
  2438f7:	je     243a9f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x42f>
  2438fd:	cmp    %r14,0x8(%rsp)
  243902:	je     243ab4 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x444>
  243908:	mov    (%r10,%r14,8),%rdi
  24390c:	mov    %edi,%r11d
  24390f:	and    $0x3e,%r11d
  243913:	mov    0x100(%r15,%r11,4),%rax
  24391b:	mov    (%rbx,%r14,8),%rsi
  24391f:	mov    %rbp,%rcx
  243922:	sub    %rsi,%rax
  243925:	jb     243929 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x2b9>
  243927:	xor    %ecx,%ecx
  243929:	add    %rcx,%rax
  24392c:	mov    0x18(%rsp),%rcx
  243931:	mulq   (%rcx,%r14,8)
  243935:	mov    %rax,%r8
  243938:	mov    %rdx,%rcx
  24393b:	mov    %r12,%rax
  24393e:	imul   %r8,%rax
  243942:	mul    %rbp
  243945:	cmp    $0x1,%r8
  243949:	sbb    $0xffffffffffffffff,%rcx
  24394d:	add    %rdx,%rcx
  243950:	mov    $0x0,%eax
  243955:	cmp    %rbp,%rcx
  243958:	jb     24395d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x2ed>
  24395a:	mov    %rbp,%rax
  24395d:	sub    %rax,%rcx
  243960:	mov    %rcx,%rax
  243963:	mul    %rcx
  243966:	mov    %rax,%r9
  243969:	mov    %rdx,%r8
  24396c:	mov    %r12,%rax
  24396f:	imul   %r9,%rax
  243973:	mul    %rbp
  243976:	cmp    $0x1,%r9
  24397a:	sbb    $0xffffffffffffffff,%r8
  24397e:	add    %rdx,%r8
  243981:	mov    $0x0,%eax
  243986:	cmp    %rbp,%r8
  243989:	jae    243a00 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x390>
  24398b:	sub    %rax,%r8
  24398e:	mov    %rbp,%rax
  243991:	sub    %rdi,%r8
  243994:	jae    243a0e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x39e>
  243996:	add    %r8,%rax
  243999:	mov    %rbp,%rdx
  24399c:	sub    (%r15,%r11,4),%rax
  2439a0:	jae    243a1c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3ac>
  2439a2:	lea    (%rdx,%rax,1),%r9
  2439a6:	mov    %rbp,%rax
  2439a9:	sub    %r9,%rdi
  2439ac:	jb     2439b0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x340>
  2439ae:	xor    %eax,%eax
  2439b0:	add    %rdi,%rax
  2439b3:	mul    %rcx
  2439b6:	mov    %rax,%r8
  2439b9:	mov    %rdx,%rcx
  2439bc:	mov    %r12,%rax
  2439bf:	imul   %r8,%rax
  2439c3:	mul    %rbp
  2439c6:	cmp    $0x1,%r8
  2439ca:	sbb    $0xffffffffffffffff,%rcx
  2439ce:	add    %rdx,%rcx
  2439d1:	mov    $0x0,%eax
  2439d6:	cmp    %rbp,%rcx
  2439d9:	jb     2439de <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x36e>
  2439db:	mov    %rbp,%rax
  2439de:	sub    %rax,%rcx
  2439e1:	mov    %rbp,%rax
  2439e4:	sub    %rsi,%rcx
  2439e7:	jb     2438d0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x260>
  2439ed:	xor    %eax,%eax
  2439ef:	jmp    2438d0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x260>
  2439f4:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  243a00:	mov    %rbp,%rax
  243a03:	sub    %rax,%r8
  243a06:	mov    %rbp,%rax
  243a09:	sub    %rdi,%r8
  243a0c:	jb     243996 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x326>
  243a0e:	xor    %eax,%eax
  243a10:	add    %r8,%rax
  243a13:	mov    %rbp,%rdx
  243a16:	sub    (%r15,%r11,4),%rax
  243a1a:	jb     2439a2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x332>
  243a1c:	xor    %edx,%edx
  243a1e:	lea    (%rdx,%rax,1),%r9
  243a22:	mov    %rbp,%rax
  243a25:	sub    %r9,%rdi
  243a28:	jae    2439ae <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x33e>
  243a2a:	jmp    2439b0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x340>
  243a2c:	nopl   0x0(%rax)
  243a30:	mov    0x20(%rsp),%rcx
  243a35:	inc    %rcx
  243a38:	mov    0x10(%rsp),%rax
  243a3d:	mov    %rcx,0x20(%rsp)
  243a42:	mov    %rcx,0x68(%rax)
  243a46:	inc    %r14
  243a49:	cmp    %r14,%r13
  243a4c:	jne    2438f2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x282>
  243a52:	jmp    243a5d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x3ed>
  243a54:	mov    0x10(%rsp),%rax
  243a59:	incq   0x60(%rax)
  243a5d:	add    $0x38,%rsp
  243a61:	pop    %rbx
  243a62:	pop    %r12
  243a64:	pop    %r13
  243a66:	pop    %r14
  243a68:	pop    %r15
  243a6a:	pop    %rbp
  243a6b:	ret
  243a6c:	lea    0x15db85(%rip),%rdx        # 3a15f8 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x168>
  243a73:	mov    %rdi,%rsi
  243a76:	call   *0x164424(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243a7c:	lea    0x15db45(%rip),%rdx        # 3a15c8 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x138>
  243a83:	mov    %r12,%rdi
  243a86:	mov    %r12,%rsi
  243a89:	call   *0x164411(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243a8f:	lea    0x15db4a(%rip),%rdx        # 3a15e0 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x150>
  243a96:	mov    %rdi,%rsi
  243a99:	call   *0x164401(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243a9f:	lea    0x15daaa(%rip),%rdx        # 3a1550 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0xc0>
  243aa6:	mov    0x28(%rsp),%rdi
  243aab:	mov    %rdi,%rsi
  243aae:	call   *0x1643ec(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243ab4:	lea    0x15daad(%rip),%rdx        # 3a1568 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0xd8>
  243abb:	mov    0x8(%rsp),%rdi
  243ac0:	mov    %rdi,%rsi
  243ac3:	call   *0x1643d7(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243ac9:	dec    %rdi
  243acc:	lea    0x15dac5(%rip),%rdx        # 3a1598 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x108>
  243ad3:	mov    0x8(%rsp),%rsi
  243ad8:	call   *0x1643c2(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243ade:	lea    0x15dacb(%rip),%rdx        # 3a15b0 <anon.9b47b1dce56bb38c8f4f0638e0cf3e2f.376.llvm.2357065098697664748+0x120>
  243ae5:	mov    %rcx,%rsi
  243ae8:	call   *0x1643b2(%rip)        # 3a7ea0 <_DYNAMIC+0x4e8>
  243aee:	int3
  243aef:	int3

