000000000025f390 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64>:
  25f390:	push   %rbp
  25f391:	push   %r15
  25f393:	push   %r14
  25f395:	push   %r13
  25f397:	push   %r12
  25f399:	push   %rbx
  25f39a:	sub    $0x58,%rsp
  25f39e:	mov    %rcx,0x50(%rsp)
  25f3a3:	mov    %rdi,%r13
  25f3a6:	mov    0x10(%rsi),%r15
  25f3aa:	test   %r15,%r15
  25f3ad:	mov    %rsi,0x8(%rsp)
  25f3b2:	mov    %rdx,0x48(%rsp)
  25f3b7:	je     25f4aa <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x11a>
  25f3bd:	mov    0x40(%rsi),%rdi
  25f3c1:	mov    0x8(%rsi),%rax
  25f3c5:	mov    0x38(%rsi),%r9
  25f3c9:	mov    0xf0(%r13),%r10
  25f3d0:	xor    %ecx,%ecx
  25f3d2:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  25f3e0:	cmp    %rcx,%rdi
  25f3e3:	je     25f77c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ec>
  25f3e9:	mov    (%rax,%rcx,8),%r11
  25f3ed:	mov    %r11d,%r8d
  25f3f0:	and    $0x3e,%r8d
  25f3f4:	mov    (%rdx,%r8,4),%r8
  25f3f8:	sub    %r11,%r8
  25f3fb:	mov    $0x0,%r11d
  25f401:	cmovb  %r10,%r11
  25f405:	add    %r11,%r8
  25f408:	mov    %r8,(%r9,%rcx,8)
  25f40c:	inc    %rcx
  25f40f:	cmp    %rcx,%r15
  25f412:	jne    25f3e0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x50>
  25f414:	mov    0xf8(%r13),%r11
  25f41b:	mov    0x100(%r13),%rcx
  25f422:	mov    0x50(%rsi),%rbx
  25f426:	mov    0x58(%rsi),%r12
  25f42a:	mov    0x60(%rsi),%r14
  25f42e:	xor    %ebp,%ebp
  25f430:	jmp    25f475 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xe5>
  25f432:	data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  25f440:	mul    %rcx
  25f443:	mov    %rax,%r8
  25f446:	mov    %r11,%rax
  25f449:	imul   %r8,%rax
  25f44d:	mov    %rdx,%rcx
  25f450:	mul    %r10
  25f453:	cmp    $0x1,%r8
  25f457:	sbb    $0xffffffffffffffff,%rcx
  25f45b:	add    %rdx,%rcx
  25f45e:	cmp    %r10,%rcx
  25f461:	mov    $0x0,%eax
  25f466:	cmovae %r10,%rax
  25f46a:	sub    %rax,%rcx
  25f46d:	inc    %rbp
  25f470:	cmp    %rbp,%r15
  25f473:	je     25f4b1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x121>
  25f475:	cmp    %rbp,%r12
  25f478:	je     25f78c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3fc>
  25f47e:	mov    %rcx,(%rbx,%rbp,8)
  25f482:	cmp    %rbp,%rdi
  25f485:	je     25f79f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x40f>
  25f48b:	mov    (%r9,%rbp,8),%rax
  25f48f:	test   %rax,%rax
  25f492:	jne    25f440 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xb0>
  25f494:	inc    %r14
  25f497:	mov    0x8(%rsp),%rax
  25f49c:	mov    %r14,0x60(%rax)
  25f4a0:	inc    %rbp
  25f4a3:	cmp    %rbp,%r15
  25f4a6:	jne    25f475 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xe5>
  25f4a8:	jmp    25f4b1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x121>
  25f4aa:	mov    0x100(%r13),%rcx
  25f4b1:	mov    %r13,%rdi
  25f4b4:	mov    %rcx,%rsi
  25f4b7:	call   *0x16a733(%rip)        # 3c9bf0 <_DYNAMIC+0xeb0>
  25f4bd:	test   $0x1,%al
  25f4bf:	je     25f764 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3d4>
  25f4c5:	mov    0x8(%rsp),%rax
  25f4ca:	mov    0x38(%rax),%rcx
  25f4ce:	mov    %rcx,0x18(%rsp)
  25f4d3:	mov    0x40(%rax),%rcx
  25f4d7:	mov    0xf0(%r13),%rbp
  25f4de:	mov    0xf8(%r13),%rax
  25f4e5:	mov    %rax,0x30(%rsp)
  25f4ea:	lea    -0x1(%r15),%rax
  25f4ee:	mov    %r15,%rdi
  25f4f1:	mov    %rcx,0x10(%rsp)
  25f4f6:	cmp    %rcx,%rax
  25f4f9:	jae    25f5c2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x232>
  25f4ff:	mov    0x8(%rsp),%rax
  25f504:	mov    0x58(%rax),%rcx
  25f508:	mov    0x50(%rax),%rsi
  25f50c:	xor    %r11d,%r11d
  25f50f:	mov    %r15,%rdi
  25f512:	mov    0x30(%rsp),%r14
  25f517:	nopw   0x0(%rax,%rax,1)
  25f520:	test   %rdi,%rdi
  25f523:	je     25f5cb <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x23b>
  25f529:	mov    0x18(%rsp),%rax
  25f52e:	mov    -0x8(%rax,%rdi,8),%r9
  25f533:	dec    %rdi
  25f536:	test   %r9,%r9
  25f539:	je     25f520 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x190>
  25f53b:	cmp    %rcx,%rdi
  25f53e:	jae    25f7ee <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x45e>
  25f544:	mov    %rdx,%rax
  25f547:	mov    %rdx,%rbx
  25f54a:	mulq   (%rsi,%rdi,8)
  25f54e:	mov    %rax,%r10
  25f551:	mov    %rdx,%r8
  25f554:	mov    %r14,%rax
  25f557:	imul   %r10,%rax
  25f55b:	mul    %rbp
  25f55e:	cmp    $0x1,%r10
  25f562:	sbb    $0xffffffffffffffff,%r8
  25f566:	add    %rdx,%r8
  25f569:	cmp    %rbp,%r8
  25f56c:	mov    %rbp,%rax
  25f56f:	cmovb  %r11,%rax
  25f573:	sub    %rax,%r8
  25f576:	mov    %r9,%rax
  25f579:	mul    %rbx
  25f57c:	mov    %rax,%r9
  25f57f:	mov    %r14,%rax
  25f582:	imul   %r9,%rax
  25f586:	mov    %rdx,%r10
  25f589:	mul    %rbp
  25f58c:	mov    %rdx,%rax
  25f58f:	mov    %r10,%rdx
  25f592:	cmp    $0x1,%r9
  25f596:	sbb    $0xffffffffffffffff,%rdx
  25f59a:	add    %rax,%rdx
  25f59d:	cmp    %rbp,%rdx
  25f5a0:	mov    %rbp,%rax
  25f5a3:	cmovb  %r11,%rax
  25f5a7:	mov    0x18(%rsp),%r9
  25f5ac:	mov    %r8,(%r9,%rdi,8)
  25f5b0:	sub    %rax,%rdx
  25f5b3:	lea    -0x1(%rdi),%rax
  25f5b7:	cmp    0x10(%rsp),%rax
  25f5bc:	jb     25f520 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x190>
  25f5c2:	test   %rdi,%rdi
  25f5c5:	jne    25f7d9 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x449>
  25f5cb:	test   %r15,%r15
  25f5ce:	je     25f76d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3dd>
  25f5d4:	mov    0x8(%rsp),%rax
  25f5d9:	mov    0x8(%rax),%rcx
  25f5dd:	mov    %rcx,0x40(%rsp)
  25f5e2:	mov    0x20(%rax),%rcx
  25f5e6:	mov    %rcx,0x38(%rsp)
  25f5eb:	mov    0x28(%rax),%rcx
  25f5ef:	mov    %rcx,0x28(%rsp)
  25f5f4:	mov    0x68(%rax),%rax
  25f5f8:	mov    %rax,0x20(%rsp)
  25f5fd:	xor    %esi,%esi
  25f5ff:	mov    0x40(%rsp),%rbx
  25f604:	mov    0x48(%rsp),%r12
  25f609:	mov    0x38(%rsp),%r13
  25f60e:	jmp    25f61c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x28c>
  25f610:	inc    %rsi
  25f613:	cmp    %rsi,%r15
  25f616:	je     25f76d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3dd>
  25f61c:	cmp    %rsi,0x28(%rsp)
  25f621:	je     25f7af <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x41f>
  25f627:	cmp    %rsi,0x10(%rsp)
  25f62c:	je     25f7c4 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x434>
  25f632:	mov    (%rbx,%rsi,8),%rdi
  25f636:	mov    %edi,%r11d
  25f639:	and    $0x3e,%r11d
  25f63d:	mov    0x100(%r12,%r11,4),%rax
  25f645:	mov    0x0(%r13,%rsi,8),%r14
  25f64a:	sub    %r14,%rax
  25f64d:	mov    $0x0,%ecx
  25f652:	cmovb  %rbp,%rcx
  25f656:	add    %rcx,%rax
  25f659:	mov    0x18(%rsp),%rcx
  25f65e:	mulq   (%rcx,%rsi,8)
  25f662:	mov    %rax,%r8
  25f665:	mov    %rdx,%rcx
  25f668:	mov    0x30(%rsp),%r10
  25f66d:	mov    %r10,%rax
  25f670:	imul   %r8,%rax
  25f674:	mul    %rbp
  25f677:	cmp    $0x1,%r8
  25f67b:	sbb    $0xffffffffffffffff,%rcx
  25f67f:	add    %rdx,%rcx
  25f682:	cmp    %rbp,%rcx
  25f685:	mov    %rbp,%rax
  25f688:	mov    $0x0,%edx
  25f68d:	cmovb  %rdx,%rax
  25f691:	sub    %rax,%rcx
  25f694:	mov    %rcx,%rax
  25f697:	mul    %rcx
  25f69a:	mov    %rax,%r8
  25f69d:	mov    %rdx,%r9
  25f6a0:	mov    %r10,%rax
  25f6a3:	imul   %r8,%rax
  25f6a7:	mul    %rbp
  25f6aa:	cmp    $0x1,%r8
  25f6ae:	sbb    $0xffffffffffffffff,%r9
  25f6b2:	add    %rdx,%r9
  25f6b5:	cmp    %rbp,%r9
  25f6b8:	mov    %rbp,%rax
  25f6bb:	mov    $0x0,%edx
  25f6c0:	cmovb  %rdx,%rax
  25f6c4:	sub    %rax,%r9
  25f6c7:	sub    %rdi,%r9
  25f6ca:	mov    $0x0,%eax
  25f6cf:	cmovb  %rbp,%rax
  25f6d3:	add    %r9,%rax
  25f6d6:	sub    (%r12,%r11,4),%rax
  25f6da:	mov    $0x0,%edx
  25f6df:	cmovb  %rbp,%rdx
  25f6e3:	lea    (%rdx,%rax,1),%r9
  25f6e7:	sub    %r9,%rdi
  25f6ea:	mov    $0x0,%eax
  25f6ef:	cmovb  %rbp,%rax
  25f6f3:	add    %rdi,%rax
  25f6f6:	mul    %rcx
  25f6f9:	mov    %rax,%rcx
  25f6fc:	mov    %rdx,%r8
  25f6ff:	imul   %rax,%r10
  25f703:	mov    %r10,%rax
  25f706:	mul    %rbp
  25f709:	cmp    $0x1,%rcx
  25f70d:	sbb    $0xffffffffffffffff,%r8
  25f711:	add    %rdx,%r8
  25f714:	cmp    %rbp,%r8
  25f717:	mov    %rbp,%rax
  25f71a:	mov    $0x0,%ecx
  25f71f:	cmovb  %rcx,%rax
  25f723:	sub    %rax,%r8
  25f726:	sub    %r14,%r8
  25f729:	mov    $0x0,%eax
  25f72e:	cmovb  %rbp,%rax
  25f732:	add    %r8,%rax
  25f735:	mov    %r9,(%rbx,%rsi,8)
  25f739:	mov    %rax,0x0(%r13,%rsi,8)
  25f73e:	test   %r9,0x50(%rsp)
  25f743:	jne    25f610 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x280>
  25f749:	mov    0x20(%rsp),%rcx
  25f74e:	inc    %rcx
  25f751:	mov    0x8(%rsp),%rax
  25f756:	mov    %rcx,0x20(%rsp)
  25f75b:	mov    %rcx,0x68(%rax)
  25f75f:	jmp    25f610 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x280>
  25f764:	mov    0x8(%rsp),%rax
  25f769:	incq   0x60(%rax)
  25f76d:	add    $0x58,%rsp
  25f771:	pop    %rbx
  25f772:	pop    %r12
  25f774:	pop    %r13
  25f776:	pop    %r14
  25f778:	pop    %r15
  25f77a:	pop    %rbp
  25f77b:	ret
  25f77c:	lea    0x162c1d(%rip),%rdx        # 3c23a0 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0xd8>
  25f783:	mov    %rdi,%rsi
  25f786:	call   *0x169aa4(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f78c:	lea    0x162bdd(%rip),%rdx        # 3c2370 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0xa8>
  25f793:	mov    %r12,%rdi
  25f796:	mov    %r12,%rsi
  25f799:	call   *0x169a91(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f79f:	lea    0x162be2(%rip),%rdx        # 3c2388 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0xc0>
  25f7a6:	mov    %rdi,%rsi
  25f7a9:	call   *0x169a81(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f7af:	lea    0x162b5a(%rip),%rdx        # 3c2310 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x48>
  25f7b6:	mov    0x28(%rsp),%rdi
  25f7bb:	mov    %rdi,%rsi
  25f7be:	call   *0x169a6c(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f7c4:	lea    0x162b5d(%rip),%rdx        # 3c2328 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x60>
  25f7cb:	mov    0x10(%rsp),%rdi
  25f7d0:	mov    %rdi,%rsi
  25f7d3:	call   *0x169a57(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f7d9:	dec    %rdi
  25f7dc:	lea    0x162b5d(%rip),%rdx        # 3c2340 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x78>
  25f7e3:	mov    0x10(%rsp),%rsi
  25f7e8:	call   *0x169a42(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f7ee:	lea    0x162b63(%rip),%rdx        # 3c2358 <anon.3c11ed13effa8077175ba73b06d1f882.64.llvm.7896473603554930564+0x90>
  25f7f5:	mov    %rcx,%rsi
  25f7f8:	call   *0x169a32(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  25f7fe:	int3
  25f7ff:	int3

