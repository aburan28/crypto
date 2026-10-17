00000000002c4320 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64>:
  2c4320:	push   %rbp
  2c4321:	push   %r15
  2c4323:	push   %r14
  2c4325:	push   %r13
  2c4327:	push   %r12
  2c4329:	push   %rbx
  2c432a:	sub    $0x38,%rsp
  2c432e:	mov    %rcx,0x30(%rsp)
  2c4333:	mov    %rdx,%r15
  2c4336:	mov    %rdi,%rbx
  2c4339:	mov    0x10(%rsi),%r13
  2c433d:	test   %r13,%r13
  2c4340:	mov    %rsi,0x10(%rsp)
  2c4345:	mov    %rdi,0x8(%rsp)
  2c434a:	je     2c443a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x11a>
  2c4350:	mov    0x40(%rsi),%rdi
  2c4354:	mov    0x8(%rsi),%rax
  2c4358:	mov    0x38(%rsi),%r9
  2c435c:	mov    0xf0(%rbx),%r10
  2c4363:	xor    %ecx,%ecx
  2c4365:	data16 cs nopw 0x0(%rax,%rax,1)
  2c4370:	cmp    %rcx,%rdi
  2c4373:	je     2c471c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3fc>
  2c4379:	mov    (%rax,%rcx,8),%rdx
  2c437d:	mov    %edx,%r8d
  2c4380:	and    $0x3e,%r8d
  2c4384:	mov    (%r15,%r8,4),%r8
  2c4388:	sub    %rdx,%r8
  2c438b:	mov    $0x0,%edx
  2c4390:	cmovb  %r10,%rdx
  2c4394:	add    %r8,%rdx
  2c4397:	mov    %rdx,(%r9,%rcx,8)
  2c439b:	inc    %rcx
  2c439e:	cmp    %rcx,%r13
  2c43a1:	jne    2c4370 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x50>
  2c43a3:	mov    0xf8(%rbx),%r11
  2c43aa:	mov    0x100(%rbx),%rcx
  2c43b1:	mov    0x50(%rsi),%rbx
  2c43b5:	mov    0x58(%rsi),%r12
  2c43b9:	mov    0x60(%rsi),%r14
  2c43bd:	xor    %ebp,%ebp
  2c43bf:	jmp    2c4405 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xe5>
  2c43c1:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  2c43d0:	mul    %rcx
  2c43d3:	mov    %rax,%r8
  2c43d6:	mov    %r11,%rax
  2c43d9:	imul   %r8,%rax
  2c43dd:	mov    %rdx,%rcx
  2c43e0:	mul    %r10
  2c43e3:	cmp    $0x1,%r8
  2c43e7:	sbb    $0xffffffffffffffff,%rcx
  2c43eb:	add    %rdx,%rcx
  2c43ee:	cmp    %r10,%rcx
  2c43f1:	mov    $0x0,%eax
  2c43f6:	cmovae %r10,%rax
  2c43fa:	sub    %rax,%rcx
  2c43fd:	inc    %rbp
  2c4400:	cmp    %rbp,%r13
  2c4403:	je     2c4441 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x121>
  2c4405:	cmp    %rbp,%r12
  2c4408:	je     2c472c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x40c>
  2c440e:	mov    %rcx,(%rbx,%rbp,8)
  2c4412:	cmp    %rbp,%rdi
  2c4415:	je     2c473f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x41f>
  2c441b:	mov    (%r9,%rbp,8),%rax
  2c441f:	test   %rax,%rax
  2c4422:	jne    2c43d0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xb0>
  2c4424:	inc    %r14
  2c4427:	mov    0x10(%rsp),%rax
  2c442c:	mov    %r14,0x60(%rax)
  2c4430:	inc    %rbp
  2c4433:	cmp    %rbp,%r13
  2c4436:	jne    2c4405 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xe5>
  2c4438:	jmp    2c4441 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x121>
  2c443a:	mov    0x100(%rbx),%rcx
  2c4441:	mov    0x8(%rsp),%rbx
  2c4446:	mov    %rbx,%rdi
  2c4449:	mov    %rcx,%rsi
  2c444c:	call   *0x105d1e(%rip)        # 3ca170 <_DYNAMIC+0xf40>
  2c4452:	test   $0x1,%al
  2c4454:	je     2c4704 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3e4>
  2c445a:	mov    0x10(%rsp),%rsi
  2c445f:	mov    0x38(%rsi),%rax
  2c4463:	mov    %rax,0x18(%rsp)
  2c4468:	mov    0x40(%rsi),%rcx
  2c446c:	mov    0xf0(%rbx),%rbp
  2c4473:	mov    0xf8(%rbx),%r12
  2c447a:	lea    -0x1(%r13),%rax
  2c447e:	mov    %r13,%rdi
  2c4481:	mov    %rcx,0x8(%rsp)
  2c4486:	cmp    %rcx,%rax
  2c4489:	jae    2c4542 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x222>
  2c448f:	mov    0x58(%rsi),%rcx
  2c4493:	mov    0x50(%rsi),%rsi
  2c4497:	xor    %r11d,%r11d
  2c449a:	mov    %r13,%rdi
  2c449d:	nopl   (%rax)
  2c44a0:	test   %rdi,%rdi
  2c44a3:	je     2c454b <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x22b>
  2c44a9:	mov    0x18(%rsp),%rax
  2c44ae:	mov    -0x8(%rax,%rdi,8),%r9
  2c44b3:	dec    %rdi
  2c44b6:	test   %r9,%r9
  2c44b9:	je     2c44a0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x180>
  2c44bb:	cmp    %rcx,%rdi
  2c44be:	jae    2c478e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x46e>
  2c44c4:	mov    %rdx,%rax
  2c44c7:	mov    %rdx,%rbx
  2c44ca:	mulq   (%rsi,%rdi,8)
  2c44ce:	mov    %rax,%r10
  2c44d1:	mov    %rdx,%r8
  2c44d4:	mov    %r12,%rax
  2c44d7:	imul   %r10,%rax
  2c44db:	mul    %rbp
  2c44de:	cmp    $0x1,%r10
  2c44e2:	sbb    $0xffffffffffffffff,%r8
  2c44e6:	add    %rdx,%r8
  2c44e9:	cmp    %rbp,%r8
  2c44ec:	mov    %rbp,%rax
  2c44ef:	cmovb  %r11,%rax
  2c44f3:	sub    %rax,%r8
  2c44f6:	mov    %r9,%rax
  2c44f9:	mul    %rbx
  2c44fc:	mov    %rax,%r9
  2c44ff:	mov    %r12,%rax
  2c4502:	imul   %r9,%rax
  2c4506:	mov    %rdx,%r10
  2c4509:	mul    %rbp
  2c450c:	mov    %rdx,%rax
  2c450f:	mov    %r10,%rdx
  2c4512:	cmp    $0x1,%r9
  2c4516:	sbb    $0xffffffffffffffff,%rdx
  2c451a:	add    %rax,%rdx
  2c451d:	cmp    %rbp,%rdx
  2c4520:	mov    %rbp,%rax
  2c4523:	cmovb  %r11,%rax
  2c4527:	mov    0x18(%rsp),%r9
  2c452c:	mov    %r8,(%r9,%rdi,8)
  2c4530:	sub    %rax,%rdx
  2c4533:	lea    -0x1(%rdi),%rax
  2c4537:	cmp    0x8(%rsp),%rax
  2c453c:	jb     2c44a0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x180>
  2c4542:	test   %rdi,%rdi
  2c4545:	jne    2c4779 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x459>
  2c454b:	test   %r13,%r13
  2c454e:	je     2c470d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ed>
  2c4554:	mov    0x10(%rsp),%rax
  2c4559:	mov    0x8(%rax),%r10
  2c455d:	mov    0x20(%rax),%rbx
  2c4561:	mov    0x28(%rax),%rcx
  2c4565:	mov    %rcx,0x28(%rsp)
  2c456a:	mov    0x68(%rax),%rax
  2c456e:	mov    %rax,0x20(%rsp)
  2c4573:	xor    %r14d,%r14d
  2c4576:	jmp    2c45a2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x282>
  2c4578:	nopl   0x0(%rax,%rax,1)
  2c4580:	add    %rcx,%rax
  2c4583:	mov    %r9,(%r10,%r14,8)
  2c4587:	mov    %rax,(%rbx,%r14,8)
  2c458b:	test   %r9,0x30(%rsp)
  2c4590:	je     2c46e0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3c0>
  2c4596:	inc    %r14
  2c4599:	cmp    %r14,%r13
  2c459c:	je     2c470d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ed>
  2c45a2:	cmp    %r14,0x28(%rsp)
  2c45a7:	je     2c474f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x42f>
  2c45ad:	cmp    %r14,0x8(%rsp)
  2c45b2:	je     2c4764 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x444>
  2c45b8:	mov    (%r10,%r14,8),%rdi
  2c45bc:	mov    %edi,%r11d
  2c45bf:	and    $0x3e,%r11d
  2c45c3:	mov    0x100(%r15,%r11,4),%rax
  2c45cb:	mov    (%rbx,%r14,8),%rsi
  2c45cf:	mov    %rbp,%rcx
  2c45d2:	sub    %rsi,%rax
  2c45d5:	jb     2c45d9 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x2b9>
  2c45d7:	xor    %ecx,%ecx
  2c45d9:	add    %rcx,%rax
  2c45dc:	mov    0x18(%rsp),%rcx
  2c45e1:	mulq   (%rcx,%r14,8)
  2c45e5:	mov    %rax,%r8
  2c45e8:	mov    %rdx,%rcx
  2c45eb:	mov    %r12,%rax
  2c45ee:	imul   %r8,%rax
  2c45f2:	mul    %rbp
  2c45f5:	cmp    $0x1,%r8
  2c45f9:	sbb    $0xffffffffffffffff,%rcx
  2c45fd:	add    %rdx,%rcx
  2c4600:	mov    $0x0,%eax
  2c4605:	cmp    %rbp,%rcx
  2c4608:	jb     2c460d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x2ed>
  2c460a:	mov    %rbp,%rax
  2c460d:	sub    %rax,%rcx
  2c4610:	mov    %rcx,%rax
  2c4613:	mul    %rcx
  2c4616:	mov    %rax,%r9
  2c4619:	mov    %rdx,%r8
  2c461c:	mov    %r12,%rax
  2c461f:	imul   %r9,%rax
  2c4623:	mul    %rbp
  2c4626:	cmp    $0x1,%r9
  2c462a:	sbb    $0xffffffffffffffff,%r8
  2c462e:	add    %rdx,%r8
  2c4631:	mov    $0x0,%eax
  2c4636:	cmp    %rbp,%r8
  2c4639:	jae    2c46b0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x390>
  2c463b:	sub    %rax,%r8
  2c463e:	mov    %rbp,%rax
  2c4641:	sub    %rdi,%r8
  2c4644:	jae    2c46be <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x39e>
  2c4646:	add    %r8,%rax
  2c4649:	mov    %rbp,%rdx
  2c464c:	sub    (%r15,%r11,4),%rax
  2c4650:	jae    2c46cc <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ac>
  2c4652:	lea    (%rdx,%rax,1),%r9
  2c4656:	mov    %rbp,%rax
  2c4659:	sub    %r9,%rdi
  2c465c:	jb     2c4660 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x340>
  2c465e:	xor    %eax,%eax
  2c4660:	add    %rdi,%rax
  2c4663:	mul    %rcx
  2c4666:	mov    %rax,%r8
  2c4669:	mov    %rdx,%rcx
  2c466c:	mov    %r12,%rax
  2c466f:	imul   %r8,%rax
  2c4673:	mul    %rbp
  2c4676:	cmp    $0x1,%r8
  2c467a:	sbb    $0xffffffffffffffff,%rcx
  2c467e:	add    %rdx,%rcx
  2c4681:	mov    $0x0,%eax
  2c4686:	cmp    %rbp,%rcx
  2c4689:	jb     2c468e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x36e>
  2c468b:	mov    %rbp,%rax
  2c468e:	sub    %rax,%rcx
  2c4691:	mov    %rbp,%rax
  2c4694:	sub    %rsi,%rcx
  2c4697:	jb     2c4580 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x260>
  2c469d:	xor    %eax,%eax
  2c469f:	jmp    2c4580 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x260>
  2c46a4:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  2c46b0:	mov    %rbp,%rax
  2c46b3:	sub    %rax,%r8
  2c46b6:	mov    %rbp,%rax
  2c46b9:	sub    %rdi,%r8
  2c46bc:	jb     2c4646 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x326>
  2c46be:	xor    %eax,%eax
  2c46c0:	add    %r8,%rax
  2c46c3:	mov    %rbp,%rdx
  2c46c6:	sub    (%r15,%r11,4),%rax
  2c46ca:	jb     2c4652 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x332>
  2c46cc:	xor    %edx,%edx
  2c46ce:	lea    (%rdx,%rax,1),%r9
  2c46d2:	mov    %rbp,%rax
  2c46d5:	sub    %r9,%rdi
  2c46d8:	jae    2c465e <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x33e>
  2c46da:	jmp    2c4660 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x340>
  2c46dc:	nopl   0x0(%rax)
  2c46e0:	mov    0x20(%rsp),%rcx
  2c46e5:	inc    %rcx
  2c46e8:	mov    0x10(%rsp),%rax
  2c46ed:	mov    %rcx,0x20(%rsp)
  2c46f2:	mov    %rcx,0x68(%rax)
  2c46f6:	inc    %r14
  2c46f9:	cmp    %r14,%r13
  2c46fc:	jne    2c45a2 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x282>
  2c4702:	jmp    2c470d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x3ed>
  2c4704:	mov    0x10(%rsp),%rax
  2c4709:	incq   0x60(%rax)
  2c470d:	add    $0x38,%rsp
  2c4711:	pop    %rbx
  2c4712:	pop    %r12
  2c4714:	pop    %r13
  2c4716:	pop    %r14
  2c4718:	pop    %r15
  2c471a:	pop    %rbp
  2c471b:	ret
  2c471c:	lea    0xff625(%rip),%rdx        # 3c3d48 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0xd8>
  2c4723:	mov    %rdi,%rsi
  2c4726:	call   *0x104ff4(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c472c:	lea    0xff5e5(%rip),%rdx        # 3c3d18 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0xa8>
  2c4733:	mov    %r12,%rdi
  2c4736:	mov    %r12,%rsi
  2c4739:	call   *0x104fe1(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c473f:	lea    0xff5ea(%rip),%rdx        # 3c3d30 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0xc0>
  2c4746:	mov    %rdi,%rsi
  2c4749:	call   *0x104fd1(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c474f:	lea    0xff562(%rip),%rdx        # 3c3cb8 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x48>
  2c4756:	mov    0x28(%rsp),%rdi
  2c475b:	mov    %rdi,%rsi
  2c475e:	call   *0x104fbc(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c4764:	lea    0xff565(%rip),%rdx        # 3c3cd0 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x60>
  2c476b:	mov    0x8(%rsp),%rdi
  2c4770:	mov    %rdi,%rsi
  2c4773:	call   *0x104fa7(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c4779:	dec    %rdi
  2c477c:	lea    0xff565(%rip),%rdx        # 3c3ce8 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x78>
  2c4783:	mov    0x8(%rsp),%rsi
  2c4788:	call   *0x104f92(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c478e:	lea    0xff56b(%rip),%rdx        # 3c3d00 <anon.842e8684e504a58b9e8cc66d76e73277.32.llvm.18170108012268971648+0x90>
  2c4795:	mov    %rcx,%rsi
  2c4798:	call   *0x104f82(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  2c479e:	int3
  2c479f:	int3

