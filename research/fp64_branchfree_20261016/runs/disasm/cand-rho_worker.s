  2031f0:	push   %rbp
  2031f1:	push   %r15
  2031f3:	push   %r14
  2031f5:	push   %r13
  2031f7:	push   %r12
  2031f9:	push   %rbx
  2031fa:	sub    $0x3a8,%rsp
  203201:	mov    %r9,%r13
  203204:	mov    %rdi,%rbp
  203207:	mov    0x3e8(%rsp),%r15
  20320f:	movabs $0x9e3779b97f4a7c15,%rax
  203219:	movabs $0x7fffffffffffffff,%r12
  203223:	mov    %rsi,0x1a0(%rsp)
  20322b:	mov    %rdx,0x1a8(%rsp)
  203233:	mov    %rcx,0x170(%rsp)
  20323b:	mov    %r8,0x178(%rsp)
  203243:	movdqu 0x18(%rdi),%xmm0
  203248:	movdqu 0x28(%rdi),%xmm1
  20324d:	movdqu 0x38(%rdi),%xmm2
  203252:	movdqa %xmm2,0x230(%rsp)
  20325b:	movdqa %xmm1,0x220(%rsp)
  203264:	movdqa %xmm0,0x210(%rsp)
  20326d:	mov    0x58(%rdi),%rbx
  203271:	mov    %rbx,0xb8(%rsp)
  203279:	mov    0x8(%r15),%r14
  20327d:	mov    0x3f0(%rsp),%rsi
  203285:	imul   %rax,%rsi
  203289:	add    %rax,%rsi
  20328c:	xor    0x10(%r15),%rsi
  203290:	lea    0x260(%rsp),%rdi
  203298:	call   26d6d0 <<rand::rngs::std::StdRng as rand_core::SeedableRng>::seed_from_u64>
  20329d:	mov    0x24(%r15),%ecx
  2032a1:	mov    $0x28,%eax
  2032a6:	shl    %cl,%rax
  2032a9:	cmp    $0x101,%rax
  2032af:	mov    $0x100,%esi
  2032b4:	cmovae %rax,%rsi
  2032b8:	mov    $0xffffffff,%eax
  2032bd:	cmp    %rax,%rsi
  2032c0:	cmovae %rax,%rsi
  2032c4:	lea    0x0(,%r14,8),%r15
  2032cc:	mov    %r14,%rax
  2032cf:	shr    $0x3d,%rax
  2032d3:	setne  %al
  2032d6:	lea    -0x7(%r12),%rdx
  2032db:	cmp    %rdx,%r15
  2032de:	seta   %dl
  2032e1:	or     %al,%dl
  2032e3:	je     2032f4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x104>
  2032e5:	xor    %r13d,%r13d
  2032e8:	mov    %r13,%rdi
  2032eb:	mov    %r15,%rsi
  2032ee:	call   *0x1c5cdc(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  2032f4:	mov    %rsi,0x1c8(%rsp)
  2032fc:	mov    %ecx,0x64(%rsp)
  203300:	mov    %rbx,0xe8(%rsp)
  203308:	mov    %rbp,0xb0(%rsp)
  203310:	mov    %r13,0x88(%rsp)
  203318:	mov    %r14,0x8(%rsp)
  20331d:	test   %r15,%r15
  203320:	je     203418 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x228>
  203326:	call   *0x1c5cac(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  20332c:	mov    $0x8,%r13d
  203332:	mov    $0x8,%esi
  203337:	mov    %r15,%rdi
  20333a:	call   *0x1c6078(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  203340:	test   %rax,%rax
  203343:	mov    0x88(%rsp),%r12
  20334b:	mov    0x64(%rsp),%ebx
  20334f:	je     2032e8 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf8>
  203351:	mov    %rax,%rbp
  203354:	call   *0x1c5c7e(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  20335a:	mov    $0x8,%esi
  20335f:	mov    %r15,%rdi
  203362:	call   *0x1c6050(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  203368:	mov    %rax,0x40(%rsp)
  20336d:	test   %rax,%rax
  203370:	je     204655 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1465>
  203376:	call   *0x1c5c5c(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  20337c:	mov    $0x8,%esi
  203381:	mov    %r15,%rdi
  203384:	call   *0x1c602e(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  20338a:	mov    %rax,0x30(%rsp)
  20338f:	test   %rax,%rax
  203392:	je     204665 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1475>
  203398:	call   *0x1c5c3a(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  20339e:	mov    $0x8,%esi
  2033a3:	mov    %r15,%rdi
  2033a6:	call   *0x1c600c(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  2033ac:	mov    %rax,0x38(%rsp)
  2033b1:	test   %rax,%rax
  2033b4:	je     204675 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1485>
  2033ba:	call   *0x1c5c18(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  2033c0:	mov    $0x8,%esi
  2033c5:	mov    %r15,%rdi
  2033c8:	call   *0x1c5c12(%rip)        # 3c8fe0 <_DYNAMIC+0x2a0>
  2033ce:	mov    %rax,0x28(%rsp)
  2033d3:	test   %rax,%rax
  2033d6:	je     204685 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1495>
  2033dc:	mov    0x8(%rsp),%rcx
  2033e1:	mov    %rcx,%rdx
  2033e4:	cmp    $0x2,%rcx
  2033e8:	mov    %rbp,0x18(%rsp)
  2033ed:	mov    %rdx,0x108(%rsp)
  2033f5:	jb     20345c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x26c>
  2033f7:	lea    -0x8(%r15),%rdx
  2033fb:	mov    0x28(%rsp),%r14
  203400:	mov    %r14,%rdi
  203403:	mov    $0xff,%esi
  203408:	call   *0x1c63f2(%rip)        # 3c9800 <memset@GLIBC_2.2.5>
  20340e:	lea    (%r14,%r15,1),%rax
  203412:	add    $0xfffffffffffffff8,%rax
  203416:	jmp    203472 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x282>
  203418:	mov    $0x8,%eax
  20341d:	mov    %rax,0x38(%rsp)
  203422:	xor    %edx,%edx
  203424:	mov    %rax,0x40(%rsp)
  203429:	mov    $0x8,%ebp
  20342e:	mov    %rax,0x30(%rsp)
  203433:	mov    %rax,0x28(%rsp)
  203438:	mov    0x8(%rsp),%rcx
  20343d:	mov    0x88(%rsp),%r12
  203445:	mov    0x64(%rsp),%ebx
  203449:	cmp    $0x2,%rcx
  20344d:	mov    %rbp,0x18(%rsp)
  203452:	mov    %rdx,0x108(%rsp)
  20345a:	jae    2033f7 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x207>
  20345c:	mov    0x28(%rsp),%rax
  203461:	test   %rcx,%rcx
  203464:	mov    0x3e8(%rsp),%r14
  20346c:	je     203581 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x391>
  203472:	movq   $0xffffffffffffffff,(%rax)
  203479:	mov    0x8(%rsp),%rax
  20347e:	lea    0x0(,%rax,4),%r13
  203486:	call   *0x1c5b4c(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  20348c:	mov    $0x4,%esi
  203491:	mov    %r13,%rdi
  203494:	call   *0x1c5f1e(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  20349a:	mov    %rax,0x48(%rsp)
  20349f:	test   %rax,%rax
  2034a2:	je     2045f2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1402>
  2034a8:	test   %r15,%r15
  2034ab:	je     203500 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x310>
  2034ad:	call   *0x1c5b25(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  2034b3:	mov    $0x8,%esi
  2034b8:	mov    %r15,%rdi
  2034bb:	call   *0x1c5ef7(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  2034c1:	mov    %rax,0x58(%rsp)
  2034c6:	test   %rax,%rax
  2034c9:	je     204695 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14a5>
  2034cf:	call   *0x1c5b03(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  2034d5:	mov    $0x8,%esi
  2034da:	mov    %r15,%rdi
  2034dd:	call   *0x1c5ed5(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  2034e3:	mov    %rax,0x78(%rsp)
  2034e8:	test   %rax,%rax
  2034eb:	je     2046a5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14b5>
  2034f1:	mov    0x8(%rsp),%r14
  2034f6:	mov    %r14,0xa0(%rsp)
  2034fe:	jmp    203520 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x330>
  203500:	mov    $0x8,%eax
  203505:	mov    %rax,0x78(%rsp)
  20350a:	movq   $0x0,0xa0(%rsp)
  203516:	mov    %rax,0x58(%rsp)
  20351b:	mov    0x8(%rsp),%r14
  203520:	call   *0x1c5ab2(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  203526:	mov    $0x1,%esi
  20352b:	mov    %r14,%rdi
  20352e:	call   *0x1c5e84(%rip)        # 3c93b8 <_DYNAMIC+0x678>
  203534:	mov    %rax,0x98(%rsp)
  20353c:	test   %rax,%rax
  20353f:	je     204605 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1415>
  203545:	test   %r15,%r15
  203548:	mov    0x3e8(%rsp),%r14
  203550:	je     203578 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x388>
  203552:	call   *0x1c5a80(%rip)        # 3c8fd8 <_DYNAMIC+0x298>
  203558:	mov    $0x8,%esi
  20355d:	mov    %r15,%rdi
  203560:	call   *0x1c5a7a(%rip)        # 3c8fe0 <_DYNAMIC+0x2a0>
  203566:	test   %rax,%rax
  203569:	je     2046b5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14c5>
  20356f:	xor    %edx,%edx
  203571:	mov    0x8(%rsp),%rcx
  203576:	jmp    2035bc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x3cc>
  203578:	xor    %edx,%edx
  20357a:	mov    $0x8,%eax
  20357f:	jmp    2035ba <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x3ca>
  203581:	mov    $0x8,%eax
  203586:	mov    %rax,0x78(%rsp)
  20358b:	mov    $0x1,%dl
  20358d:	mov    $0x4,%eax
  203592:	mov    %rax,0x48(%rsp)
  203597:	mov    $0x1,%eax
  20359c:	mov    %rax,0x98(%rsp)
  2035a4:	movq   $0x0,0xa0(%rsp)
  2035b0:	mov    $0x8,%eax
  2035b5:	mov    %rax,0x58(%rsp)
  2035ba:	xor    %ecx,%ecx
  2035bc:	mov    %rcx,0x158(%rsp)
  2035c4:	mov    %rax,0x160(%rsp)
  2035cc:	movq   $0x0,0x168(%rsp)
  2035d8:	lea    0xb8(%rsp),%rax
  2035e0:	mov    %rax,0x1e8(%rsp)
  2035e8:	mov    0xb0(%rsp),%rax
  2035f0:	mov    %rax,0x1f0(%rsp)
  2035f8:	lea    0x1a0(%rsp),%rax
  203600:	mov    %rax,0x1f8(%rsp)
  203608:	lea    0x170(%rsp),%rax
  203610:	mov    %rax,0x200(%rsp)
  203618:	mov    %r14,0x208(%rsp)
  203620:	mov    %rdx,0x20(%rsp)
  203625:	test   %dl,%dl
  203627:	mov    0x8(%rsp),%r14
  20362c:	jne    203688 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x498>
  20362e:	xor    %r15d,%r15d
  203631:	lea    0x260(%rsp),%r13
  203639:	nopl   0x0(%rax)
  203640:	sub    $0x8,%rsp
  203644:	lea    0x1f0(%rsp),%rdi
  20364c:	mov    %r15,%rsi
  20364f:	mov    %r13,%rdx
  203652:	mov    %rbp,%rcx
  203655:	mov    %r14,%r8
  203658:	mov    0x48(%rsp),%r9
  20365d:	push   %r14
  20365f:	push   0x58(%rsp)
  203663:	push   %r14
  203665:	push   0x48(%rsp)
  203669:	push   %r14
  20366b:	push   0x68(%rsp)
  20366f:	push   %r14
  203671:	push   0x70(%rsp)
  203675:	push   %r14
  203677:	call   201940 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  20367c:	add    $0x50,%rsp
  203680:	inc    %r15
  203683:	cmp    %r15,%r14
  203686:	jne    203640 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x450>
  203688:	mov    0x3f8(%rsp),%r13
  203690:	movzbl 0x60(%r13),%eax
  203695:	test   %al,%al
  203697:	je     2037dc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x5ec>
  20369d:	movq   $0x0,0x80(%rsp)
  2036a9:	movq   $0x0,0x68(%rsp)
  2036b2:	mov    0x8(%rsp),%rbx
  2036b7:	mov    0x20(%rsp),%r14
  2036bc:	mov    0x80(%rsp),%rax
  2036c4:	lock add %rax,0x50(%r13)
  2036c9:	mov    0x68(%rsp),%rax
  2036ce:	lock add %rax,0x58(%r13)
  2036d3:	mov    0x158(%rsp),%rsi
  2036db:	test   %rsi,%rsi
  2036de:	je     2036f7 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x507>
  2036e0:	mov    0x160(%rsp),%rdi
  2036e8:	shl    $0x3,%rsi
  2036ec:	mov    $0x8,%edx
  2036f1:	call   *0x1c58b9(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  2036f7:	test   %r14b,%r14b
  2036fa:	jne    203712 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x522>
  2036fc:	mov    $0x1,%edx
  203701:	mov    0x98(%rsp),%rdi
  203709:	mov    %rbx,%rsi
  20370c:	call   *0x1c589e(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  203712:	mov    0xa0(%rsp),%r14
  20371a:	test   %r14,%r14
  20371d:	je     20374d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x55d>
  20371f:	shl    $0x3,%r14
  203723:	mov    0x1c5886(%rip),%rbx        # 3c8fb0 <_DYNAMIC+0x270>
  20372a:	mov    $0x8,%edx
  20372f:	mov    0x78(%rsp),%rdi
  203734:	mov    %r14,%rsi
  203737:	call   *%rbx
  203739:	mov    $0x8,%edx
  20373e:	mov    0x58(%rsp),%rdi
  203743:	mov    %r14,%rsi
  203746:	call   *%rbx
  203748:	mov    0x8(%rsp),%rbx
  20374d:	test   %rbx,%rbx
  203750:	je     203769 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x579>
  203752:	shl    $0x2,%rbx
  203756:	mov    $0x4,%edx
  20375b:	mov    0x48(%rsp),%rdi
  203760:	mov    %rbx,%rsi
  203763:	call   *0x1c5847(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  203769:	mov    0x108(%rsp),%r14
  203771:	test   %r14,%r14
  203774:	je     2037ca <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x5da>
  203776:	shl    $0x3,%r14
  20377a:	mov    0x1c582f(%rip),%rbx        # 3c8fb0 <_DYNAMIC+0x270>
  203781:	mov    $0x8,%edx
  203786:	mov    0x28(%rsp),%rdi
  20378b:	mov    %r14,%rsi
  20378e:	call   *%rbx
  203790:	mov    $0x8,%edx
  203795:	mov    0x38(%rsp),%rdi
  20379a:	mov    %r14,%rsi
  20379d:	call   *%rbx
  20379f:	mov    $0x8,%edx
  2037a4:	mov    0x30(%rsp),%rdi
  2037a9:	mov    %r14,%rsi
  2037ac:	call   *%rbx
  2037ae:	mov    $0x8,%edx
  2037b3:	mov    0x40(%rsp),%rdi
  2037b8:	mov    %r14,%rsi
  2037bb:	call   *%rbx
  2037bd:	mov    $0x8,%edx
  2037c2:	mov    %rbp,%rdi
  2037c5:	mov    %r14,%rsi
  2037c8:	call   *%rbx
  2037ca:	add    $0x3a8,%rsp
  2037d1:	pop    %rbx
  2037d2:	pop    %r12
  2037d4:	pop    %r13
  2037d6:	pop    %r14
  2037d8:	pop    %r15
  2037da:	pop    %rbp
  2037db:	ret
  2037dc:	mov    0x3e0(%rsp),%rax
  2037e4:	mov    %rax,0xe0(%rsp)
  2037ec:	mov    0x3e8(%rsp),%rcx
  2037f4:	mov    0x18(%rcx),%rax
  2037f8:	mov    %rax,0x1d8(%rsp)
  203800:	xor    %eax,%eax
  203802:	sub    0x20(%rcx),%eax
  203805:	and    $0x3f,%eax
  203808:	mov    %rax,0x1e0(%rsp)
  203810:	movzbl 0x28(%rcx),%eax
  203814:	mov    %al,0x17(%rsp)
  203818:	mov    %ebx,%eax
  20381a:	neg    %eax
  20381c:	and    $0x3f,%eax
  20381f:	mov    %rax,0x1d0(%rsp)
  203827:	lea    0x8(%r13),%rax
  20382b:	mov    %rax,0x1b8(%rsp)
  203833:	lea    0x28(%r13),%rax
  203837:	mov    %rax,0x1c0(%rsp)
  20383f:	movq   $0x0,0x80(%rsp)
  20384b:	movq   $0x0,0x68(%rsp)
  203854:	mov    0x8(%rsp),%rbx
  203859:	mov    0x20(%rsp),%r14
  20385e:	jmp    203899 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x6a9>
  203860:	mov    0x3f8(%rsp),%r13
  203868:	lock add %rax,0x50(%r13)
  20386d:	mov    0x68(%rsp),%rax
  203872:	lock add %rax,0x58(%r13)
  203877:	movq   $0x0,0x68(%rsp)
  203880:	movq   $0x0,0x80(%rsp)
  20388c:	movzbl 0x60(%r13),%eax
  203891:	test   %al,%al
  203893:	jne    2036bc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4cc>
  203899:	mov    0x1d8(%rsp),%rcx
  2038a1:	test   %rcx,%rcx
  2038a4:	je     2038b3 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x6c3>
  2038a6:	mov    0x50(%r13),%rax
  2038aa:	cmp    %rcx,%rax
  2038ad:	ja     2036bc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4cc>
  2038b3:	test   %r14b,%r14b
  2038b6:	jne    203942 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x752>
  2038bc:	mov    0x210(%rsp),%rax
  2038c4:	mov    0x220(%rsp),%rdx
  2038cc:	xor    %esi,%esi
  2038ce:	xchg   %ax,%ax
  2038d0:	mov    0x0(%rbp,%rsi,8),%r13
  2038d5:	movabs $0x9e3779b97f4a7c15,%rcx
  2038df:	imul   %rcx,%r13
  2038e3:	mov    0x1e0(%rsp),%rcx
  2038eb:	shr    %cl,%r13
  2038ee:	mov    0x78(%rsp),%rcx
  2038f3:	mov    %r13,(%rcx,%rsi,8)
  2038f7:	cmp    0x3e0(%rsp),%r13
  2038ff:	jae    204539 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1349>
  203905:	shl    $0x5,%r13
  203909:	mov    (%r12,%r13,1),%rcx
  20390d:	sub    0x0(%rbp,%rsi,8),%rcx
  203912:	mov    $0x0,%edi
  203917:	cmovb  %rax,%rdi
  20391b:	add    %rdi,%rcx
  20391e:	test   %rcx,%rcx
  203921:	cmove  %rdx,%rcx
  203925:	mov    0x98(%rsp),%rdi
  20392d:	sete   (%rdi,%rsi,1)
  203931:	mov    0x58(%rsp),%rdi
  203936:	mov    %rcx,(%rdi,%rsi,8)
  20393a:	inc    %rsi
  20393d:	cmp    %rsi,%rbx
  203940:	jne    2038d0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x6e0>
  203942:	lea    0x210(%rsp),%rdi
  20394a:	mov    0x58(%rsp),%rsi
  20394f:	mov    %rbx,%rdx
  203952:	lea    0x158(%rsp),%rcx
  20395a:	call   *0x1c5c80(%rip)        # 3c95e0 <_DYNAMIC+0x8a0>
  203960:	test   %r14b,%r14b
  203963:	jne    204405 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1215>
  203969:	xor    %r13d,%r13d
  20396c:	mov    0x210(%rsp),%r14
  203974:	mov    0x218(%rsp),%rax
  20397c:	mov    %rax,0x50(%rsp)
  203981:	mov    0x238(%rsp),%rax
  203989:	mov    %rax,0xa8(%rsp)
  203991:	mov    0x40(%rsp),%rax
  203996:	lea    (%rax,%r13,8),%rax
  20399a:	mov    %rax,0x90(%rsp)
  2039a2:	lea    0x0(%rbp,%r13,8),%r15
  2039a7:	mov    0x30(%rsp),%rax
  2039ac:	lea    (%rax,%r13,8),%rax
  2039b0:	mov    %rax,0xf0(%rsp)
  2039b8:	mov    0x38(%rsp),%rax
  2039bd:	lea    (%rax,%r13,8),%rax
  2039c1:	mov    %rax,0x70(%rsp)
  2039c6:	jmp    203a38 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x848>
  2039c8:	nopl   0x0(%rax,%rax,1)
  2039d0:	sub    $0x8,%rsp
  2039d4:	lea    0x1f0(%rsp),%rdi
  2039dc:	mov    %r13,%rsi
  2039df:	lea    0x268(%rsp),%rdx
  2039e7:	mov    %rbp,%rcx
  2039ea:	mov    %rbx,%r8
  2039ed:	mov    0x48(%rsp),%r9
  2039f2:	push   %rbx
  2039f3:	push   0x58(%rsp)
  2039f7:	push   %rbx
  2039f8:	push   0x48(%rsp)
  2039fc:	push   %rbx
  2039fd:	push   0x68(%rsp)
  203a01:	push   %rbx
  203a02:	push   0x70(%rsp)
  203a06:	push   %rbx
  203a07:	call   201940 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  203a0c:	add    $0x50,%rsp
  203a10:	inc    %r13
  203a13:	addq   $0x8,0x90(%rsp)
  203a1c:	add    $0x8,%r15
  203a20:	addq   $0x8,0xf0(%rsp)
  203a29:	addq   $0x8,0x70(%rsp)
  203a2f:	cmp    %rbx,%r13
  203a32:	jae    204400 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1210>
  203a38:	cmp    %rbx,%r13
  203a3b:	jae    204542 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1352>
  203a41:	mov    0x98(%rsp),%rax
  203a49:	cmpb   $0x0,(%rax,%r13,1)
  203a4e:	jne    2039d0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x7e0>
  203a50:	mov    0x78(%rsp),%rax
  203a55:	mov    (%rax,%r13,8),%r10
  203a59:	cmp    0x3e0(%rsp),%r10
  203a61:	jae    20456e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x137e>
  203a67:	shl    $0x5,%r10
  203a6b:	mov    0x40(%rsp),%rax
  203a70:	mov    (%rax,%r13,8),%r9
  203a74:	mov    0x8(%r12,%r10,1),%rax
  203a79:	sub    %r9,%rax
  203a7c:	mov    $0x0,%ecx
  203a81:	cmovb  %r14,%rcx
  203a85:	add    %rcx,%rax
  203a88:	mov    0x58(%rsp),%rcx
  203a8d:	mulq   (%rcx,%r13,8)
  203a91:	mov    %rax,%rsi
  203a94:	mov    %rdx,%rdi
  203a97:	mov    0x50(%rsp),%rbx
  203a9c:	mov    %rbx,%rax
  203a9f:	imul   %rsi,%rax
  203aa3:	mul    %r14
  203aa6:	mov    %rdx,%rcx
  203aa9:	add    %rsi,%rax
  203aac:	adc    %rdi,%rcx
  203aaf:	cmp    %rcx,%r14
  203ab2:	mov    %r14,%rax
  203ab5:	mov    $0x0,%r11d
  203abb:	cmova  %r11,%rax
  203abf:	sub    %rax,%rcx
  203ac2:	mov    %rcx,%rax
  203ac5:	mul    %rcx
  203ac8:	mov    %rax,%rsi
  203acb:	mov    %rdx,%rdi
  203ace:	mov    %rbx,%rax
  203ad1:	imul   %rsi,%rax
  203ad5:	mul    %r14
  203ad8:	add    %rsi,%rax
  203adb:	adc    %rdi,%rdx
  203ade:	cmp    %rdx,%r14
  203ae1:	mov    %r14,%rax
  203ae4:	cmova  %r11,%rax
  203ae8:	sub    %rax,%rdx
  203aeb:	mov    0x0(%rbp,%r13,8),%rdi
  203af0:	sub    %rdi,%rdx
  203af3:	mov    $0x0,%eax
  203af8:	cmovb  %r14,%rax
  203afc:	add    %rdx,%rax
  203aff:	sub    (%r12,%r10,1),%rax
  203b03:	mov    $0x0,%edx
  203b08:	cmovb  %r14,%rdx
  203b0c:	lea    (%rdx,%rax,1),%r8
  203b10:	mov    %rdi,%rax
  203b13:	sub    %r8,%rax
  203b16:	mov    $0x0,%edx
  203b1b:	cmovb  %r14,%rdx
  203b1f:	add    %rdx,%rax
  203b22:	mul    %rcx
  203b25:	mov    %rax,%rcx
  203b28:	mov    %rdx,%rsi
  203b2b:	mov    %rbx,%rax
  203b2e:	imul   %rcx,%rax
  203b32:	mul    %r14
  203b35:	add    %rcx,%rax
  203b38:	adc    %rsi,%rdx
  203b3b:	cmp    %rdx,%r14
  203b3e:	mov    %r14,%rax
  203b41:	cmova  %r11,%rax
  203b45:	sub    %rax,%rdx
  203b48:	sub    %r9,%rdx
  203b4b:	mov    $0x0,%eax
  203b50:	cmovb  %r14,%rax
  203b54:	lea    (%rax,%rdx,1),%rcx
  203b58:	mov    0x30(%rsp),%rax
  203b5d:	mov    (%rax,%r13,8),%rax
  203b61:	mov    0x10(%r12,%r10,1),%rbx
  203b66:	add    %rax,%rbx
  203b69:	mov    0xb8(%rsp),%rbp
  203b71:	cmp    %rbp,%rbx
  203b74:	mov    %rbp,%rdx
  203b77:	cmovb  %r11,%rdx
  203b7b:	sub    %rdx,%rbx
  203b7e:	mov    0x38(%rsp),%rdx
  203b83:	mov    (%rdx,%r13,8),%rdx
  203b87:	mov    0x18(%r12,%r10,1),%r12
  203b8c:	add    %rdx,%r12
  203b8f:	cmp    %rbp,%r12
  203b92:	mov    %rbp,%rsi
  203b95:	cmovb  %r11,%rsi
  203b99:	sub    %rsi,%r12
  203b9c:	cmpb   $0x0,0x17(%rsp)
  203ba1:	je     203bd0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x9e0>
  203ba3:	cmp    0xa8(%rsp),%rcx
  203bab:	jbe    203bd0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x9e0>
  203bad:	mov    %r14,%rsi
  203bb0:	sub    %rcx,%rsi
  203bb3:	mov    %rbp,%rcx
  203bb6:	sub    %rbx,%rcx
  203bb9:	test   %rbx,%rbx
  203bbc:	cmovne %rcx,%rbx
  203bc0:	mov    %rbp,%rcx
  203bc3:	sub    %r12,%rcx
  203bc6:	test   %r12,%r12
  203bc9:	cmovne %rcx,%r12
  203bcd:	mov    %rsi,%rcx
  203bd0:	cmp    0x8(%rsp),%r13
  203bd5:	jae    204553 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1363>
  203bdb:	mov    0x28(%rsp),%rsi
  203be0:	cmp    (%rsi,%r13,8),%r8
  203be4:	jne    203c80 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xa90>
  203bea:	cmp    %r8,%rdi
  203bed:	jae    203bfb <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xa0b>
  203bef:	mov    %rdi,%r8
  203bf2:	mov    %r9,%rcx
  203bf5:	mov    %rax,%rbx
  203bf8:	mov    %rdx,%r12
  203bfb:	lea    0x110(%rsp),%rdi
  203c03:	mov    0xb0(%rsp),%rsi
  203c0b:	mov    %r8,%rdx
  203c0e:	call   *0x1c59dc(%rip)        # 3c95f0 <_DYNAMIC+0x8b0>
  203c14:	cmpb   $0x0,0x110(%rsp)
  203c1c:	jne    203cc0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xad0>
  203c22:	sub    $0x8,%rsp
  203c26:	lea    0x1f0(%rsp),%rdi
  203c2e:	mov    %r13,%rsi
  203c31:	lea    0x268(%rsp),%rdx
  203c39:	mov    0x20(%rsp),%rbp
  203c3e:	mov    %rbp,%rcx
  203c41:	mov    0x10(%rsp),%rbx
  203c46:	mov    %rbx,%r8
  203c49:	mov    0x48(%rsp),%r9
  203c4e:	push   %rbx
  203c4f:	push   0x58(%rsp)
  203c53:	push   %rbx
  203c54:	push   0x48(%rsp)
  203c58:	push   %rbx
  203c59:	push   0x68(%rsp)
  203c5d:	push   %rbx
  203c5e:	push   0x70(%rsp)
  203c62:	push   %rbx
  203c63:	call   201940 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  203c68:	add    $0x50,%rsp
  203c6c:	mov    0x88(%rsp),%r12
  203c74:	jmp    203a10 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x820>
  203c79:	nopl   0x0(%rax)
  203c80:	mov    0x28(%rsp),%rax
  203c85:	mov    %rdi,(%rax,%r13,8)
  203c89:	mov    0x18(%rsp),%rbp
  203c8e:	mov    %r8,0x0(%rbp,%r13,8)
  203c93:	mov    0x40(%rsp),%rax
  203c98:	mov    %rcx,(%rax,%r13,8)
  203c9c:	mov    0x30(%rsp),%rax
  203ca1:	mov    %rbx,(%rax,%r13,8)
  203ca5:	mov    0x38(%rsp),%rax
  203caa:	mov    %r12,(%rax,%r13,8)
  203cae:	mov    0x8(%rsp),%rbx
  203cb3:	mov    0x20(%rsp),%r14
  203cb8:	jmp    203d83 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb93>
  203cbd:	nopl   (%rax)
  203cc0:	mov    0x118(%rsp),%rdx
  203cc8:	mov    0x120(%rsp),%rsi
  203cd0:	add    %rbx,%rbx
  203cd3:	cmp    %rbp,%rbx
  203cd6:	mov    %rbp,%rax
  203cd9:	mov    $0x0,%edi
  203cde:	cmovb  %rdi,%rax
  203ce2:	mov    %rbx,%rcx
  203ce5:	sub    %rax,%rcx
  203ce8:	add    %r12,%r12
  203ceb:	cmp    %rbp,%r12
  203cee:	cmovb  %rdi,%rbp
  203cf2:	sub    %rbp,%r12
  203cf5:	mov    0x18(%rsp),%rbp
  203cfa:	mov    %rdx,0x0(%rbp,%r13,8)
  203cff:	cmp    0xa8(%rsp),%rsi
  203d07:	seta   %dl
  203d0a:	sub    %rsi,%r14
  203d0d:	test   %dl,0x17(%rsp)
  203d11:	cmove  %rsi,%r14
  203d15:	mov    0x40(%rsp),%rdx
  203d1a:	mov    %r14,(%rdx,%r13,8)
  203d1e:	je     203d57 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb67>
  203d20:	mov    0xb8(%rsp),%rdx
  203d28:	sub    %rcx,%rdx
  203d2b:	cmp    %rax,%rbx
  203d2e:	cmove  %rdi,%rdx
  203d32:	mov    0x30(%rsp),%rax
  203d37:	mov    %rdx,(%rax,%r13,8)
  203d3b:	mov    $0x0,%eax
  203d40:	test   %r12,%r12
  203d43:	mov    0x20(%rsp),%r14
  203d48:	je     203d68 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb78>
  203d4a:	mov    0xb8(%rsp),%rax
  203d52:	sub    %r12,%rax
  203d55:	jmp    203d68 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb78>
  203d57:	mov    0x30(%rsp),%rax
  203d5c:	mov    %rcx,(%rax,%r13,8)
  203d60:	mov    %r12,%rax
  203d63:	mov    0x20(%rsp),%r14
  203d68:	mov    0x38(%rsp),%rcx
  203d6d:	mov    %rax,(%rcx,%r13,8)
  203d71:	mov    0x28(%rsp),%rax
  203d76:	movq   $0xffffffffffffffff,(%rax,%r13,8)
  203d7e:	mov    0x8(%rsp),%rbx
  203d83:	mov    0x48(%rsp),%rcx
  203d88:	mov    (%rcx,%r13,4),%eax
  203d8c:	inc    %eax
  203d8e:	mov    %eax,(%rcx,%r13,4)
  203d92:	cmpl   $0x0,0x64(%rsp)
  203d97:	mov    0x88(%rsp),%r12
  203d9f:	je     203e40 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xc50>
  203da5:	mov    (%r15),%rdx
  203da8:	rol    $0x17,%rdx
  203dac:	movabs $0xd6e8feb86659fd93,%rcx
  203db6:	imul   %rcx,%rdx
  203dba:	mov    0x1d0(%rsp),%rcx
  203dc2:	shr    %cl,%rdx
  203dc5:	test   %rdx,%rdx
  203dc8:	je     203e40 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xc50>
  203dca:	cmp    0x1c8(%rsp),%eax
  203dd1:	jbe    203e20 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xc30>
  203dd3:	sub    $0x8,%rsp
  203dd7:	lea    0x1f0(%rsp),%rdi
  203ddf:	mov    %r13,%rsi
  203de2:	lea    0x268(%rsp),%rdx
  203dea:	mov    %rbp,%rcx
  203ded:	mov    %rbx,%r8
  203df0:	mov    0x48(%rsp),%r9
  203df5:	push   %rbx
  203df6:	push   0x58(%rsp)
  203dfa:	push   %rbx
  203dfb:	push   0x48(%rsp)
  203dff:	push   %rbx
  203e00:	push   0x68(%rsp)
  203e04:	push   %rbx
  203e05:	push   0x70(%rsp)
  203e09:	push   %rbx
  203e0a:	call   201940 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  203e0f:	add    $0x50,%rsp
  203e13:	data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  203e20:	inc    %r13
  203e23:	cmp    %rbx,%r13
  203e26:	jb     20396c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x77c>
  203e2c:	jmp    204405 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1215>
  203e31:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  203e40:	xor    %eax,%eax
  203e42:	mov    0x3f8(%rsp),%rcx
  203e4a:	mov    $0x1,%edx
  203e4f:	lock cmpxchg %edx,(%rcx)
  203e53:	jne    20438e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x119e>
  203e59:	mov    0x1c5848(%rip),%rax        # 3c96a8 <_DYNAMIC+0x968>
  203e60:	mov    (%rax),%rax
  203e63:	movabs $0x7fffffffffffffff,%rcx
  203e6d:	test   %rcx,%rax
  203e70:	jne    2043a1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11b1>
  203e76:	movl   $0x0,0x50(%rsp)
  203e7e:	mov    0x3f8(%rsp),%rcx
  203e86:	movzbl 0x4(%rcx),%eax
  203e8a:	test   %al,%al
  203e8c:	jne    2044f9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1309>
  203e92:	mov    0x68(%rsp),%rbx
  203e97:	inc    %rbx
  203e9a:	cmpq   $0x0,0x20(%rcx)
  203e9f:	lea    0xc0(%rsp),%r12
  203ea7:	mov    %rbx,0x68(%rsp)
  203eac:	je     203f5d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xd6d>
  203eb2:	mov    0x1c0(%rsp),%rdi
  203eba:	mov    %r15,%rsi
  203ebd:	call   250110 <<std::hash::random::RandomState as core::hash::BuildHasher>::hash_one::<&crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolMono>>
  203ec2:	mov    %rax,%rsi
  203ec5:	shr    $0x39,%rsi
  203ec9:	mov    0x3f8(%rsp),%rdx
  203ed1:	mov    0x8(%rdx),%rcx
  203ed5:	mov    0x10(%rdx),%rdx
  203ed9:	movd   %esi,%xmm0
  203edd:	punpcklbw %xmm0,%xmm0
  203ee1:	pshuflw $0x0,%xmm0,%xmm0
  203ee6:	pshufd $0x44,%xmm0,%xmm0
  203eeb:	xor    %esi,%esi
  203eed:	and    %rdx,%rax
  203ef0:	movdqu (%rcx,%rax,1),%xmm1
  203ef5:	movdqa %xmm1,%xmm2
  203ef9:	pcmpeqb %xmm0,%xmm2
  203efd:	pmovmskb %xmm2,%edi
  203f01:	test   %edi,%edi
  203f03:	je     203f40 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xd50>
  203f05:	mov    (%r15),%r8
  203f08:	tzcnt  %edi,%r10d
  203f0d:	add    %rax,%r10
  203f10:	and    %rdx,%r10
  203f13:	shl    $0x5,%r10
  203f17:	mov    %rcx,%r9
  203f1a:	sub    %r10,%r9
  203f1d:	cmp    -0x20(%r9),%r8
  203f21:	je     203fb8 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xdc8>
  203f27:	lea    -0x1(%rdi),%r9d
  203f2b:	and    %di,%r9w
  203f2f:	mov    %r9d,%edi
  203f32:	jne    203f08 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xd18>
  203f34:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  203f40:	pcmpeqd %xmm2,%xmm2
  203f44:	pcmpeqb %xmm2,%xmm1
  203f48:	pmovmskb %xmm1,%edi
  203f4c:	test   %edi,%edi
  203f4e:	jne    203f5d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xd6d>
  203f50:	add    %rsi,%rax
  203f53:	add    $0x10,%rax
  203f57:	add    $0x10,%rsi
  203f5b:	jmp    203eed <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xcfd>
  203f5d:	mov    (%r15),%rdx
  203f60:	mov    0x90(%rsp),%rax
  203f68:	mov    (%rax),%rax
  203f6b:	mov    0xf0(%rsp),%rcx
  203f73:	mov    (%rcx),%rcx
  203f76:	mov    0x70(%rsp),%rsi
  203f7b:	mov    (%rsi),%rsi
  203f7e:	mov    %rax,0xc0(%rsp)
  203f86:	mov    %rcx,0xc8(%rsp)
  203f8e:	mov    %rsi,0xd0(%rsp)
  203f96:	lea    0x110(%rsp),%rdi
  203f9e:	mov    0x1b8(%rsp),%rsi
  203fa6:	mov    %r12,%rcx
  203fa9:	call   29d900 <<hashbrown::map::HashMap<u64, (u64, u64, u64), std::hash::random::RandomState>>::insert>
  203fae:	mov    0x20(%rsp),%r14
  203fb3:	jmp    204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  203fb8:	mov    -0x18(%r9),%r12
  203fbc:	mov    -0x10(%r9),%rbx
  203fc0:	mov    -0x8(%r9),%rcx
  203fc4:	mov    0xf0(%rsp),%rax
  203fcc:	mov    (%rax),%rbp
  203fcf:	mov    0x70(%rsp),%rax
  203fd4:	mov    (%rax),%r15
  203fd7:	cmp    %rbx,%rbp
  203fda:	jne    203fe5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xdf5>
  203fdc:	cmp    %r15,%rcx
  203fdf:	je     204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  203fe5:	mov    %rcx,0x70(%rsp)
  203fea:	mov    0x90(%rsp),%rax
  203ff2:	mov    (%rax),%r14
  203ff5:	movdqa 0x1a0(%rsp),%xmm0
  203ffe:	movdqa %xmm0,0xf0(%rsp)
  204007:	mov    0x170(%rsp),%rax
  20400f:	mov    %rax,0x90(%rsp)
  204017:	mov    0x178(%rsp),%rax
  20401f:	mov    %rax,0xa8(%rsp)
  204027:	lea    0x110(%rsp),%rdi
  20402f:	mov    0xe8(%rsp),%rsi
  204037:	call   *0x1c5593(%rip)        # 3c95d0 <_DYNAMIC+0x890>
  20403d:	cmp    %r12,%r14
  204040:	jne    20408d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xe9d>
  204042:	sub    %rbx,%rbp
  204045:	mov    0x110(%rsp),%rax
  20404d:	mov    $0x0,%ecx
  204052:	cmovb  %rax,%rcx
  204056:	add    %rbp,%rcx
  204059:	mov    0x70(%rsp),%rsi
  20405e:	sub    %r15,%rsi
  204061:	mov    $0x0,%edx
  204066:	cmovae %rdx,%rax
  20406a:	lea    (%rsi,%rax,1),%r15
  20406e:	mov    0xe8(%rsp),%r11
  204076:	mov    0x20(%rsp),%r14
  20407b:	test   %r15,%r15
  20407e:	lea    0xc0(%rsp),%rbp
  204086:	jne    2040dd <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xeed>
  204088:	jmp    204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  20408d:	add    %rbx,%rbp
  204090:	mov    0xe8(%rsp),%r11
  204098:	mov    %r11,%rax
  20409b:	sub    %rbp,%rax
  20409e:	mov    %r11,%rcx
  2040a1:	mov    $0x0,%edx
  2040a6:	cmova  %rdx,%rcx
  2040aa:	cmp    %rcx,%rbp
  2040ad:	lea    (%rax,%rcx,1),%rcx
  2040b1:	cmove  %rdx,%rcx
  2040b5:	add    0x70(%rsp),%r15
  2040ba:	cmp    %r11,%r15
  2040bd:	mov    %r11,%rax
  2040c0:	cmovb  %rdx,%rax
  2040c4:	sub    %rax,%r15
  2040c7:	mov    0x20(%rsp),%r14
  2040cc:	test   %r15,%r15
  2040cf:	lea    0xc0(%rsp),%rbp
  2040d7:	je     204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  2040dd:	cmp    $0x1,%r15
  2040e1:	sete   %dil
  2040e5:	mov    $0x1,%edx
  2040ea:	cmp    $0x1,%r11
  2040ee:	je     2041dd <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xfed>
  2040f4:	mov    $0x0,%eax
  2040f9:	cmp    $0x1,%r15
  2040fd:	mov    0xb0(%rsp),%rbx
  204105:	je     2041e7 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xff7>
  20410b:	mov    %r11,%rsi
  20410e:	xor    %eax,%eax
  204110:	test   $0x1,%r15b
  204114:	jne    20414f <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf5f>
  204116:	cs nopw 0x0(%rax,%rax,1)
  204120:	mov    %r15,%rdi
  204123:	shr    $1,%rdi
  204126:	mov    %edx,%r8d
  204129:	mov    %rdx,%r9
  20412c:	shr    $1,%r9
  20412f:	mov    %rdx,%r10
  204132:	and    %r11,%r10
  204135:	xor    %r11,%rdx
  204138:	shr    $1,%rdx
  20413b:	add    %r10,%rdx
  20413e:	test   $0x1,%r8b
  204142:	cmove  %r9,%rdx
  204146:	test   $0x2,%r15b
  20414a:	mov    %rdi,%r15
  20414d:	je     204120 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf30>
  20414f:	test   $0x1,%sil
  204153:	jne    20418f <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf9f>
  204155:	data16 cs nopw 0x0(%rax,%rax,1)
  204160:	mov    %rsi,%rdi
  204163:	shr    $1,%rdi
  204166:	mov    %eax,%r8d
  204169:	mov    %rax,%r9
  20416c:	shr    $1,%r9
  20416f:	mov    %rax,%r10
  204172:	and    %r11,%r10
  204175:	xor    %r11,%rax
  204178:	shr    $1,%rax
  20417b:	add    %r10,%rax
  20417e:	test   $0x1,%r8b
  204182:	cmove  %r9,%rax
  204186:	test   $0x2,%sil
  20418a:	mov    %rdi,%rsi
  20418d:	je     204160 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf70>
  20418f:	mov    %r15,%rdi
  204192:	sub    %rsi,%rdi
  204195:	jae    2041b0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xfc0>
  204197:	sub    %r15,%rsi
  20419a:	sub    %rdx,%rax
  20419d:	mov    $0x0,%edi
  2041a2:	cmovb  %r11,%rdi
  2041a6:	add    %rdi,%rax
  2041a9:	jmp    2041c3 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xfd3>
  2041ab:	nopl   0x0(%rax,%rax,1)
  2041b0:	sub    %rax,%rdx
  2041b3:	mov    $0x0,%r8d
  2041b9:	cmovb  %r11,%r8
  2041bd:	add    %r8,%rdx
  2041c0:	mov    %rdi,%r15
  2041c3:	cmp    $0x1,%r15
  2041c7:	sete   %dil
  2041cb:	cmp    $0x1,%rsi
  2041cf:	je     2041e7 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xff7>
  2041d1:	cmp    $0x1,%r15
  2041d5:	jne    204110 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf20>
  2041db:	jmp    2041e7 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xff7>
  2041dd:	xor    %eax,%eax
  2041df:	mov    0xb0(%rsp),%rbx
  2041e7:	test   %r11,%r11
  2041ea:	je     2046c5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14d5>
  2041f0:	test   %dil,%dil
  2041f3:	cmovne %rdx,%rax
  2041f7:	mul    %rcx
  2041fa:	mov    %rax,%rdi
  2041fd:	mov    %rdx,%rsi
  204200:	mov    %r11,%rdx
  204203:	xor    %ecx,%ecx
  204205:	call   *0x1c5535(%rip)        # 3c9740 <_DYNAMIC+0xa00>
  20420b:	mov    %rax,%r15
  20420e:	movdqa 0xf0(%rsp),%xmm0
  204217:	movdqu %xmm0,0x148(%rsp)
  204220:	movq   $0x1,0x140(%rsp)
  20422c:	movq   $0x0,0x180(%rsp)
  204238:	or     %rax,%rdx
  20423b:	mov    %rax,%r12
  20423e:	jne    2042b3 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10c3>
  204240:	jmp    204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  204245:	data16 cs nopw 0x0(%rax,%rax,1)
  204250:	mov    0x150(%rsp),%rax
  204258:	mov    %rax,0xd0(%rsp)
  204260:	movaps 0x140(%rsp),%xmm0
  204268:	movaps %xmm0,0xc0(%rsp)
  204270:	lea    0x248(%rsp),%rdi
  204278:	mov    %rbx,%rsi
  20427b:	mov    %rbp,%rdx
  20427e:	lea    0x140(%rsp),%rcx
  204286:	call   *0x1c5374(%rip)        # 3c9600 <_DYNAMIC+0x8c0>
  20428c:	movdqu 0x248(%rsp),%xmm0
  204295:	movdqa %xmm0,0x140(%rsp)
  20429e:	mov    0x258(%rsp),%rax
  2042a6:	mov    %rax,0x150(%rsp)
  2042ae:	shr    $1,%r12
  2042b1:	je     2042fc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x110c>
  2042b3:	test   $0x1,%r12b
  2042b7:	je     204250 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1060>
  2042b9:	mov    %rbp,%rdi
  2042bc:	mov    %rbx,%rsi
  2042bf:	lea    0x180(%rsp),%rdx
  2042c7:	lea    0x140(%rsp),%rcx
  2042cf:	call   *0x1c532b(%rip)        # 3c9600 <_DYNAMIC+0x8c0>
  2042d5:	mov    0xd0(%rsp),%rax
  2042dd:	mov    %rax,0x190(%rsp)
  2042e5:	movdqu 0xc0(%rsp),%xmm0
  2042ee:	movdqa %xmm0,0x180(%rsp)
  2042f7:	jmp    204250 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1060>
  2042fc:	testb  $0x1,0x180(%rsp)
  204304:	je     204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  204306:	mov    0x90(%rsp),%rax
  20430e:	cmp    %rax,0x188(%rsp)
  204316:	jne    204330 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1140>
  204318:	mov    0xa8(%rsp),%rax
  204320:	cmp    %rax,0x190(%rsp)
  204328:	je     20443e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x124e>
  20432e:	xchg   %ax,%ax
  204330:	cmpb   $0x0,0x50(%rsp)
  204335:	mov    0x18(%rsp),%rbp
  20433a:	mov    0x8(%rsp),%rbx
  20433f:	mov    0x88(%rsp),%r12
  204347:	jne    204362 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1172>
  204349:	mov    0x1c5358(%rip),%rax        # 3c96a8 <_DYNAMIC+0x968>
  204350:	mov    (%rax),%rax
  204353:	movabs $0x7fffffffffffffff,%rcx
  20435d:	test   %rcx,%rax
  204360:	jne    2043d6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11e6>
  204362:	xor    %eax,%eax
  204364:	mov    0x3f8(%rsp),%rcx
  20436c:	xchg   %eax,(%rcx)
  20436e:	cmp    $0x2,%eax
  204371:	je     2043c6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11d6>
  204373:	mov    0x48(%rsp),%rax
  204378:	movl   $0x0,(%rax,%r13,4)
  204380:	inc    %r13
  204383:	cmp    %rbx,%r13
  204386:	jb     20396c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x77c>
  20438c:	jmp    204405 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1215>
  20438e:	mov    0x3f8(%rsp),%rdi
  204396:	call   *0x1c5554(%rip)        # 3c98f0 <_DYNAMIC+0xbb0>
  20439c:	jmp    203e59 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xc69>
  2043a1:	call   *0x1c5311(%rip)        # 3c96b8 <_DYNAMIC+0x978>
  2043a7:	xor    $0x1,%al
  2043a9:	mov    %eax,0x50(%rsp)
  2043ad:	mov    0x3f8(%rsp),%rcx
  2043b5:	movzbl 0x4(%rcx),%eax
  2043b9:	test   %al,%al
  2043bb:	je     203e92 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xca2>
  2043c1:	jmp    2044f9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1309>
  2043c6:	mov    0x3f8(%rsp),%rdi
  2043ce:	call   *0x1c52dc(%rip)        # 3c96b0 <_DYNAMIC+0x970>
  2043d4:	jmp    204373 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1183>
  2043d6:	call   *0x1c52dc(%rip)        # 3c96b8 <_DYNAMIC+0x978>
  2043dc:	test   %al,%al
  2043de:	jne    204362 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1172>
  2043e0:	mov    0x3f8(%rsp),%rax
  2043e8:	movb   $0x1,0x4(%rax)
  2043ec:	jmp    204362 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1172>
  2043f1:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  204400:	mov    0x20(%rsp),%r14
  204405:	mov    0x80(%rsp),%rax
  20440d:	add    %rbx,%rax
  204410:	cmp    $0xffff,%rax
  204416:	ja     203860 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x670>
  20441c:	mov    %rax,0x80(%rsp)
  204424:	mov    0x3f8(%rsp),%r13
  20442c:	movzbl 0x60(%r13),%eax
  204431:	test   %al,%al
  204433:	je     203899 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x6a9>
  204439:	jmp    2036bc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4cc>
  20443e:	mov    0x3f8(%rsp),%r13
  204446:	lea    0x38(%r13),%rsi
  20444a:	lea    0x110(%rsp),%rdi
  204452:	mov    0x50(%rsp),%ebp
  204456:	call   25b760 <<std::sync::poison::mutex::Mutex<core::option::Option<alloc::vec::Vec<num_bigint::biguint::BigUint>>>>::lock>
  20445b:	movabs $0x7fffffffffffffff,%rbx
  204465:	cmpl   $0x1,0x110(%rsp)
  20446d:	je     204598 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13a8>
  204473:	mov    0x118(%rsp),%r12
  20447b:	cmpb   $0x0,0x120(%rsp)
  204483:	movq   $0x1,0x8(%r12)
  20448c:	mov    %r15,0x10(%r12)
  204491:	jne    2044a6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12b6>
  204493:	mov    0x1c520e(%rip),%rax        # 3c96a8 <_DYNAMIC+0x968>
  20449a:	mov    (%rax),%rax
  20449d:	test   %rbx,%rax
  2044a0:	jne    20461a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x142a>
  2044a6:	xor    %eax,%eax
  2044a8:	xchg   %eax,(%r12)
  2044ac:	cmp    $0x2,%eax
  2044af:	je     2045e4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13f4>
  2044b5:	movb   $0x1,0x60(%r13)
  2044ba:	test   %bpl,%bpl
  2044bd:	jne    2044d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12e2>
  2044bf:	mov    0x1c51e2(%rip),%rax        # 3c96a8 <_DYNAMIC+0x968>
  2044c6:	mov    (%rax),%rax
  2044c9:	test   %rbx,%rax
  2044cc:	jne    204633 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1443>
  2044d2:	xor    %eax,%eax
  2044d4:	xchg   %eax,0x0(%r13)
  2044d8:	cmp    $0x2,%eax
  2044db:	mov    0x18(%rsp),%rbp
  2044e0:	mov    0x8(%rsp),%rbx
  2044e5:	jne    2036bc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4cc>
  2044eb:	mov    %r13,%rdi
  2044ee:	call   *0x1c51bc(%rip)        # 3c96b0 <_DYNAMIC+0x970>
  2044f4:	jmp    2036bc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4cc>
  2044f9:	mov    %rcx,0x110(%rsp)
  204501:	mov    0x50(%rsp),%eax
  204505:	mov    %al,0x118(%rsp)
  20450c:	lea    -0x1e9028(%rip),%rdi        # 1b4eb <anon.9bc0095d74cae2c1d1f2191eec18534c.121.llvm.14577323480697803226+0x52>
  204513:	lea    0x1bd52e(%rip),%rcx        # 3c1a48 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x138>
  20451a:	lea    0x1bd5ef(%rip),%r8        # 3c1b10 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0x78>
  204521:	lea    0x110(%rsp),%rdx
  204529:	mov    $0x2b,%esi
  20452e:	call   *0x1c4b44(%rip)        # 3c9078 <_DYNAMIC+0x338>
  204534:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204539:	lea    0x1bd600(%rip),%rdx        # 3c1b40 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0xa8>
  204540:	jmp    204582 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1392>
  204542:	mov    %rbx,0xe0(%rsp)
  20454a:	lea    0x1bd577(%rip),%rdx        # 3c1ac8 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0x30>
  204551:	jmp    20457d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x138d>
  204553:	mov    0x8(%rsp),%rbx
  204558:	mov    %rbx,0xe0(%rsp)
  204560:	lea    0x1bd591(%rip),%rdx        # 3c1af8 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0x60>
  204567:	mov    0x18(%rsp),%rbp
  20456c:	jmp    20457d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x138d>
  20456e:	mov    %r10,%r13
  204571:	lea    0x1bd568(%rip),%rdx        # 3c1ae0 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0x48>
  204578:	mov    0x8(%rsp),%rbx
  20457d:	mov    0x20(%rsp),%r14
  204582:	mov    %r13,%rdi
  204585:	mov    0xe0(%rsp),%rsi
  20458d:	call   *0x1c4c9d(%rip)        # 3c9230 <_DYNAMIC+0x4f0>
  204593:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204598:	mov    0x118(%rsp),%rax
  2045a0:	movzbl 0x120(%rsp),%ecx
  2045a8:	mov    %rax,0xc0(%rsp)
  2045b0:	mov    %cl,0xc8(%rsp)
  2045b7:	lea    -0x1e90d3(%rip),%rdi        # 1b4eb <anon.9bc0095d74cae2c1d1f2191eec18534c.121.llvm.14577323480697803226+0x52>
  2045be:	lea    0x1bd463(%rip),%rcx        # 3c1a28 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x118>
  2045c5:	lea    0x1bd55c(%rip),%r8        # 3c1b28 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0x90>
  2045cc:	lea    0xc0(%rsp),%rdx
  2045d4:	mov    $0x2b,%esi
  2045d9:	call   *0x1c4a99(%rip)        # 3c9078 <_DYNAMIC+0x338>
  2045df:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  2045e4:	mov    %r12,%rdi
  2045e7:	call   *0x1c50c3(%rip)        # 3c96b0 <_DYNAMIC+0x970>
  2045ed:	jmp    2044b5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12c5>
  2045f2:	mov    $0x4,%edi
  2045f7:	mov    %r13,%rsi
  2045fa:	call   *0x1c49d0(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  204600:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204605:	mov    $0x1,%edi
  20460a:	mov    0x8(%rsp),%rsi
  20460f:	call   *0x1c49bb(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  204615:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  20461a:	call   *0x1c5098(%rip)        # 3c96b8 <_DYNAMIC+0x978>
  204620:	test   %al,%al
  204622:	jne    2044a6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12b6>
  204628:	movb   $0x1,0x4(%r12)
  20462e:	jmp    2044a6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12b6>
  204633:	mov    0x18(%rsp),%rbp
  204638:	mov    0x8(%rsp),%rbx
  20463d:	call   *0x1c5075(%rip)        # 3c96b8 <_DYNAMIC+0x978>
  204643:	test   %al,%al
  204645:	jne    2044d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12e2>
  20464b:	movb   $0x1,0x4(%r13)
  204650:	jmp    2044d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x12e2>
  204655:	mov    $0x8,%edi
  20465a:	mov    %r15,%rsi
  20465d:	call   *0x1c496d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  204663:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204665:	mov    $0x8,%edi
  20466a:	mov    %r15,%rsi
  20466d:	call   *0x1c495d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  204673:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204675:	mov    $0x8,%edi
  20467a:	mov    %r15,%rsi
  20467d:	call   *0x1c494d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  204683:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204685:	mov    $0x8,%edi
  20468a:	mov    %r15,%rsi
  20468d:	call   *0x1c493d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  204693:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  204695:	mov    $0x8,%edi
  20469a:	mov    %r15,%rsi
  20469d:	call   *0x1c492d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  2046a3:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  2046a5:	mov    $0x8,%edi
  2046aa:	mov    %r15,%rsi
  2046ad:	call   *0x1c491d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  2046b3:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  2046b5:	mov    $0x8,%edi
  2046ba:	mov    %r15,%rsi
  2046bd:	call   *0x1c490d(%rip)        # 3c8fd0 <_DYNAMIC+0x290>
  2046c3:	jmp    2046d2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e2>
  2046c5:	lea    0x1bd4d4(%rip),%rdi        # 3c1ba0 <anon.2ead17125980934146c792e9bfa8463b.125.llvm.15942493848686059347+0x108>
  2046cc:	call   *0x1c5046(%rip)        # 3c9718 <_DYNAMIC+0x9d8>
  2046d2:	ud2
  2046d4:	jmp    2047e4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15f4>
  2046d9:	jmp    2047e4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15f4>
  2046de:	jmp    2047ad <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15bd>
  2046e3:	mov    %rax,%r15
  2046e6:	mov    0x18(%rsp),%rbp
  2046eb:	mov    0x8(%rsp),%rbx
  2046f0:	jmp    20485e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x166e>
  2046f5:	mov    %rax,%r15
  2046f8:	mov    0x8(%rsp),%rbx
  2046fd:	mov    %rbx,%r14
  204700:	mov    0x18(%rsp),%rbp
  204705:	jmp    204899 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16a9>
  20470a:	mov    %rax,%r15
  20470d:	mov    0x18(%rsp),%rbp
  204712:	mov    0x8(%rsp),%rbx
  204717:	jmp    2048b5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16c5>
  20471c:	mov    %rax,%r15
  20471f:	mov    0x8(%rsp),%rbx
  204724:	test   %rbx,%rbx
  204727:	jne    2048f1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1701>
  20472d:	jmp    20494e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x175e>
  204732:	mov    %rax,%r15
  204735:	mov    0x8(%rsp),%rbx
  20473a:	test   %rbx,%rbx
  20473d:	jne    204909 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1719>
  204743:	jmp    20494e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x175e>
  204748:	mov    %rax,%r15
  20474b:	mov    0x8(%rsp),%rbx
  204750:	test   %rbx,%rbx
  204753:	jne    204921 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1731>
  204759:	jmp    20494e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x175e>
  20475e:	mov    %rax,%r15
  204761:	mov    0x8(%rsp),%rbx
  204766:	test   %rbx,%rbx
  204769:	jne    204939 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1749>
  20476f:	jmp    20494e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x175e>
  204774:	jmp    2047e4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15f4>
  204776:	mov    %rax,%r15
  204779:	mov    0x18(%rsp),%rbp
  20477e:	mov    0x8(%rsp),%rbx
  204783:	jmp    204874 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1684>
  204788:	mov    %rax,%r15
  20478b:	mov    0x18(%rsp),%rbp
  204790:	jmp    2048cc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16dc>
  204795:	mov    %rax,%r15
  204798:	lea    0xc0(%rsp),%rdi
  2047a0:	call   200210 <core::ptr::drop_glue::<std::sync::poison::PoisonError<std::sync::poison::mutex::MutexGuard<core::option::Option<u64>>>>>
  2047a5:	jmp    2047b0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15c0>
  2047a7:	call   *0x1c480b(%rip)        # 3c8fb8 <_DYNAMIC+0x278>
  2047ad:	mov    %rax,%r15
  2047b0:	movzbl 0x50(%rsp),%esi
  2047b5:	mov    0x3f8(%rsp),%rdi
  2047bd:	call   2002f0 <core::ptr::drop_glue::<std::sync::poison::mutex::MutexGuard<std::collections::hash::map::HashMap<u64, (u64, u64, u64)>>>>
  2047c2:	xor    %r14d,%r14d
  2047c5:	mov    0x18(%rsp),%rbp
  2047ca:	mov    0x8(%rsp),%rbx
  2047cf:	jmp    204835 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1645>
  2047d1:	call   *0x1c47e1(%rip)        # 3c8fb8 <_DYNAMIC+0x278>
  2047d7:	jmp    2047ec <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15fc>
  2047d9:	mov    %rax,%r15
  2047dc:	mov    %r14,%rbx
  2047df:	xor    %r14d,%r14d
  2047e2:	jmp    204835 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1645>
  2047e4:	mov    %rax,%r15
  2047e7:	xor    %r14d,%r14d
  2047ea:	jmp    204835 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1645>
  2047ec:	mov    %rax,%r15
  2047ef:	jmp    204835 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1645>
  2047f1:	mov    %rax,%r15
  2047f4:	mov    0x110(%rsp),%r12
  2047fc:	cmpb   $0x0,0x118(%rsp)
  204804:	jne    204823 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1633>
  204806:	mov    0x1c4e9b(%rip),%rax        # 3c96a8 <_DYNAMIC+0x968>
  20480d:	mov    (%rax),%rax
  204810:	movabs $0x7fffffffffffffff,%rcx
  20481a:	test   %rcx,%rax
  20481d:	jne    204964 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1774>
  204823:	xor    %r14d,%r14d
  204826:	xor    %eax,%eax
  204828:	xchg   %eax,(%r12)
  20482c:	cmp    $0x2,%eax
  20482f:	je     204956 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1766>
  204835:	mov    0x158(%rsp),%rsi
  20483d:	test   %rsi,%rsi
  204840:	je     204859 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1669>
  204842:	mov    0x160(%rsp),%rdi
  20484a:	shl    $0x3,%rsi
  20484e:	mov    $0x8,%edx
  204853:	call   *0x1c4757(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  204859:	test   %r14b,%r14b
  20485c:	jne    204874 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1684>
  20485e:	mov    $0x1,%edx
  204863:	mov    0x98(%rsp),%rdi
  20486b:	mov    %rbx,%rsi
  20486e:	call   *0x1c473c(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  204874:	mov    0xa0(%rsp),%r14
  20487c:	test   %r14,%r14
  20487f:	je     2048b0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16c0>
  204881:	lea    0x0(,%r14,8),%rsi
  204889:	mov    $0x8,%edx
  20488e:	mov    0x78(%rsp),%rdi
  204893:	call   *0x1c4717(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  204899:	shl    $0x3,%r14
  20489d:	mov    $0x8,%edx
  2048a2:	mov    0x58(%rsp),%rdi
  2048a7:	mov    %r14,%rsi
  2048aa:	call   *0x1c4700(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  2048b0:	test   %rbx,%rbx
  2048b3:	je     2048cc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16dc>
  2048b5:	shl    $0x2,%rbx
  2048b9:	mov    $0x4,%edx
  2048be:	mov    0x48(%rsp),%rdi
  2048c3:	mov    %rbx,%rsi
  2048c6:	call   *0x1c46e4(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  2048cc:	mov    0x108(%rsp),%rbx
  2048d4:	test   %rbx,%rbx
  2048d7:	je     20494e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x175e>
  2048d9:	lea    0x0(,%rbx,8),%rsi
  2048e1:	mov    $0x8,%edx
  2048e6:	mov    0x28(%rsp),%rdi
  2048eb:	call   *0x1c46bf(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  2048f1:	lea    0x0(,%rbx,8),%rsi
  2048f9:	mov    $0x8,%edx
  2048fe:	mov    0x38(%rsp),%rdi
  204903:	call   *0x1c46a7(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  204909:	lea    0x0(,%rbx,8),%rsi
  204911:	mov    $0x8,%edx
  204916:	mov    0x30(%rsp),%rdi
  20491b:	call   *0x1c468f(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  204921:	lea    0x0(,%rbx,8),%rsi
  204929:	mov    $0x8,%edx
  20492e:	mov    0x40(%rsp),%rdi
  204933:	call   *0x1c4677(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  204939:	shl    $0x3,%rbx
  20493d:	mov    $0x8,%edx
  204942:	mov    %rbp,%rdi
  204945:	mov    %rbx,%rsi
  204948:	call   *0x1c4662(%rip)        # 3c8fb0 <_DYNAMIC+0x270>
  20494e:	mov    %r15,%rdi
  204951:	call   3bec00 <_Unwind_Resume@plt>
  204956:	mov    %r12,%rdi
  204959:	call   *0x1c4d51(%rip)        # 3c96b0 <_DYNAMIC+0x970>
  20495f:	jmp    2047c2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x15d2>
  204964:	call   *0x1c4d4e(%rip)        # 3c96b8 <_DYNAMIC+0x978>
  20496a:	test   %al,%al
  20496c:	mov    0x18(%rsp),%rbp
  204971:	mov    0x8(%rsp),%rbx
  204976:	jne    204823 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1633>
  20497c:	movb   $0x1,0x4(%r12)
  204982:	jmp    204823 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1633>
  204987:	call   *0x1c462b(%rip)        # 3c8fb8 <_DYNAMIC+0x278>
