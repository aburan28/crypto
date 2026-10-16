000000000022cf90 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf>:
  22cf90:	push   %rbp
  22cf91:	push   %r15
  22cf93:	push   %r14
  22cf95:	push   %r13
  22cf97:	push   %r12
  22cf99:	push   %rbx
  22cf9a:	sub    $0x38,%rsp
  22cf9e:	mov    %rcx,0x30(%rsp)
  22cfa3:	mov    %rdx,%r15
  22cfa6:	mov    %rsi,%rbp
  22cfa9:	mov    %rdi,%r12
  22cfac:	mov    0x10(%rsi),%r13
  22cfb0:	test   %r13,%r13
  22cfb3:	je     22d000 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x70>
  22cfb5:	mov    0x40(%rbp),%rdi
  22cfb9:	mov    0x8(%rbp),%rax
  22cfbd:	mov    0x38(%rbp),%rcx
  22cfc1:	mov    0x120(%r12),%r8
  22cfc9:	lea    -0x1(%r13),%rdx
  22cfcd:	cmp    %rdx,%rdi
  22cfd0:	cmovb  %rdi,%rdx
  22cfd4:	inc    %rdx
  22cfd7:	cmp    $0x11,%rdx
  22cfdb:	setb   %sil
  22cfdf:	mov    %rcx,%r9
  22cfe2:	sub    %rax,%r9
  22cfe5:	cmp    $0x80,%r9
  22cfec:	setb   %r9b
  22cff0:	or     %sil,%r9b
  22cff3:	cmp    $0x1,%r9b
  22cff7:	jne    22d00d <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x7d>
  22cff9:	xor    %edx,%edx
  22cffb:	jmp    22d120 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x190>
  22d000:	mov    0x100(%r12),%rsi
  22d008:	jmp    22d1e1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x251>
  22d00d:	mov    %edx,%esi
  22d00f:	and    $0xf,%esi
  22d012:	mov    $0x10,%r9d
  22d018:	cmovne %rsi,%r9
  22d01c:	sub    %r9,%rdx
  22d01f:	vpbroadcastq %r8,%ymm0
  22d025:	xor    %esi,%esi
  22d027:	vpbroadcastq -0x21b680(%rip),%ymm1        # 119b0 <anon.b68ad25efc32a4f1d321fe93ba3d9fb4.405.llvm.17162556477609682337+0x10>
  22d030:	vmovdqu (%rax,%rsi,8),%ymm2
  22d035:	vpsrlq $0x1,%ymm2,%ymm3
  22d03a:	vpand  %ymm1,%ymm3,%ymm3
  22d03e:	kxnorw %k0,%k0,%k1
  22d042:	vpxor  %xmm4,%xmm4,%xmm4
  22d046:	vpgatherqq (%r15,%ymm3,8),%ymm4{%k1}
  22d04d:	vmovdqu 0x20(%rax,%rsi,8),%ymm3
  22d053:	vpsrlq $0x1,%ymm3,%ymm5
  22d058:	vpand  %ymm1,%ymm5,%ymm5
  22d05c:	kxnorw %k0,%k0,%k1
  22d060:	vpxor  %xmm6,%xmm6,%xmm6
  22d064:	vpgatherqq (%r15,%ymm5,8),%ymm6{%k1}
  22d06b:	vmovdqu 0x40(%rax,%rsi,8),%ymm5
  22d071:	vpsrlq $0x1,%ymm5,%ymm7
  22d076:	vpand  %ymm1,%ymm7,%ymm7
  22d07a:	kxnorw %k0,%k0,%k1
  22d07e:	vpxor  %xmm8,%xmm8,%xmm8
  22d083:	vpgatherqq (%r15,%ymm7,8),%ymm8{%k1}
  22d08a:	vmovdqu 0x60(%rax,%rsi,8),%ymm7
  22d090:	vpsrlq $0x1,%ymm7,%ymm9
  22d095:	vpand  %ymm1,%ymm9,%ymm9
  22d099:	kxnorw %k0,%k0,%k1
  22d09d:	vpxor  %xmm10,%xmm10,%xmm10
  22d0a2:	vpgatherqq (%r15,%ymm9,8),%ymm10{%k1}
  22d0a9:	vpsubq %ymm2,%ymm4,%ymm9
  22d0ad:	vpsubq %ymm3,%ymm6,%ymm11
  22d0b1:	vpsubq %ymm5,%ymm8,%ymm12
  22d0b5:	vpsubq %ymm7,%ymm10,%ymm13
  22d0b9:	vpcmpltuq %ymm2,%ymm4,%k1
  22d0c0:	vpcmpltuq %ymm3,%ymm6,%k2
  22d0c7:	vpcmpltuq %ymm5,%ymm8,%k3
  22d0ce:	vpcmpltuq %ymm7,%ymm10,%k4
  22d0d5:	vpaddq %ymm0,%ymm9,%ymm9{%k1}
  22d0db:	vpaddq %ymm0,%ymm11,%ymm11{%k2}
  22d0e1:	vpaddq %ymm0,%ymm12,%ymm12{%k3}
  22d0e7:	vpaddq %ymm0,%ymm13,%ymm13{%k4}
  22d0ed:	vmovdqu %ymm9,(%rcx,%rsi,8)
  22d0f2:	vmovdqu %ymm11,0x20(%rcx,%rsi,8)
  22d0f8:	vmovdqu %ymm12,0x40(%rcx,%rsi,8)
  22d0fe:	vmovdqu %ymm13,0x60(%rcx,%rsi,8)
  22d104:	add    $0x10,%rsi
  22d108:	cmp    %rsi,%rdx
  22d10b:	jne    22d030 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0xa0>
  22d111:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  22d120:	cmp    %rdx,%rdi
  22d123:	je     22d469 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x4d9>
  22d129:	mov    (%rax,%rdx,8),%rsi
  22d12d:	mov    %esi,%r9d
  22d130:	and    $0x3e,%r9d
  22d134:	mov    (%r15,%r9,4),%r9
  22d138:	sub    %rsi,%r9
  22d13b:	mov    $0x0,%esi
  22d140:	cmovb  %r8,%rsi
  22d144:	add    %r9,%rsi
  22d147:	mov    %rsi,(%rcx,%rdx,8)
  22d14b:	inc    %rdx
  22d14e:	cmp    %rdx,%r13
  22d151:	jne    22d120 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x190>
  22d153:	mov    0x100(%r12),%rsi
  22d15b:	mov    0x128(%r12),%r9
  22d163:	mov    0x50(%rbp),%r10
  22d167:	mov    0x58(%rbp),%rax
  22d16b:	mov    0x60(%rbp),%r11
  22d16f:	xor    %ebx,%ebx
  22d171:	jmp    22d1b3 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x223>
  22d173:	data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  22d180:	mulx   %rsi,%r14,%rsi
  22d185:	mov    %r9,%rdx
  22d188:	imul   %r14,%rdx
  22d18c:	mulx   %r8,%rdx,%rdx
  22d191:	cmp    $0x1,%r14
  22d195:	sbb    $0xffffffffffffffff,%rsi
  22d199:	add    %rdx,%rsi
  22d19c:	cmp    %r8,%rsi
  22d19f:	mov    $0x0,%edx
  22d1a4:	cmovae %r8,%rdx
  22d1a8:	sub    %rdx,%rsi
  22d1ab:	inc    %rbx
  22d1ae:	cmp    %rbx,%r13
  22d1b1:	je     22d1e1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x251>
  22d1b3:	cmp    %rbx,%rax
  22d1b6:	je     22d47c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x4ec>
  22d1bc:	mov    %rsi,(%r10,%rbx,8)
  22d1c0:	cmp    %rbx,%rdi
  22d1c3:	je     22d492 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x502>
  22d1c9:	mov    (%rcx,%rbx,8),%rdx
  22d1cd:	test   %rdx,%rdx
  22d1d0:	jne    22d180 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x1f0>
  22d1d2:	inc    %r11
  22d1d5:	mov    %r11,0x60(%rbp)
  22d1d9:	inc    %rbx
  22d1dc:	cmp    %rbx,%r13
  22d1df:	jne    22d1b3 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x223>
  22d1e1:	mov    %r12,%rdi
  22d1e4:	vzeroupper
  22d1e7:	call   *0x187f4b(%rip)        # 3b5138 <_DYNAMIC+0xcd0>
  22d1ed:	test   $0x1,%al
  22d1ef:	je     22d456 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x4c6>
  22d1f5:	mov    0x38(%rbp),%rax
  22d1f9:	mov    %rax,0x10(%rsp)
  22d1fe:	mov    0x40(%rbp),%rcx
  22d202:	mov    0x120(%r12),%r10
  22d20a:	mov    0x128(%r12),%r11
  22d212:	lea    -0x1(%r13),%rax
  22d216:	mov    %r13,%rdi
  22d219:	mov    %rcx,0x8(%rsp)
  22d21e:	cmp    %rcx,%rax
  22d221:	jae    22d2d8 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x348>
  22d227:	mov    0x58(%rbp),%rax
  22d22b:	mov    0x50(%rbp),%rbx
  22d22f:	xor    %r14d,%r14d
  22d232:	mov    %r13,%rdi
  22d235:	data16 cs nopw 0x0(%rax,%rax,1)
  22d240:	test   %rdi,%rdi
  22d243:	je     22d2e1 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x351>
  22d249:	mov    0x10(%rsp),%rcx
  22d24e:	mov    -0x8(%rcx,%rdi,8),%rcx
  22d253:	dec    %rdi
  22d256:	test   %rcx,%rcx
  22d259:	je     22d240 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x2b0>
  22d25b:	cmp    %rax,%rdi
  22d25e:	jae    22d4e4 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x554>
  22d264:	mov    %rdx,%r8
  22d267:	mulx   (%rbx,%rdi,8),%rsi,%r9
  22d26d:	mov    %r11,%rdx
  22d270:	imul   %rsi,%rdx
  22d274:	mulx   %r10,%rdx,%rdx
  22d279:	cmp    $0x1,%rsi
  22d27d:	sbb    $0xffffffffffffffff,%r9
  22d281:	add    %rdx,%r9
  22d284:	cmp    %r10,%r9
  22d287:	mov    %r10,%rdx
  22d28a:	cmovb  %r14,%rdx
  22d28e:	sub    %rdx,%r9
  22d291:	mov    0x10(%rsp),%rdx
  22d296:	mov    %r9,(%rdx,%rdi,8)
  22d29a:	mov    %rcx,%rdx
  22d29d:	mulx   %r8,%rcx,%r8
  22d2a2:	mov    %r11,%rdx
  22d2a5:	imul   %rcx,%rdx
  22d2a9:	mulx   %r10,%rsi,%rsi
  22d2ae:	mov    %r8,%rdx
  22d2b1:	cmp    $0x1,%rcx
  22d2b5:	sbb    $0xffffffffffffffff,%rdx
  22d2b9:	add    %rsi,%rdx
  22d2bc:	cmp    %r10,%rdx
  22d2bf:	mov    %r10,%rcx
  22d2c2:	cmovb  %r14,%rcx
  22d2c6:	sub    %rcx,%rdx
  22d2c9:	lea    -0x1(%rdi),%rcx
  22d2cd:	cmp    0x8(%rsp),%rcx
  22d2d2:	jb     22d240 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x2b0>
  22d2d8:	test   %rdi,%rdi
  22d2db:	jne    22d4cf <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x53f>
  22d2e1:	mov    %rbp,0x18(%rsp)
  22d2e6:	test   %r13,%r13
  22d2e9:	je     22d45a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x4ca>
  22d2ef:	mov    0x18(%rsp),%rax
  22d2f4:	mov    0x8(%rax),%rcx
  22d2f8:	mov    0x20(%rax),%r8
  22d2fc:	mov    0x28(%rax),%rdx
  22d300:	mov    %rdx,0x28(%rsp)
  22d305:	mov    0x68(%rax),%rax
  22d309:	mov    %rax,0x20(%rsp)
  22d30e:	xor    %ebp,%ebp
  22d310:	xor    %ebx,%ebx
  22d312:	jmp    22d32c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x39c>
  22d314:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  22d320:	inc    %rbx
  22d323:	cmp    %rbx,%r13
  22d326:	je     22d45a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x4ca>
  22d32c:	cmp    %rbx,0x28(%rsp)
  22d331:	je     22d4a5 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x515>
  22d337:	cmp    %rbx,0x8(%rsp)
  22d33c:	je     22d4ba <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x52a>
  22d342:	mov    (%rcx,%rbx,8),%r9
  22d346:	mov    %r9d,%esi
  22d349:	and    $0x3e,%esi
  22d34c:	mov    0x100(%r15,%rsi,4),%rax
  22d354:	mov    (%r8,%rbx,8),%r14
  22d358:	sub    %r14,%rax
  22d35b:	mov    $0x0,%edx
  22d360:	cmovb  %r10,%rdx
  22d364:	add    %rax,%rdx
  22d367:	mov    0x10(%rsp),%rax
  22d36c:	mulx   (%rax,%rbx,8),%r12,%rax
  22d372:	mov    %r11,%rdx
  22d375:	imul   %r12,%rdx
  22d379:	mulx   %r10,%rdx,%rdx
  22d37e:	cmp    $0x1,%r12
  22d382:	sbb    $0xffffffffffffffff,%rax
  22d386:	add    %rdx,%rax
  22d389:	cmp    %r10,%rax
  22d38c:	mov    %r10,%rdx
  22d38f:	cmovb  %rbp,%rdx
  22d393:	sub    %rdx,%rax
  22d396:	mov    %rax,%rdx
  22d399:	mulx   %rax,%r12,%rdi
  22d39e:	mov    %r11,%rdx
  22d3a1:	imul   %r12,%rdx
  22d3a5:	cmp    $0x1,%r12
  22d3a9:	mulx   %r10,%rdx,%rdx
  22d3ae:	sbb    $0xffffffffffffffff,%rdi
  22d3b2:	add    %rdx,%rdi
  22d3b5:	cmp    %r10,%rdi
  22d3b8:	mov    %r10,%rdx
  22d3bb:	cmovb  %rbp,%rdx
  22d3bf:	sub    %rdx,%rdi
  22d3c2:	sub    %r9,%rdi
  22d3c5:	mov    $0x0,%edx
  22d3ca:	cmovb  %r10,%rdx
  22d3ce:	add    %rdi,%rdx
  22d3d1:	sub    (%r15,%rsi,4),%rdx
  22d3d5:	mov    $0x0,%esi
  22d3da:	cmovb  %r10,%rsi
  22d3de:	add    %rdx,%rsi
  22d3e1:	sub    %rsi,%r9
  22d3e4:	mov    $0x0,%edx
  22d3e9:	cmovb  %r10,%rdx
  22d3ed:	add    %r9,%rdx
  22d3f0:	mulx   %rax,%rax,%rdi
  22d3f5:	mov    %r11,%rdx
  22d3f8:	imul   %rax,%rdx
  22d3fc:	cmp    $0x1,%rax
  22d400:	mulx   %r10,%rax,%rax
  22d405:	sbb    $0xffffffffffffffff,%rdi
  22d409:	add    %rax,%rdi
  22d40c:	cmp    %r10,%rdi
  22d40f:	mov    %r10,%rax
  22d412:	cmovb  %rbp,%rax
  22d416:	sub    %rax,%rdi
  22d419:	sub    %r14,%rdi
  22d41c:	mov    $0x0,%eax
  22d421:	cmovb  %r10,%rax
  22d425:	add    %rdi,%rax
  22d428:	mov    %rsi,(%rcx,%rbx,8)
  22d42c:	mov    %rax,(%r8,%rbx,8)
  22d430:	test   %rsi,0x30(%rsp)
  22d435:	jne    22d320 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x390>
  22d43b:	mov    0x20(%rsp),%rdx
  22d440:	inc    %rdx
  22d443:	mov    0x18(%rsp),%rax
  22d448:	mov    %rdx,0x20(%rsp)
  22d44d:	mov    %rdx,0x68(%rax)
  22d451:	jmp    22d320 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64bf+0x390>
  22d456:	incq   0x60(%rbp)
  22d45a:	add    $0x38,%rsp
  22d45e:	pop    %rbx
  22d45f:	pop    %r12
  22d461:	pop    %r13
  22d463:	pop    %r14
  22d465:	pop    %r15
  22d467:	pop    %rbp
  22d468:	ret
  22d469:	lea    0x180798(%rip),%rdx        # 3adc08 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x3a0>
  22d470:	mov    %rdi,%rsi
  22d473:	vzeroupper
  22d476:	call   *0x1874ec(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d47c:	lea    0x180755(%rip),%rdx        # 3adbd8 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x370>
  22d483:	mov    %rax,%rdi
  22d486:	mov    %rax,%rsi
  22d489:	vzeroupper
  22d48c:	call   *0x1874d6(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d492:	lea    0x180757(%rip),%rdx        # 3adbf0 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x388>
  22d499:	mov    %rdi,%rsi
  22d49c:	vzeroupper
  22d49f:	call   *0x1874c3(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d4a5:	lea    0x1806b4(%rip),%rdx        # 3adb60 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x2f8>
  22d4ac:	mov    0x28(%rsp),%rdi
  22d4b1:	mov    %rdi,%rsi
  22d4b4:	call   *0x1874ae(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d4ba:	lea    0x1806b7(%rip),%rdx        # 3adb78 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x310>
  22d4c1:	mov    0x8(%rsp),%rdi
  22d4c6:	mov    %rdi,%rsi
  22d4c9:	call   *0x187499(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d4cf:	dec    %rdi
  22d4d2:	lea    0x1806cf(%rip),%rdx        # 3adba8 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x340>
  22d4d9:	mov    0x8(%rsp),%rsi
  22d4de:	call   *0x187484(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d4e4:	lea    0x1806d5(%rip),%rdx        # 3adbc0 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x358>
  22d4eb:	mov    %rax,%rsi
  22d4ee:	call   *0x187474(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22d4f4:	int3
  22d4f5:	int3
  22d4f6:	int3
  22d4f7:	int3
  22d4f8:	int3
  22d4f9:	int3
  22d4fa:	int3
  22d4fb:	int3
  22d4fc:	int3
  22d4fd:	int3
  22d4fe:	int3
  22d4ff:	int3

