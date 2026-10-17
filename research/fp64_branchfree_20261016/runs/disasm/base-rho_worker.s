  1f65b0:	push   %rbp
  1f65b1:	push   %r15
  1f65b3:	push   %r14
  1f65b5:	push   %r13
  1f65b7:	push   %r12
  1f65b9:	push   %rbx
  1f65ba:	sub    $0x3a8,%rsp
  1f65c1:	mov    %r9,%r13
  1f65c4:	mov    %rdi,%rbp
  1f65c7:	mov    0x3e8(%rsp),%r15
  1f65cf:	movabs $0x9e3779b97f4a7c15,%rax
  1f65d9:	movabs $0x7fffffffffffffff,%r12
  1f65e3:	mov    %rsi,0x1b0(%rsp)
  1f65eb:	mov    %rdx,0x1b8(%rsp)
  1f65f3:	mov    %rcx,0x180(%rsp)
  1f65fb:	mov    %r8,0x188(%rsp)
  1f6603:	movdqu 0x18(%rdi),%xmm0
  1f6608:	movdqu 0x28(%rdi),%xmm1
  1f660d:	movdqu 0x38(%rdi),%xmm2
  1f6612:	movdqa %xmm2,0x230(%rsp)
  1f661b:	movdqa %xmm1,0x220(%rsp)
  1f6624:	movdqa %xmm0,0x210(%rsp)
  1f662d:	mov    0x58(%rdi),%rbx
  1f6631:	mov    %rbx,0xb8(%rsp)
  1f6639:	mov    0x8(%r15),%r14
  1f663d:	mov    0x3f0(%rsp),%rsi
  1f6645:	imul   %rax,%rsi
  1f6649:	add    %rax,%rsi
  1f664c:	xor    0x10(%r15),%rsi
  1f6650:	lea    0x260(%rsp),%rdi
  1f6658:	call   2d25e0 <<rand::rngs::std::StdRng as rand_core::SeedableRng>::seed_from_u64>
  1f665d:	mov    0x24(%r15),%ecx
  1f6661:	mov    $0x28,%eax
  1f6666:	shl    %cl,%rax
  1f6669:	cmp    $0x101,%rax
  1f666f:	mov    $0x100,%esi
  1f6674:	cmovae %rax,%rsi
  1f6678:	mov    $0xffffffff,%eax
  1f667d:	cmp    %rax,%rsi
  1f6680:	cmovae %rax,%rsi
  1f6684:	lea    0x0(,%r14,8),%r15
  1f668c:	mov    %r14,%rax
  1f668f:	shr    $0x3d,%rax
  1f6693:	setne  %al
  1f6696:	lea    -0x7(%r12),%rdx
  1f669b:	cmp    %rdx,%r15
  1f669e:	seta   %dl
  1f66a1:	or     %al,%dl
  1f66a3:	je     1f66b4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x104>
  1f66a5:	xor    %r13d,%r13d
  1f66a8:	mov    %r13,%rdi
  1f66ab:	mov    %r15,%rsi
  1f66ae:	call   *0x1d2e0c(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f66b4:	mov    %rsi,0x1d0(%rsp)
  1f66bc:	mov    %ecx,0x6c(%rsp)
  1f66c0:	mov    %rbx,0xe8(%rsp)
  1f66c8:	mov    %rbp,0xb0(%rsp)
  1f66d0:	mov    %r13,0x80(%rsp)
  1f66d8:	mov    %r14,0x8(%rsp)
  1f66dd:	test   %r15,%r15
  1f66e0:	je     1f67d8 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x228>
  1f66e6:	call   *0x1d2ddc(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f66ec:	mov    $0x8,%r13d
  1f66f2:	mov    $0x8,%esi
  1f66f7:	mov    %r15,%rdi
  1f66fa:	call   *0x1d31a8(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f6700:	test   %rax,%rax
  1f6703:	mov    0x80(%rsp),%r12
  1f670b:	mov    0x6c(%rsp),%ebx
  1f670f:	je     1f66a8 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf8>
  1f6711:	mov    %rax,%rbp
  1f6714:	call   *0x1d2dae(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f671a:	mov    $0x8,%esi
  1f671f:	mov    %r15,%rdi
  1f6722:	call   *0x1d3180(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f6728:	mov    %rax,0x40(%rsp)
  1f672d:	test   %rax,%rax
  1f6730:	je     1f7a79 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14c9>
  1f6736:	call   *0x1d2d8c(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f673c:	mov    $0x8,%esi
  1f6741:	mov    %r15,%rdi
  1f6744:	call   *0x1d315e(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f674a:	mov    %rax,0x30(%rsp)
  1f674f:	test   %rax,%rax
  1f6752:	je     1f7a89 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14d9>
  1f6758:	call   *0x1d2d6a(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f675e:	mov    $0x8,%esi
  1f6763:	mov    %r15,%rdi
  1f6766:	call   *0x1d313c(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f676c:	mov    %rax,0x38(%rsp)
  1f6771:	test   %rax,%rax
  1f6774:	je     1f7a99 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14e9>
  1f677a:	call   *0x1d2d48(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f6780:	mov    $0x8,%esi
  1f6785:	mov    %r15,%rdi
  1f6788:	call   *0x1d2d42(%rip)        # 3c94d0 <_DYNAMIC+0x2a0>
  1f678e:	mov    %rax,0x28(%rsp)
  1f6793:	test   %rax,%rax
  1f6796:	je     1f7aa9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14f9>
  1f679c:	mov    0x8(%rsp),%rcx
  1f67a1:	mov    %rcx,%rdx
  1f67a4:	cmp    $0x2,%rcx
  1f67a8:	mov    %rbp,0x18(%rsp)
  1f67ad:	mov    %rdx,0x108(%rsp)
  1f67b5:	jb     1f681c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x26c>
  1f67b7:	lea    -0x8(%r15),%rdx
  1f67bb:	mov    0x28(%rsp),%r14
  1f67c0:	mov    %r14,%rdi
  1f67c3:	mov    $0xff,%esi
  1f67c8:	call   *0x1d33aa(%rip)        # 3c9b78 <memset@GLIBC_2.2.5>
  1f67ce:	lea    (%r14,%r15,1),%rax
  1f67d2:	add    $0xfffffffffffffff8,%rax
  1f67d6:	jmp    1f6832 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x282>
  1f67d8:	mov    $0x8,%eax
  1f67dd:	mov    %rax,0x38(%rsp)
  1f67e2:	xor    %edx,%edx
  1f67e4:	mov    %rax,0x40(%rsp)
  1f67e9:	mov    $0x8,%ebp
  1f67ee:	mov    %rax,0x30(%rsp)
  1f67f3:	mov    %rax,0x28(%rsp)
  1f67f8:	mov    0x8(%rsp),%rcx
  1f67fd:	mov    0x80(%rsp),%r12
  1f6805:	mov    0x6c(%rsp),%ebx
  1f6809:	cmp    $0x2,%rcx
  1f680d:	mov    %rbp,0x18(%rsp)
  1f6812:	mov    %rdx,0x108(%rsp)
  1f681a:	jae    1f67b7 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x207>
  1f681c:	mov    0x28(%rsp),%rax
  1f6821:	test   %rcx,%rcx
  1f6824:	mov    0x3e8(%rsp),%r14
  1f682c:	je     1f6941 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x391>
  1f6832:	movq   $0xffffffffffffffff,(%rax)
  1f6839:	mov    0x8(%rsp),%rax
  1f683e:	lea    0x0(,%rax,4),%r13
  1f6846:	call   *0x1d2c7c(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f684c:	mov    $0x4,%esi
  1f6851:	mov    %r13,%rdi
  1f6854:	call   *0x1d304e(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f685a:	mov    %rax,0x48(%rsp)
  1f685f:	test   %rax,%rax
  1f6862:	je     1f7a0f <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x145f>
  1f6868:	test   %r15,%r15
  1f686b:	je     1f68c0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x310>
  1f686d:	call   *0x1d2c55(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f6873:	mov    $0x8,%esi
  1f6878:	mov    %r15,%rdi
  1f687b:	call   *0x1d3027(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f6881:	mov    %rax,0x60(%rsp)
  1f6886:	test   %rax,%rax
  1f6889:	je     1f7ab9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1509>
  1f688f:	call   *0x1d2c33(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f6895:	mov    $0x8,%esi
  1f689a:	mov    %r15,%rdi
  1f689d:	call   *0x1d3005(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f68a3:	mov    %rax,0x78(%rsp)
  1f68a8:	test   %rax,%rax
  1f68ab:	je     1f7ac9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1519>
  1f68b1:	mov    0x8(%rsp),%r14
  1f68b6:	mov    %r14,0xa8(%rsp)
  1f68be:	jmp    1f68e0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x330>
  1f68c0:	mov    $0x8,%eax
  1f68c5:	mov    %rax,0x78(%rsp)
  1f68ca:	movq   $0x0,0xa8(%rsp)
  1f68d6:	mov    %rax,0x60(%rsp)
  1f68db:	mov    0x8(%rsp),%r14
  1f68e0:	call   *0x1d2be2(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f68e6:	mov    $0x1,%esi
  1f68eb:	mov    %r14,%rdi
  1f68ee:	call   *0x1d2fb4(%rip)        # 3c98a8 <_DYNAMIC+0x678>
  1f68f4:	mov    %rax,0xa0(%rsp)
  1f68fc:	test   %rax,%rax
  1f68ff:	je     1f7a22 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1472>
  1f6905:	test   %r15,%r15
  1f6908:	mov    0x3e8(%rsp),%r14
  1f6910:	je     1f6938 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x388>
  1f6912:	call   *0x1d2bb0(%rip)        # 3c94c8 <_DYNAMIC+0x298>
  1f6918:	mov    $0x8,%esi
  1f691d:	mov    %r15,%rdi
  1f6920:	call   *0x1d2baa(%rip)        # 3c94d0 <_DYNAMIC+0x2a0>
  1f6926:	test   %rax,%rax
  1f6929:	je     1f7ad9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1529>
  1f692f:	xor    %edx,%edx
  1f6931:	mov    0x8(%rsp),%rcx
  1f6936:	jmp    1f697c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x3cc>
  1f6938:	xor    %edx,%edx
  1f693a:	mov    $0x8,%eax
  1f693f:	jmp    1f697a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x3ca>
  1f6941:	mov    $0x8,%eax
  1f6946:	mov    %rax,0x78(%rsp)
  1f694b:	mov    $0x1,%dl
  1f694d:	mov    $0x4,%eax
  1f6952:	mov    %rax,0x48(%rsp)
  1f6957:	mov    $0x1,%eax
  1f695c:	mov    %rax,0xa0(%rsp)
  1f6964:	movq   $0x0,0xa8(%rsp)
  1f6970:	mov    $0x8,%eax
  1f6975:	mov    %rax,0x60(%rsp)
  1f697a:	xor    %ecx,%ecx
  1f697c:	mov    %rcx,0x168(%rsp)
  1f6984:	mov    %rax,0x170(%rsp)
  1f698c:	movq   $0x0,0x178(%rsp)
  1f6998:	lea    0xb8(%rsp),%rax
  1f69a0:	mov    %rax,0x1e8(%rsp)
  1f69a8:	mov    0xb0(%rsp),%rax
  1f69b0:	mov    %rax,0x1f0(%rsp)
  1f69b8:	lea    0x1b0(%rsp),%rax
  1f69c0:	mov    %rax,0x1f8(%rsp)
  1f69c8:	lea    0x180(%rsp),%rax
  1f69d0:	mov    %rax,0x200(%rsp)
  1f69d8:	mov    %r14,0x208(%rsp)
  1f69e0:	mov    %rdx,0x20(%rsp)
  1f69e5:	test   %dl,%dl
  1f69e7:	mov    0x8(%rsp),%r14
  1f69ec:	jne    1f6a48 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x498>
  1f69ee:	xor    %r15d,%r15d
  1f69f1:	lea    0x260(%rsp),%r13
  1f69f9:	nopl   0x0(%rax)
  1f6a00:	sub    $0x8,%rsp
  1f6a04:	lea    0x1f0(%rsp),%rdi
  1f6a0c:	mov    %r15,%rsi
  1f6a0f:	mov    %r13,%rdx
  1f6a12:	mov    %rbp,%rcx
  1f6a15:	mov    %r14,%r8
  1f6a18:	mov    0x48(%rsp),%r9
  1f6a1d:	push   %r14
  1f6a1f:	push   0x58(%rsp)
  1f6a23:	push   %r14
  1f6a25:	push   0x48(%rsp)
  1f6a29:	push   %r14
  1f6a2b:	push   0x68(%rsp)
  1f6a2f:	push   %r14
  1f6a31:	push   0x70(%rsp)
  1f6a35:	push   %r14
  1f6a37:	call   1f4d10 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  1f6a3c:	add    $0x50,%rsp
  1f6a40:	inc    %r15
  1f6a43:	cmp    %r15,%r14
  1f6a46:	jne    1f6a00 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x450>
  1f6a48:	mov    0x3f8(%rsp),%rcx
  1f6a50:	movzbl 0x60(%rcx),%eax
  1f6a54:	test   %al,%al
  1f6a56:	je     1f6b92 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x5e2>
  1f6a5c:	xor    %r15d,%r15d
  1f6a5f:	movq   $0x0,0x70(%rsp)
  1f6a68:	mov    0x8(%rsp),%rbx
  1f6a6d:	mov    0x20(%rsp),%r14
  1f6a72:	mov    0x3f8(%rsp),%rax
  1f6a7a:	lock add %r15,0x50(%rax)
  1f6a7f:	mov    0x70(%rsp),%rcx
  1f6a84:	lock add %rcx,0x58(%rax)
  1f6a89:	mov    0x168(%rsp),%rsi
  1f6a91:	test   %rsi,%rsi
  1f6a94:	je     1f6aad <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4fd>
  1f6a96:	mov    0x170(%rsp),%rdi
  1f6a9e:	shl    $0x3,%rsi
  1f6aa2:	mov    $0x8,%edx
  1f6aa7:	call   *0x1d29f3(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f6aad:	test   %r14b,%r14b
  1f6ab0:	jne    1f6ac8 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x518>
  1f6ab2:	mov    $0x1,%edx
  1f6ab7:	mov    0xa0(%rsp),%rdi
  1f6abf:	mov    %rbx,%rsi
  1f6ac2:	call   *0x1d29d8(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f6ac8:	mov    0xa8(%rsp),%r14
  1f6ad0:	test   %r14,%r14
  1f6ad3:	je     1f6b03 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x553>
  1f6ad5:	shl    $0x3,%r14
  1f6ad9:	mov    0x1d29c0(%rip),%rbx        # 3c94a0 <_DYNAMIC+0x270>
  1f6ae0:	mov    $0x8,%edx
  1f6ae5:	mov    0x78(%rsp),%rdi
  1f6aea:	mov    %r14,%rsi
  1f6aed:	call   *%rbx
  1f6aef:	mov    $0x8,%edx
  1f6af4:	mov    0x60(%rsp),%rdi
  1f6af9:	mov    %r14,%rsi
  1f6afc:	call   *%rbx
  1f6afe:	mov    0x8(%rsp),%rbx
  1f6b03:	test   %rbx,%rbx
  1f6b06:	je     1f6b1f <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x56f>
  1f6b08:	shl    $0x2,%rbx
  1f6b0c:	mov    $0x4,%edx
  1f6b11:	mov    0x48(%rsp),%rdi
  1f6b16:	mov    %rbx,%rsi
  1f6b19:	call   *0x1d2981(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f6b1f:	mov    0x108(%rsp),%r14
  1f6b27:	test   %r14,%r14
  1f6b2a:	je     1f6b80 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x5d0>
  1f6b2c:	shl    $0x3,%r14
  1f6b30:	mov    0x1d2969(%rip),%rbx        # 3c94a0 <_DYNAMIC+0x270>
  1f6b37:	mov    $0x8,%edx
  1f6b3c:	mov    0x28(%rsp),%rdi
  1f6b41:	mov    %r14,%rsi
  1f6b44:	call   *%rbx
  1f6b46:	mov    $0x8,%edx
  1f6b4b:	mov    0x38(%rsp),%rdi
  1f6b50:	mov    %r14,%rsi
  1f6b53:	call   *%rbx
  1f6b55:	mov    $0x8,%edx
  1f6b5a:	mov    0x30(%rsp),%rdi
  1f6b5f:	mov    %r14,%rsi
  1f6b62:	call   *%rbx
  1f6b64:	mov    $0x8,%edx
  1f6b69:	mov    0x40(%rsp),%rdi
  1f6b6e:	mov    %r14,%rsi
  1f6b71:	call   *%rbx
  1f6b73:	mov    $0x8,%edx
  1f6b78:	mov    %rbp,%rdi
  1f6b7b:	mov    %r14,%rsi
  1f6b7e:	call   *%rbx
  1f6b80:	add    $0x3a8,%rsp
  1f6b87:	pop    %rbx
  1f6b88:	pop    %r12
  1f6b8a:	pop    %r13
  1f6b8c:	pop    %r14
  1f6b8e:	pop    %r15
  1f6b90:	pop    %rbp
  1f6b91:	ret
  1f6b92:	mov    0x3e0(%rsp),%rax
  1f6b9a:	mov    %rax,0xe0(%rsp)
  1f6ba2:	mov    0x3e8(%rsp),%rsi
  1f6baa:	mov    0x18(%rsi),%rdx
  1f6bae:	xor    %eax,%eax
  1f6bb0:	sub    0x20(%rsi),%eax
  1f6bb3:	and    $0x3f,%eax
  1f6bb6:	mov    %rax,0x1e0(%rsp)
  1f6bbe:	movzbl 0x28(%rsi),%eax
  1f6bc2:	mov    %al,0x17(%rsp)
  1f6bc6:	mov    %ebx,%eax
  1f6bc8:	neg    %eax
  1f6bca:	and    $0x3f,%eax
  1f6bcd:	mov    %rax,0x1d8(%rsp)
  1f6bd5:	lea    0x8(%rcx),%rax
  1f6bd9:	mov    %rax,0x1c0(%rsp)
  1f6be1:	lea    0x28(%rcx),%rax
  1f6be5:	mov    %rax,0x1c8(%rsp)
  1f6bed:	xor    %r15d,%r15d
  1f6bf0:	movq   $0x0,0x70(%rsp)
  1f6bf9:	mov    0x8(%rsp),%rbx
  1f6bfe:	mov    0x20(%rsp),%r14
  1f6c03:	mov    %rdx,0x110(%rsp)
  1f6c0b:	jmp    1f6c24 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x674>
  1f6c0d:	nopl   (%rax)
  1f6c10:	mov    0x110(%rsp),%rdx
  1f6c18:	movzbl 0x60(%rax),%eax
  1f6c1c:	test   %al,%al
  1f6c1e:	jne    1f6a72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4c2>
  1f6c24:	test   %rdx,%rdx
  1f6c27:	je     1f6c3e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x68e>
  1f6c29:	mov    0x3f8(%rsp),%rax
  1f6c31:	mov    0x50(%rax),%rax
  1f6c35:	cmp    %rdx,%rax
  1f6c38:	ja     1f6a72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4c2>
  1f6c3e:	test   %r14b,%r14b
  1f6c41:	je     1f7760 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11b0>
  1f6c47:	lea    0x210(%rsp),%rdi
  1f6c4f:	mov    0x60(%rsp),%rsi
  1f6c54:	mov    %rbx,%rdx
  1f6c57:	lea    0x168(%rsp),%rcx
  1f6c5f:	call   *0x1d2e6b(%rip)        # 3c9ad0 <_DYNAMIC+0x8a0>
  1f6c65:	test   %r14b,%r14b
  1f6c68:	jne    1f7815 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1265>
  1f6c6e:	xor    %r13d,%r13d
  1f6c71:	mov    %r15,0x88(%rsp)
  1f6c79:	mov    0x210(%rsp),%r14
  1f6c81:	mov    0x218(%rsp),%rax
  1f6c89:	mov    %rax,0x58(%rsp)
  1f6c8e:	mov    0x238(%rsp),%rax
  1f6c96:	mov    %rax,0x118(%rsp)
  1f6c9e:	mov    0x40(%rsp),%rax
  1f6ca3:	lea    (%rax,%r13,8),%rax
  1f6ca7:	mov    %rax,0x90(%rsp)
  1f6caf:	lea    0x0(%rbp,%r13,8),%rax
  1f6cb4:	mov    %rax,0x50(%rsp)
  1f6cb9:	mov    0x30(%rsp),%rax
  1f6cbe:	lea    (%rax,%r13,8),%rax
  1f6cc2:	mov    %rax,0x98(%rsp)
  1f6cca:	mov    0x38(%rsp),%rax
  1f6ccf:	lea    (%rax,%r13,8),%rax
  1f6cd3:	mov    %rax,0xf0(%rsp)
  1f6cdb:	jmp    1f6d4d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x79d>
  1f6cdd:	nopl   (%rax)
  1f6ce0:	sub    $0x8,%rsp
  1f6ce4:	lea    0x1f0(%rsp),%rdi
  1f6cec:	mov    %r13,%rsi
  1f6cef:	lea    0x268(%rsp),%rdx
  1f6cf7:	mov    %rbp,%rcx
  1f6cfa:	mov    %rbx,%r8
  1f6cfd:	mov    0x48(%rsp),%r9
  1f6d02:	push   %rbx
  1f6d03:	push   0x58(%rsp)
  1f6d07:	push   %rbx
  1f6d08:	push   0x48(%rsp)
  1f6d0c:	push   %rbx
  1f6d0d:	push   0x68(%rsp)
  1f6d11:	push   %rbx
  1f6d12:	push   0x70(%rsp)
  1f6d16:	push   %rbx
  1f6d17:	call   1f4d10 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  1f6d1c:	add    $0x50,%rsp
  1f6d20:	inc    %r13
  1f6d23:	addq   $0x8,0x90(%rsp)
  1f6d2c:	addq   $0x8,0x50(%rsp)
  1f6d32:	addq   $0x8,0x98(%rsp)
  1f6d3b:	addq   $0x8,0xf0(%rsp)
  1f6d44:	cmp    %rbx,%r13
  1f6d47:	jae    1f7810 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1260>
  1f6d4d:	cmp    %rbx,%r13
  1f6d50:	jae    1f7969 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13b9>
  1f6d56:	mov    0xa0(%rsp),%rax
  1f6d5e:	cmpb   $0x0,(%rax,%r13,1)
  1f6d63:	jne    1f6ce0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x730>
  1f6d69:	mov    0x78(%rsp),%rax
  1f6d6e:	mov    (%rax,%r13,8),%r11
  1f6d72:	cmp    0x3e0(%rsp),%r11
  1f6d7a:	jae    1f7990 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13e0>
  1f6d80:	shl    $0x5,%r11
  1f6d84:	mov    0x40(%rsp),%rax
  1f6d89:	mov    (%rax,%r13,8),%r9
  1f6d8d:	mov    0x8(%r12,%r11,1),%rax
  1f6d92:	mov    $0x0,%r10d
  1f6d98:	mov    %r14,%rcx
  1f6d9b:	sub    %r9,%rax
  1f6d9e:	jb     1f6da5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x7f5>
  1f6da0:	mov    $0x0,%ecx
  1f6da5:	add    %rcx,%rax
  1f6da8:	mov    0x60(%rsp),%rcx
  1f6dad:	mulq   (%rcx,%r13,8)
  1f6db1:	mov    %rax,%rcx
  1f6db4:	mov    %rdx,%rdi
  1f6db7:	mov    0x58(%rsp),%rbx
  1f6dbc:	mov    %rbx,%rax
  1f6dbf:	imul   %rcx,%rax
  1f6dc3:	mul    %r14
  1f6dc6:	mov    %rdx,%rsi
  1f6dc9:	add    %rcx,%rax
  1f6dcc:	adc    %rdi,%rsi
  1f6dcf:	mov    $0x0,%eax
  1f6dd4:	cmp    %rsi,%r14
  1f6dd7:	ja     1f6ddc <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x82c>
  1f6dd9:	mov    %r14,%rax
  1f6ddc:	sub    %rax,%rsi
  1f6ddf:	mov    %rsi,%rax
  1f6de2:	mul    %rsi
  1f6de5:	mov    %rax,%rcx
  1f6de8:	mov    %rdx,%rdi
  1f6deb:	mov    %rbx,%rax
  1f6dee:	imul   %rcx,%rax
  1f6df2:	mul    %r14
  1f6df5:	add    %rcx,%rax
  1f6df8:	adc    %rdi,%rdx
  1f6dfb:	mov    $0x0,%eax
  1f6e00:	cmp    %rdx,%r14
  1f6e03:	ja     1f6e08 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x858>
  1f6e05:	mov    %r14,%rax
  1f6e08:	sub    %rax,%rdx
  1f6e0b:	mov    0x0(%rbp,%r13,8),%rdi
  1f6e10:	mov    %r14,%rax
  1f6e13:	sub    %rdi,%rdx
  1f6e16:	jb     1f6e1d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x86d>
  1f6e18:	mov    $0x0,%eax
  1f6e1d:	add    %rdx,%rax
  1f6e20:	mov    %r14,%rcx
  1f6e23:	sub    (%r12,%r11,1),%rax
  1f6e27:	jb     1f6e2e <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x87e>
  1f6e29:	mov    $0x0,%ecx
  1f6e2e:	lea    (%rcx,%rax,1),%r8
  1f6e32:	mov    %rdi,%rax
  1f6e35:	mov    %r14,%rcx
  1f6e38:	sub    %r8,%rax
  1f6e3b:	jb     1f6e42 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x892>
  1f6e3d:	mov    $0x0,%ecx
  1f6e42:	add    %rcx,%rax
  1f6e45:	mul    %rsi
  1f6e48:	mov    %rax,%rcx
  1f6e4b:	mov    %rdx,%rsi
  1f6e4e:	mov    %rbx,%rax
  1f6e51:	imul   %rcx,%rax
  1f6e55:	mul    %r14
  1f6e58:	add    %rcx,%rax
  1f6e5b:	adc    %rsi,%rdx
  1f6e5e:	mov    $0x0,%eax
  1f6e63:	cmp    %rdx,%r14
  1f6e66:	ja     1f6e6b <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x8bb>
  1f6e68:	mov    %r14,%rax
  1f6e6b:	sub    %rax,%rdx
  1f6e6e:	mov    %r14,%rcx
  1f6e71:	sub    %r9,%rdx
  1f6e74:	jb     1f6e7b <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x8cb>
  1f6e76:	mov    $0x0,%ecx
  1f6e7b:	mov    0x30(%rsp),%rax
  1f6e80:	mov    (%rax,%r13,8),%rax
  1f6e84:	mov    0x10(%r12,%r11,1),%rbx
  1f6e89:	add    %rax,%rbx
  1f6e8c:	mov    0xb8(%rsp),%r15
  1f6e94:	mov    $0x0,%ebp
  1f6e99:	cmp    %r15,%rbx
  1f6e9c:	jb     1f6ea1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x8f1>
  1f6e9e:	mov    %r15,%rbp
  1f6ea1:	mov    0x38(%rsp),%rsi
  1f6ea6:	mov    (%rsi,%r13,8),%rsi
  1f6eaa:	mov    0x18(%r12,%r11,1),%r12
  1f6eaf:	add    %rsi,%r12
  1f6eb2:	cmp    %r15,%r12
  1f6eb5:	jb     1f6eba <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x90a>
  1f6eb7:	mov    %r15,%r10
  1f6eba:	add    %rdx,%rcx
  1f6ebd:	sub    %rbp,%rbx
  1f6ec0:	sub    %r10,%r12
  1f6ec3:	cmpb   $0x0,0x17(%rsp)
  1f6ec8:	je     1f6fac <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x9fc>
  1f6ece:	cmp    0x118(%rsp),%rcx
  1f6ed6:	mov    0x18(%rsp),%rbp
  1f6edb:	jbe    1f6f00 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x950>
  1f6edd:	mov    %r14,%rdx
  1f6ee0:	sub    %rcx,%rdx
  1f6ee3:	mov    %r15,%rcx
  1f6ee6:	sub    %rbx,%rcx
  1f6ee9:	test   %rbx,%rbx
  1f6eec:	cmovne %rcx,%rbx
  1f6ef0:	mov    %r15,%rcx
  1f6ef3:	sub    %r12,%rcx
  1f6ef6:	test   %r12,%r12
  1f6ef9:	cmovne %rcx,%r12
  1f6efd:	mov    %rdx,%rcx
  1f6f00:	cmp    0x8(%rsp),%r13
  1f6f05:	jae    1f797a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13ca>
  1f6f0b:	mov    0x28(%rsp),%rdx
  1f6f10:	cmp    (%rdx,%r13,8),%r8
  1f6f14:	jne    1f6fd0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xa20>
  1f6f1a:	cmp    %r8,%rdi
  1f6f1d:	jae    1f6f2b <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x97b>
  1f6f1f:	mov    %rdi,%r8
  1f6f22:	mov    %r9,%rcx
  1f6f25:	mov    %rax,%rbx
  1f6f28:	mov    %rsi,%r12
  1f6f2b:	lea    0x120(%rsp),%rdi
  1f6f33:	mov    0xb0(%rsp),%rsi
  1f6f3b:	mov    %r8,%rdx
  1f6f3e:	call   *0x1d2b9c(%rip)        # 3c9ae0 <_DYNAMIC+0x8b0>
  1f6f44:	cmpb   $0x0,0x120(%rsp)
  1f6f4c:	jne    1f7010 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xa60>
  1f6f52:	sub    $0x8,%rsp
  1f6f56:	lea    0x1f0(%rsp),%rdi
  1f6f5e:	mov    %r13,%rsi
  1f6f61:	lea    0x268(%rsp),%rdx
  1f6f69:	mov    %rbp,%rcx
  1f6f6c:	mov    0x10(%rsp),%rbx
  1f6f71:	mov    %rbx,%r8
  1f6f74:	mov    0x48(%rsp),%r9
  1f6f79:	push   %rbx
  1f6f7a:	push   0x58(%rsp)
  1f6f7e:	push   %rbx
  1f6f7f:	push   0x48(%rsp)
  1f6f83:	push   %rbx
  1f6f84:	push   0x68(%rsp)
  1f6f88:	push   %rbx
  1f6f89:	push   0x70(%rsp)
  1f6f8d:	push   %rbx
  1f6f8e:	call   1f4d10 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  1f6f93:	add    $0x50,%rsp
  1f6f97:	mov    0x80(%rsp),%r12
  1f6f9f:	mov    0x88(%rsp),%r15
  1f6fa7:	jmp    1f6d20 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x770>
  1f6fac:	mov    0x18(%rsp),%rbp
  1f6fb1:	cmp    0x8(%rsp),%r13
  1f6fb6:	jb     1f6f0b <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x95b>
  1f6fbc:	jmp    1f797a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13ca>
  1f6fc1:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  1f6fd0:	mov    0x28(%rsp),%rax
  1f6fd5:	mov    %rdi,(%rax,%r13,8)
  1f6fd9:	mov    %r8,0x0(%rbp,%r13,8)
  1f6fde:	mov    0x40(%rsp),%rax
  1f6fe3:	mov    %rcx,(%rax,%r13,8)
  1f6fe7:	mov    0x30(%rsp),%rax
  1f6fec:	mov    %rbx,(%rax,%r13,8)
  1f6ff0:	mov    0x38(%rsp),%rax
  1f6ff5:	mov    %r12,(%rax,%r13,8)
  1f6ff9:	mov    0x8(%rsp),%rbx
  1f6ffe:	mov    0x20(%rsp),%r14
  1f7003:	jmp    1f70ce <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb1e>
  1f7008:	nopl   0x0(%rax,%rax,1)
  1f7010:	mov    0x128(%rsp),%rdx
  1f7018:	mov    0x130(%rsp),%rsi
  1f7020:	add    %rbx,%rbx
  1f7023:	cmp    %r15,%rbx
  1f7026:	mov    %r15,%rax
  1f7029:	mov    $0x0,%edi
  1f702e:	cmovb  %rdi,%rax
  1f7032:	mov    %rbx,%rcx
  1f7035:	sub    %rax,%rcx
  1f7038:	add    %r12,%r12
  1f703b:	cmp    %r15,%r12
  1f703e:	cmovb  %rdi,%r15
  1f7042:	sub    %r15,%r12
  1f7045:	mov    %rdx,0x0(%rbp,%r13,8)
  1f704a:	cmp    0x118(%rsp),%rsi
  1f7052:	seta   %dl
  1f7055:	sub    %rsi,%r14
  1f7058:	test   %dl,0x17(%rsp)
  1f705c:	cmove  %rsi,%r14
  1f7060:	mov    0x40(%rsp),%rdx
  1f7065:	mov    %r14,(%rdx,%r13,8)
  1f7069:	je     1f70a2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xaf2>
  1f706b:	mov    0xb8(%rsp),%rdx
  1f7073:	sub    %rcx,%rdx
  1f7076:	cmp    %rax,%rbx
  1f7079:	cmove  %rdi,%rdx
  1f707d:	mov    0x30(%rsp),%rax
  1f7082:	mov    %rdx,(%rax,%r13,8)
  1f7086:	mov    $0x0,%eax
  1f708b:	test   %r12,%r12
  1f708e:	mov    0x20(%rsp),%r14
  1f7093:	je     1f70b3 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb03>
  1f7095:	mov    0xb8(%rsp),%rax
  1f709d:	sub    %r12,%rax
  1f70a0:	jmp    1f70b3 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xb03>
  1f70a2:	mov    0x30(%rsp),%rax
  1f70a7:	mov    %rcx,(%rax,%r13,8)
  1f70ab:	mov    %r12,%rax
  1f70ae:	mov    0x20(%rsp),%r14
  1f70b3:	mov    0x38(%rsp),%rcx
  1f70b8:	mov    %rax,(%rcx,%r13,8)
  1f70bc:	mov    0x28(%rsp),%rax
  1f70c1:	movq   $0xffffffffffffffff,(%rax,%r13,8)
  1f70c9:	mov    0x8(%rsp),%rbx
  1f70ce:	mov    0x88(%rsp),%r15
  1f70d6:	mov    0x48(%rsp),%rcx
  1f70db:	mov    (%rcx,%r13,4),%eax
  1f70df:	inc    %eax
  1f70e1:	mov    %eax,(%rcx,%r13,4)
  1f70e5:	cmpl   $0x0,0x6c(%rsp)
  1f70ea:	mov    0x80(%rsp),%r12
  1f70f2:	je     1f7190 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xbe0>
  1f70f8:	mov    0x50(%rsp),%rcx
  1f70fd:	mov    (%rcx),%rdx
  1f7100:	rol    $0x17,%rdx
  1f7104:	movabs $0xd6e8feb86659fd93,%rcx
  1f710e:	imul   %rcx,%rdx
  1f7112:	mov    0x1d8(%rsp),%rcx
  1f711a:	shr    %cl,%rdx
  1f711d:	test   %rdx,%rdx
  1f7120:	je     1f7190 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xbe0>
  1f7122:	cmp    0x1d0(%rsp),%eax
  1f7129:	jbe    1f7170 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xbc0>
  1f712b:	sub    $0x8,%rsp
  1f712f:	lea    0x1f0(%rsp),%rdi
  1f7137:	mov    %r13,%rsi
  1f713a:	lea    0x268(%rsp),%rdx
  1f7142:	mov    %rbp,%rcx
  1f7145:	mov    %rbx,%r8
  1f7148:	mov    0x48(%rsp),%r9
  1f714d:	push   %rbx
  1f714e:	push   0x58(%rsp)
  1f7152:	push   %rbx
  1f7153:	push   0x48(%rsp)
  1f7157:	push   %rbx
  1f7158:	push   0x68(%rsp)
  1f715c:	push   %rbx
  1f715d:	push   0x70(%rsp)
  1f7161:	push   %rbx
  1f7162:	call   1f4d10 <crypto_lib::cryptanalysis::prime_fast::rho_worker::{closure#0}>
  1f7167:	add    $0x50,%rsp
  1f716b:	nopl   0x0(%rax,%rax,1)
  1f7170:	inc    %r13
  1f7173:	cmp    %rbx,%r13
  1f7176:	jb     1f6c79 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x6c9>
  1f717c:	jmp    1f7815 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1265>
  1f7181:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  1f7190:	xor    %eax,%eax
  1f7192:	mov    0x3f8(%rsp),%rcx
  1f719a:	mov    $0x1,%edx
  1f719f:	lock cmpxchg %edx,(%rcx)
  1f71a3:	jne    1f76f9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1149>
  1f71a9:	mov    0x1d2b20(%rip),%rax        # 3c9cd0 <_DYNAMIC+0xaa0>
  1f71b0:	mov    (%rax),%rax
  1f71b3:	movabs $0x7fffffffffffffff,%rcx
  1f71bd:	test   %rcx,%rax
  1f71c0:	jne    1f770c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x115c>
  1f71c6:	movl   $0x0,0x58(%rsp)
  1f71ce:	mov    0x3f8(%rsp),%rcx
  1f71d6:	movzbl 0x4(%rcx),%eax
  1f71da:	test   %al,%al
  1f71dc:	jne    1f7920 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1370>
  1f71e2:	mov    0x70(%rsp),%r15
  1f71e7:	inc    %r15
  1f71ea:	cmpq   $0x0,0x20(%rcx)
  1f71ef:	lea    0xc0(%rsp),%r12
  1f71f7:	mov    %r15,0x70(%rsp)
  1f71fc:	je     1f72ad <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xcfd>
  1f7202:	mov    %rcx,%rbx
  1f7205:	mov    0x1c8(%rsp),%rdi
  1f720d:	mov    0x50(%rsp),%rsi
  1f7212:	call   2b50a0 <<std::hash::random::RandomState as core::hash::BuildHasher>::hash_one::<&crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolMono>>
  1f7217:	mov    %rax,%rsi
  1f721a:	shr    $0x39,%rsi
  1f721e:	mov    0x8(%rbx),%rcx
  1f7222:	mov    0x10(%rbx),%rdx
  1f7226:	movd   %esi,%xmm0
  1f722a:	punpcklbw %xmm0,%xmm0
  1f722e:	pshuflw $0x0,%xmm0,%xmm0
  1f7233:	pshufd $0x44,%xmm0,%xmm0
  1f7238:	xor    %esi,%esi
  1f723a:	and    %rdx,%rax
  1f723d:	movdqu (%rcx,%rax,1),%xmm1
  1f7242:	movdqa %xmm1,%xmm2
  1f7246:	pcmpeqb %xmm0,%xmm2
  1f724a:	pmovmskb %xmm2,%edi
  1f724e:	test   %edi,%edi
  1f7250:	je     1f7290 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xce0>
  1f7252:	mov    0x50(%rsp),%r8
  1f7257:	mov    (%r8),%r8
  1f725a:	tzcnt  %edi,%r10d
  1f725f:	add    %rax,%r10
  1f7262:	and    %rdx,%r10
  1f7265:	shl    $0x5,%r10
  1f7269:	mov    %rcx,%r9
  1f726c:	sub    %r10,%r9
  1f726f:	cmp    -0x20(%r9),%r8
  1f7273:	je     1f7310 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xd60>
  1f7279:	lea    -0x1(%rdi),%r9d
  1f727d:	and    %di,%r9w
  1f7281:	mov    %r9d,%edi
  1f7284:	jne    1f725a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xcaa>
  1f7286:	cs nopw 0x0(%rax,%rax,1)
  1f7290:	pcmpeqd %xmm2,%xmm2
  1f7294:	pcmpeqb %xmm2,%xmm1
  1f7298:	pmovmskb %xmm1,%edi
  1f729c:	test   %edi,%edi
  1f729e:	jne    1f72ad <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xcfd>
  1f72a0:	add    %rsi,%rax
  1f72a3:	add    $0x10,%rax
  1f72a7:	add    $0x10,%rsi
  1f72ab:	jmp    1f723a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xc8a>
  1f72ad:	mov    0x50(%rsp),%rax
  1f72b2:	mov    (%rax),%rdx
  1f72b5:	mov    0x90(%rsp),%rax
  1f72bd:	mov    (%rax),%rax
  1f72c0:	mov    0x98(%rsp),%rcx
  1f72c8:	mov    (%rcx),%rcx
  1f72cb:	mov    0xf0(%rsp),%rsi
  1f72d3:	mov    (%rsi),%rsi
  1f72d6:	mov    %rax,0xc0(%rsp)
  1f72de:	mov    %rcx,0xc8(%rsp)
  1f72e6:	mov    %rsi,0xd0(%rsp)
  1f72ee:	lea    0x120(%rsp),%rdi
  1f72f6:	mov    0x1c0(%rsp),%rsi
  1f72fe:	mov    %r12,%rcx
  1f7301:	call   2d3780 <<hashbrown::map::HashMap<u64, (u64, u64, u64), std::hash::random::RandomState>>::insert>
  1f7306:	mov    0x20(%rsp),%r14
  1f730b:	jmp    1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f7310:	mov    -0x18(%r9),%r12
  1f7314:	mov    -0x10(%r9),%rbx
  1f7318:	mov    -0x8(%r9),%rcx
  1f731c:	mov    0x98(%rsp),%rax
  1f7324:	mov    (%rax),%rbp
  1f7327:	mov    0xf0(%rsp),%rax
  1f732f:	mov    (%rax),%r15
  1f7332:	cmp    %rbx,%rbp
  1f7335:	jne    1f7340 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xd90>
  1f7337:	cmp    %r15,%rcx
  1f733a:	je     1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f7340:	mov    %rcx,0x50(%rsp)
  1f7345:	mov    0x90(%rsp),%rax
  1f734d:	mov    (%rax),%r14
  1f7350:	movdqa 0x1b0(%rsp),%xmm0
  1f7359:	movdqa %xmm0,0xf0(%rsp)
  1f7362:	mov    0x180(%rsp),%rax
  1f736a:	mov    %rax,0x98(%rsp)
  1f7372:	mov    0x188(%rsp),%rax
  1f737a:	mov    %rax,0x90(%rsp)
  1f7382:	lea    0x120(%rsp),%rdi
  1f738a:	mov    0xe8(%rsp),%rsi
  1f7392:	call   *0x1d2728(%rip)        # 3c9ac0 <_DYNAMIC+0x890>
  1f7398:	cmp    %r12,%r14
  1f739b:	jne    1f73e8 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xe38>
  1f739d:	sub    %rbx,%rbp
  1f73a0:	mov    0x120(%rsp),%rax
  1f73a8:	mov    $0x0,%ecx
  1f73ad:	cmovb  %rax,%rcx
  1f73b1:	add    %rbp,%rcx
  1f73b4:	mov    0x50(%rsp),%rsi
  1f73b9:	sub    %r15,%rsi
  1f73bc:	mov    $0x0,%edx
  1f73c1:	cmovae %rdx,%rax
  1f73c5:	lea    (%rsi,%rax,1),%r15
  1f73c9:	mov    0xe8(%rsp),%r11
  1f73d1:	mov    0x20(%rsp),%r14
  1f73d6:	test   %r15,%r15
  1f73d9:	lea    0xc0(%rsp),%rbp
  1f73e1:	jne    1f7438 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xe88>
  1f73e3:	jmp    1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f73e8:	add    %rbx,%rbp
  1f73eb:	mov    0xe8(%rsp),%r11
  1f73f3:	mov    %r11,%rax
  1f73f6:	sub    %rbp,%rax
  1f73f9:	mov    %r11,%rcx
  1f73fc:	mov    $0x0,%edx
  1f7401:	cmova  %rdx,%rcx
  1f7405:	cmp    %rcx,%rbp
  1f7408:	lea    (%rax,%rcx,1),%rcx
  1f740c:	cmove  %rdx,%rcx
  1f7410:	add    0x50(%rsp),%r15
  1f7415:	cmp    %r11,%r15
  1f7418:	mov    %r11,%rax
  1f741b:	cmovb  %rdx,%rax
  1f741f:	sub    %rax,%r15
  1f7422:	mov    0x20(%rsp),%r14
  1f7427:	test   %r15,%r15
  1f742a:	lea    0xc0(%rsp),%rbp
  1f7432:	je     1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f7438:	cmp    $0x1,%r15
  1f743c:	sete   %dil
  1f7440:	mov    $0x1,%edx
  1f7445:	cmp    $0x1,%r11
  1f7449:	je     1f753d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf8d>
  1f744f:	mov    $0x0,%eax
  1f7454:	cmp    $0x1,%r15
  1f7458:	mov    0xb0(%rsp),%rbx
  1f7460:	je     1f7547 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf97>
  1f7466:	mov    %r11,%rsi
  1f7469:	xor    %eax,%eax
  1f746b:	nopl   0x0(%rax,%rax,1)
  1f7470:	test   $0x1,%r15b
  1f7474:	jne    1f74af <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xeff>
  1f7476:	cs nopw 0x0(%rax,%rax,1)
  1f7480:	mov    %r15,%rdi
  1f7483:	shr    $1,%rdi
  1f7486:	mov    %edx,%r8d
  1f7489:	mov    %rdx,%r9
  1f748c:	shr    $1,%r9
  1f748f:	mov    %rdx,%r10
  1f7492:	and    %r11,%r10
  1f7495:	xor    %r11,%rdx
  1f7498:	shr    $1,%rdx
  1f749b:	add    %r10,%rdx
  1f749e:	test   $0x1,%r8b
  1f74a2:	cmove  %r9,%rdx
  1f74a6:	test   $0x2,%r15b
  1f74aa:	mov    %rdi,%r15
  1f74ad:	je     1f7480 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xed0>
  1f74af:	test   $0x1,%sil
  1f74b3:	jne    1f74ef <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf3f>
  1f74b5:	data16 cs nopw 0x0(%rax,%rax,1)
  1f74c0:	mov    %rsi,%rdi
  1f74c3:	shr    $1,%rdi
  1f74c6:	mov    %eax,%r8d
  1f74c9:	mov    %rax,%r9
  1f74cc:	shr    $1,%r9
  1f74cf:	mov    %rax,%r10
  1f74d2:	and    %r11,%r10
  1f74d5:	xor    %r11,%rax
  1f74d8:	shr    $1,%rax
  1f74db:	add    %r10,%rax
  1f74de:	test   $0x1,%r8b
  1f74e2:	cmove  %r9,%rax
  1f74e6:	test   $0x2,%sil
  1f74ea:	mov    %rdi,%rsi
  1f74ed:	je     1f74c0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf10>
  1f74ef:	mov    %r15,%rdi
  1f74f2:	sub    %rsi,%rdi
  1f74f5:	jae    1f7510 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf60>
  1f74f7:	sub    %r15,%rsi
  1f74fa:	sub    %rdx,%rax
  1f74fd:	mov    $0x0,%edi
  1f7502:	cmovb  %r11,%rdi
  1f7506:	add    %rdi,%rax
  1f7509:	jmp    1f7523 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf73>
  1f750b:	nopl   0x0(%rax,%rax,1)
  1f7510:	sub    %rax,%rdx
  1f7513:	mov    $0x0,%r8d
  1f7519:	cmovb  %r11,%r8
  1f751d:	add    %r8,%rdx
  1f7520:	mov    %rdi,%r15
  1f7523:	cmp    $0x1,%r15
  1f7527:	sete   %dil
  1f752b:	cmp    $0x1,%rsi
  1f752f:	je     1f7547 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf97>
  1f7531:	cmp    $0x1,%r15
  1f7535:	jne    1f7470 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xec0>
  1f753b:	jmp    1f7547 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xf97>
  1f753d:	xor    %eax,%eax
  1f753f:	mov    0xb0(%rsp),%rbx
  1f7547:	test   %r11,%r11
  1f754a:	je     1f7ae9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1539>
  1f7550:	test   %dil,%dil
  1f7553:	cmovne %rdx,%rax
  1f7557:	mul    %rcx
  1f755a:	mov    %rax,%rdi
  1f755d:	mov    %rdx,%rsi
  1f7560:	mov    %r11,%rdx
  1f7563:	xor    %ecx,%ecx
  1f7565:	call   *0x1d264d(%rip)        # 3c9bb8 <_DYNAMIC+0x988>
  1f756b:	mov    %rax,%r15
  1f756e:	movdqa 0xf0(%rsp),%xmm0
  1f7577:	movdqu %xmm0,0x158(%rsp)
  1f7580:	movq   $0x1,0x150(%rsp)
  1f758c:	movq   $0x0,0x190(%rsp)
  1f7598:	or     %rax,%rdx
  1f759b:	mov    %rax,%r12
  1f759e:	jne    1f7613 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1063>
  1f75a0:	jmp    1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f75a5:	data16 cs nopw 0x0(%rax,%rax,1)
  1f75b0:	mov    0x160(%rsp),%rax
  1f75b8:	mov    %rax,0xd0(%rsp)
  1f75c0:	movaps 0x150(%rsp),%xmm0
  1f75c8:	movaps %xmm0,0xc0(%rsp)
  1f75d0:	lea    0x248(%rsp),%rdi
  1f75d8:	mov    %rbx,%rsi
  1f75db:	mov    %rbp,%rdx
  1f75de:	lea    0x150(%rsp),%rcx
  1f75e6:	call   *0x1d2504(%rip)        # 3c9af0 <_DYNAMIC+0x8c0>
  1f75ec:	movdqu 0x248(%rsp),%xmm0
  1f75f5:	movdqa %xmm0,0x150(%rsp)
  1f75fe:	mov    0x258(%rsp),%rax
  1f7606:	mov    %rax,0x160(%rsp)
  1f760e:	shr    $1,%r12
  1f7611:	je     1f765c <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10ac>
  1f7613:	test   $0x1,%r12b
  1f7617:	je     1f75b0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1000>
  1f7619:	mov    %rbp,%rdi
  1f761c:	mov    %rbx,%rsi
  1f761f:	lea    0x190(%rsp),%rdx
  1f7627:	lea    0x150(%rsp),%rcx
  1f762f:	call   *0x1d24bb(%rip)        # 3c9af0 <_DYNAMIC+0x8c0>
  1f7635:	mov    0xd0(%rsp),%rax
  1f763d:	mov    %rax,0x1a0(%rsp)
  1f7645:	movdqu 0xc0(%rsp),%xmm0
  1f764e:	movdqa %xmm0,0x190(%rsp)
  1f7657:	jmp    1f75b0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1000>
  1f765c:	testb  $0x1,0x190(%rsp)
  1f7664:	je     1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f7666:	mov    0x98(%rsp),%rax
  1f766e:	cmp    %rax,0x198(%rsp)
  1f7676:	jne    1f7690 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x10e0>
  1f7678:	mov    0x90(%rsp),%rax
  1f7680:	cmp    %rax,0x1a0(%rsp)
  1f7688:	je     1f784d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x129d>
  1f768e:	xchg   %ax,%ax
  1f7690:	cmpb   $0x0,0x58(%rsp)
  1f7695:	mov    0x18(%rsp),%rbp
  1f769a:	mov    0x8(%rsp),%rbx
  1f769f:	mov    0x80(%rsp),%r12
  1f76a7:	mov    0x88(%rsp),%r15
  1f76af:	jne    1f76ca <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x111a>
  1f76b1:	mov    0x1d2618(%rip),%rax        # 3c9cd0 <_DYNAMIC+0xaa0>
  1f76b8:	mov    (%rax),%rax
  1f76bb:	movabs $0x7fffffffffffffff,%rcx
  1f76c5:	test   %rcx,%rax
  1f76c8:	jne    1f7741 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1191>
  1f76ca:	xor    %eax,%eax
  1f76cc:	mov    0x3f8(%rsp),%rcx
  1f76d4:	xchg   %eax,(%rcx)
  1f76d6:	cmp    $0x2,%eax
  1f76d9:	je     1f7731 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1181>
  1f76db:	mov    0x48(%rsp),%rax
  1f76e0:	movl   $0x0,(%rax,%r13,4)
  1f76e8:	inc    %r13
  1f76eb:	cmp    %rbx,%r13
  1f76ee:	jb     1f6c79 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x6c9>
  1f76f4:	jmp    1f7815 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1265>
  1f76f9:	mov    0x3f8(%rsp),%rdi
  1f7701:	call   *0x1d2661(%rip)        # 3c9d68 <_DYNAMIC+0xb38>
  1f7707:	jmp    1f71a9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xbf9>
  1f770c:	call   *0x1d25ce(%rip)        # 3c9ce0 <_DYNAMIC+0xab0>
  1f7712:	xor    $0x1,%al
  1f7714:	mov    %eax,0x58(%rsp)
  1f7718:	mov    0x3f8(%rsp),%rcx
  1f7720:	movzbl 0x4(%rcx),%eax
  1f7724:	test   %al,%al
  1f7726:	je     1f71e2 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0xc32>
  1f772c:	jmp    1f7920 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1370>
  1f7731:	mov    0x3f8(%rsp),%rdi
  1f7739:	call   *0x1d2599(%rip)        # 3c9cd8 <_DYNAMIC+0xaa8>
  1f773f:	jmp    1f76db <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x112b>
  1f7741:	call   *0x1d2599(%rip)        # 3c9ce0 <_DYNAMIC+0xab0>
  1f7747:	test   %al,%al
  1f7749:	jne    1f76ca <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x111a>
  1f774f:	mov    0x3f8(%rsp),%rax
  1f7757:	movb   $0x1,0x4(%rax)
  1f775b:	jmp    1f76ca <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x111a>
  1f7760:	mov    0x210(%rsp),%rax
  1f7768:	mov    0x220(%rsp),%rdx
  1f7770:	xor    %esi,%esi
  1f7772:	jmp    1f77a1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11f1>
  1f7774:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  1f7780:	mov    0xa0(%rsp),%rdi
  1f7788:	sete   (%rdi,%rsi,1)
  1f778c:	mov    0x60(%rsp),%rdi
  1f7791:	mov    %r9,(%rdi,%rsi,8)
  1f7795:	mov    %rcx,%rsi
  1f7798:	cmp    %rcx,%rbx
  1f779b:	je     1f6c47 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x697>
  1f77a1:	mov    0x0(%rbp,%rsi,8),%r13
  1f77a6:	movabs $0x9e3779b97f4a7c15,%rcx
  1f77b0:	imul   %rcx,%r13
  1f77b4:	mov    0x1e0(%rsp),%rcx
  1f77bc:	shr    %cl,%r13
  1f77bf:	mov    0x78(%rsp),%rcx
  1f77c4:	mov    %r13,(%rcx,%rsi,8)
  1f77c8:	cmp    0x3e0(%rsp),%r13
  1f77d0:	jae    1f7960 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13b0>
  1f77d6:	shl    $0x5,%r13
  1f77da:	mov    (%r12,%r13,1),%r9
  1f77de:	mov    0x0(%rbp,%rsi,8),%rdi
  1f77e3:	mov    %rax,%r8
  1f77e6:	cmp    %rdi,%r9
  1f77e9:	jb     1f77ee <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x123e>
  1f77eb:	xor    %r8d,%r8d
  1f77ee:	lea    0x1(%rsi),%rcx
  1f77f2:	add    %r9,%r8
  1f77f5:	mov    %rdx,%r9
  1f77f8:	sub    %rdi,%r8
  1f77fb:	je     1f7780 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11d0>
  1f77fd:	mov    %r8,%r9
  1f7800:	jmp    1f7780 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x11d0>
  1f7805:	data16 cs nopw 0x0(%rax,%rax,1)
  1f7810:	mov    0x20(%rsp),%r14
  1f7815:	add    %rbx,%r15
  1f7818:	mov    0x3f8(%rsp),%rax
  1f7820:	cmp    $0xffff,%r15
  1f7827:	jbe    1f6c10 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x660>
  1f782d:	lock add %r15,0x50(%rax)
  1f7832:	mov    0x70(%rsp),%rcx
  1f7837:	lock add %rcx,0x58(%rax)
  1f783c:	movq   $0x0,0x70(%rsp)
  1f7845:	xor    %r15d,%r15d
  1f7848:	jmp    1f6c10 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x660>
  1f784d:	mov    0x3f8(%rsp),%rax
  1f7855:	lea    0x38(%rax),%rsi
  1f7859:	lea    0x120(%rsp),%rdi
  1f7861:	call   2c06f0 <<std::sync::poison::mutex::Mutex<core::option::Option<alloc::vec::Vec<num_bigint::biguint::BigUint>>>>::lock>
  1f7866:	movabs $0x7fffffffffffffff,%rbx
  1f7870:	cmpl   $0x1,0x120(%rsp)
  1f7878:	je     1f79b5 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1405>
  1f787e:	mov    0x128(%rsp),%r12
  1f7886:	cmpb   $0x0,0x130(%rsp)
  1f788e:	movq   $0x1,0x8(%r12)
  1f7897:	mov    %r15,0x10(%r12)
  1f789c:	jne    1f78b1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1301>
  1f789e:	mov    0x1d242b(%rip),%rax        # 3c9cd0 <_DYNAMIC+0xaa0>
  1f78a5:	mov    (%rax),%rax
  1f78a8:	test   %rbx,%rax
  1f78ab:	jne    1f7a37 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1487>
  1f78b1:	xor    %eax,%eax
  1f78b3:	xchg   %eax,(%r12)
  1f78b7:	cmp    $0x2,%eax
  1f78ba:	je     1f7a01 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1451>
  1f78c0:	mov    0x3f8(%rsp),%rax
  1f78c8:	movb   $0x1,0x60(%rax)
  1f78cc:	cmpb   $0x0,0x58(%rsp)
  1f78d1:	jne    1f78e6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1336>
  1f78d3:	mov    0x1d23f6(%rip),%rax        # 3c9cd0 <_DYNAMIC+0xaa0>
  1f78da:	mov    (%rax),%rax
  1f78dd:	test   %rbx,%rax
  1f78e0:	jne    1f7a50 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x14a0>
  1f78e6:	xor    %eax,%eax
  1f78e8:	mov    0x3f8(%rsp),%rcx
  1f78f0:	xchg   %eax,(%rcx)
  1f78f2:	cmp    $0x2,%eax
  1f78f5:	mov    0x18(%rsp),%rbp
  1f78fa:	mov    0x8(%rsp),%rbx
  1f78ff:	mov    0x88(%rsp),%r15
  1f7907:	jne    1f6a72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4c2>
  1f790d:	mov    0x3f8(%rsp),%rdi
  1f7915:	call   *0x1d23bd(%rip)        # 3c9cd8 <_DYNAMIC+0xaa8>
  1f791b:	jmp    1f6a72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x4c2>
  1f7920:	mov    %rcx,0x120(%rsp)
  1f7928:	mov    0x58(%rsp),%eax
  1f792c:	mov    %al,0x128(%rsp)
  1f7933:	lea    -0x1dcc03(%rip),%rdi        # 1ad37 <anon.9bc0095d74cae2c1d1f2191eec18534c.121.llvm.14577323480697803226+0x52>
  1f793a:	lea    0x1ca55f(%rip),%rcx        # 3c1ea0 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x138>
  1f7941:	lea    0x1ca620(%rip),%r8        # 3c1f68 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0x78>
  1f7948:	lea    0x120(%rsp),%rdx
  1f7950:	mov    $0x2b,%esi
  1f7955:	call   *0x1d1c0d(%rip)        # 3c9568 <_DYNAMIC+0x338>
  1f795b:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7960:	lea    0x1ca631(%rip),%rdx        # 3c1f98 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0xa8>
  1f7967:	jmp    1f799f <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13ef>
  1f7969:	mov    %rbx,0xe0(%rsp)
  1f7971:	lea    0x1ca5a8(%rip),%rdx        # 3c1f20 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0x30>
  1f7978:	jmp    1f799a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13ea>
  1f797a:	mov    0x8(%rsp),%rbx
  1f797f:	mov    %rbx,0xe0(%rsp)
  1f7987:	lea    0x1ca5c2(%rip),%rdx        # 3c1f50 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0x60>
  1f798e:	jmp    1f799a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x13ea>
  1f7990:	mov    %r11,%r13
  1f7993:	lea    0x1ca59e(%rip),%rdx        # 3c1f38 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0x48>
  1f799a:	mov    0x20(%rsp),%r14
  1f799f:	mov    %r13,%rdi
  1f79a2:	mov    0xe0(%rsp),%rsi
  1f79aa:	call   *0x1d1d70(%rip)        # 3c9720 <_DYNAMIC+0x4f0>
  1f79b0:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f79b5:	mov    0x128(%rsp),%rax
  1f79bd:	movzbl 0x130(%rsp),%ecx
  1f79c5:	mov    %rax,0xc0(%rsp)
  1f79cd:	mov    %cl,0xc8(%rsp)
  1f79d4:	lea    -0x1dcca4(%rip),%rdi        # 1ad37 <anon.9bc0095d74cae2c1d1f2191eec18534c.121.llvm.14577323480697803226+0x52>
  1f79db:	lea    0x1ca49e(%rip),%rcx        # 3c1e80 <anon.9bc0095d74cae2c1d1f2191eec18534c.133.llvm.14577323480697803226+0x118>
  1f79e2:	lea    0x1ca597(%rip),%r8        # 3c1f80 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0x90>
  1f79e9:	lea    0xc0(%rsp),%rdx
  1f79f1:	mov    $0x2b,%esi
  1f79f6:	call   *0x1d1b6c(%rip)        # 3c9568 <_DYNAMIC+0x338>
  1f79fc:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7a01:	mov    %r12,%rdi
  1f7a04:	call   *0x1d22ce(%rip)        # 3c9cd8 <_DYNAMIC+0xaa8>
  1f7a0a:	jmp    1f78c0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1310>
  1f7a0f:	mov    $0x4,%edi
  1f7a14:	mov    %r13,%rsi
  1f7a17:	call   *0x1d1aa3(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7a1d:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7a22:	mov    $0x1,%edi
  1f7a27:	mov    0x8(%rsp),%rsi
  1f7a2c:	call   *0x1d1a8e(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7a32:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7a37:	call   *0x1d22a3(%rip)        # 3c9ce0 <_DYNAMIC+0xab0>
  1f7a3d:	test   %al,%al
  1f7a3f:	jne    1f78b1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1301>
  1f7a45:	movb   $0x1,0x4(%r12)
  1f7a4b:	jmp    1f78b1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1301>
  1f7a50:	mov    0x18(%rsp),%rbp
  1f7a55:	mov    0x8(%rsp),%rbx
  1f7a5a:	call   *0x1d2280(%rip)        # 3c9ce0 <_DYNAMIC+0xab0>
  1f7a60:	test   %al,%al
  1f7a62:	jne    1f78e6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1336>
  1f7a68:	mov    0x3f8(%rsp),%rax
  1f7a70:	movb   $0x1,0x4(%rax)
  1f7a74:	jmp    1f78e6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1336>
  1f7a79:	mov    $0x8,%edi
  1f7a7e:	mov    %r15,%rsi
  1f7a81:	call   *0x1d1a39(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7a87:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7a89:	mov    $0x8,%edi
  1f7a8e:	mov    %r15,%rsi
  1f7a91:	call   *0x1d1a29(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7a97:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7a99:	mov    $0x8,%edi
  1f7a9e:	mov    %r15,%rsi
  1f7aa1:	call   *0x1d1a19(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7aa7:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7aa9:	mov    $0x8,%edi
  1f7aae:	mov    %r15,%rsi
  1f7ab1:	call   *0x1d1a09(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7ab7:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7ab9:	mov    $0x8,%edi
  1f7abe:	mov    %r15,%rsi
  1f7ac1:	call   *0x1d19f9(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7ac7:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7ac9:	mov    $0x8,%edi
  1f7ace:	mov    %r15,%rsi
  1f7ad1:	call   *0x1d19e9(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7ad7:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7ad9:	mov    $0x8,%edi
  1f7ade:	mov    %r15,%rsi
  1f7ae1:	call   *0x1d19d9(%rip)        # 3c94c0 <_DYNAMIC+0x290>
  1f7ae7:	jmp    1f7af6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1546>
  1f7ae9:	lea    0x1ca508(%rip),%rdi        # 3c1ff8 <anon.c57ee1c7f25cea858f8d4563fb9f719e.125.llvm.8746345073337986050+0x108>
  1f7af0:	call   *0x1d20ca(%rip)        # 3c9bc0 <_DYNAMIC+0x990>
  1f7af6:	ud2
  1f7af8:	jmp    1f7c08 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1658>
  1f7afd:	jmp    1f7c08 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1658>
  1f7b02:	jmp    1f7bd1 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1621>
  1f7b07:	mov    %rax,%r15
  1f7b0a:	mov    0x18(%rsp),%rbp
  1f7b0f:	mov    0x8(%rsp),%rbx
  1f7b14:	jmp    1f7c82 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16d2>
  1f7b19:	mov    %rax,%r15
  1f7b1c:	mov    0x8(%rsp),%rbx
  1f7b21:	mov    %rbx,%r14
  1f7b24:	mov    0x18(%rsp),%rbp
  1f7b29:	jmp    1f7cbd <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x170d>
  1f7b2e:	mov    %rax,%r15
  1f7b31:	mov    0x18(%rsp),%rbp
  1f7b36:	mov    0x8(%rsp),%rbx
  1f7b3b:	jmp    1f7cd9 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1729>
  1f7b40:	mov    %rax,%r15
  1f7b43:	mov    0x8(%rsp),%rbx
  1f7b48:	test   %rbx,%rbx
  1f7b4b:	jne    1f7d15 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1765>
  1f7b51:	jmp    1f7d72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17c2>
  1f7b56:	mov    %rax,%r15
  1f7b59:	mov    0x8(%rsp),%rbx
  1f7b5e:	test   %rbx,%rbx
  1f7b61:	jne    1f7d2d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x177d>
  1f7b67:	jmp    1f7d72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17c2>
  1f7b6c:	mov    %rax,%r15
  1f7b6f:	mov    0x8(%rsp),%rbx
  1f7b74:	test   %rbx,%rbx
  1f7b77:	jne    1f7d45 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1795>
  1f7b7d:	jmp    1f7d72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17c2>
  1f7b82:	mov    %rax,%r15
  1f7b85:	mov    0x8(%rsp),%rbx
  1f7b8a:	test   %rbx,%rbx
  1f7b8d:	jne    1f7d5d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17ad>
  1f7b93:	jmp    1f7d72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17c2>
  1f7b98:	jmp    1f7c08 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1658>
  1f7b9a:	mov    %rax,%r15
  1f7b9d:	mov    0x18(%rsp),%rbp
  1f7ba2:	mov    0x8(%rsp),%rbx
  1f7ba7:	jmp    1f7c98 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16e8>
  1f7bac:	mov    %rax,%r15
  1f7baf:	mov    0x18(%rsp),%rbp
  1f7bb4:	jmp    1f7cf0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1740>
  1f7bb9:	mov    %rax,%r15
  1f7bbc:	lea    0xc0(%rsp),%rdi
  1f7bc4:	call   1f35e0 <core::ptr::drop_glue::<std::sync::poison::PoisonError<std::sync::poison::mutex::MutexGuard<core::option::Option<u64>>>>>
  1f7bc9:	jmp    1f7bd4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1624>
  1f7bcb:	call   *0x1d18d7(%rip)        # 3c94a8 <_DYNAMIC+0x278>
  1f7bd1:	mov    %rax,%r15
  1f7bd4:	movzbl 0x58(%rsp),%esi
  1f7bd9:	mov    0x3f8(%rsp),%rdi
  1f7be1:	call   1f36c0 <core::ptr::drop_glue::<std::sync::poison::mutex::MutexGuard<std::collections::hash::map::HashMap<u64, (u64, u64, u64)>>>>
  1f7be6:	xor    %r14d,%r14d
  1f7be9:	mov    0x18(%rsp),%rbp
  1f7bee:	mov    0x8(%rsp),%rbx
  1f7bf3:	jmp    1f7c59 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16a9>
  1f7bf5:	call   *0x1d18ad(%rip)        # 3c94a8 <_DYNAMIC+0x278>
  1f7bfb:	jmp    1f7c10 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1660>
  1f7bfd:	mov    %rax,%r15
  1f7c00:	mov    %r14,%rbx
  1f7c03:	xor    %r14d,%r14d
  1f7c06:	jmp    1f7c59 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16a9>
  1f7c08:	mov    %rax,%r15
  1f7c0b:	xor    %r14d,%r14d
  1f7c0e:	jmp    1f7c59 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16a9>
  1f7c10:	mov    %rax,%r15
  1f7c13:	jmp    1f7c59 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16a9>
  1f7c15:	mov    %rax,%r15
  1f7c18:	mov    0x120(%rsp),%r12
  1f7c20:	cmpb   $0x0,0x128(%rsp)
  1f7c28:	jne    1f7c47 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1697>
  1f7c2a:	mov    0x1d209f(%rip),%rax        # 3c9cd0 <_DYNAMIC+0xaa0>
  1f7c31:	mov    (%rax),%rax
  1f7c34:	movabs $0x7fffffffffffffff,%rcx
  1f7c3e:	test   %rcx,%rax
  1f7c41:	jne    1f7d88 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17d8>
  1f7c47:	xor    %r14d,%r14d
  1f7c4a:	xor    %eax,%eax
  1f7c4c:	xchg   %eax,(%r12)
  1f7c50:	cmp    $0x2,%eax
  1f7c53:	je     1f7d7a <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17ca>
  1f7c59:	mov    0x168(%rsp),%rsi
  1f7c61:	test   %rsi,%rsi
  1f7c64:	je     1f7c7d <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16cd>
  1f7c66:	mov    0x170(%rsp),%rdi
  1f7c6e:	shl    $0x3,%rsi
  1f7c72:	mov    $0x8,%edx
  1f7c77:	call   *0x1d1823(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7c7d:	test   %r14b,%r14b
  1f7c80:	jne    1f7c98 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x16e8>
  1f7c82:	mov    $0x1,%edx
  1f7c87:	mov    0xa0(%rsp),%rdi
  1f7c8f:	mov    %rbx,%rsi
  1f7c92:	call   *0x1d1808(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7c98:	mov    0xa8(%rsp),%r14
  1f7ca0:	test   %r14,%r14
  1f7ca3:	je     1f7cd4 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1724>
  1f7ca5:	lea    0x0(,%r14,8),%rsi
  1f7cad:	mov    $0x8,%edx
  1f7cb2:	mov    0x78(%rsp),%rdi
  1f7cb7:	call   *0x1d17e3(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7cbd:	shl    $0x3,%r14
  1f7cc1:	mov    $0x8,%edx
  1f7cc6:	mov    0x60(%rsp),%rdi
  1f7ccb:	mov    %r14,%rsi
  1f7cce:	call   *0x1d17cc(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7cd4:	test   %rbx,%rbx
  1f7cd7:	je     1f7cf0 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1740>
  1f7cd9:	shl    $0x2,%rbx
  1f7cdd:	mov    $0x4,%edx
  1f7ce2:	mov    0x48(%rsp),%rdi
  1f7ce7:	mov    %rbx,%rsi
  1f7cea:	call   *0x1d17b0(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7cf0:	mov    0x108(%rsp),%rbx
  1f7cf8:	test   %rbx,%rbx
  1f7cfb:	je     1f7d72 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x17c2>
  1f7cfd:	lea    0x0(,%rbx,8),%rsi
  1f7d05:	mov    $0x8,%edx
  1f7d0a:	mov    0x28(%rsp),%rdi
  1f7d0f:	call   *0x1d178b(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7d15:	lea    0x0(,%rbx,8),%rsi
  1f7d1d:	mov    $0x8,%edx
  1f7d22:	mov    0x38(%rsp),%rdi
  1f7d27:	call   *0x1d1773(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7d2d:	lea    0x0(,%rbx,8),%rsi
  1f7d35:	mov    $0x8,%edx
  1f7d3a:	mov    0x30(%rsp),%rdi
  1f7d3f:	call   *0x1d175b(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7d45:	lea    0x0(,%rbx,8),%rsi
  1f7d4d:	mov    $0x8,%edx
  1f7d52:	mov    0x40(%rsp),%rdi
  1f7d57:	call   *0x1d1743(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7d5d:	shl    $0x3,%rbx
  1f7d61:	mov    $0x8,%edx
  1f7d66:	mov    %rbp,%rdi
  1f7d69:	mov    %rbx,%rsi
  1f7d6c:	call   *0x1d172e(%rip)        # 3c94a0 <_DYNAMIC+0x270>
  1f7d72:	mov    %r15,%rdi
  1f7d75:	call   3bf0f0 <_Unwind_Resume@plt>
  1f7d7a:	mov    %r12,%rdi
  1f7d7d:	call   *0x1d1f55(%rip)        # 3c9cd8 <_DYNAMIC+0xaa8>
  1f7d83:	jmp    1f7be6 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1636>
  1f7d88:	call   *0x1d1f52(%rip)        # 3c9ce0 <_DYNAMIC+0xab0>
  1f7d8e:	test   %al,%al
  1f7d90:	mov    0x18(%rsp),%rbp
  1f7d95:	mov    0x8(%rsp),%rbx
  1f7d9a:	jne    1f7c47 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1697>
  1f7da0:	movb   $0x1,0x4(%r12)
  1f7da6:	jmp    1f7c47 <crypto_lib::cryptanalysis::prime_fast::rho_worker+0x1697>
  1f7dab:	call   *0x1d16f7(%rip)        # 3c94a8 <_DYNAMIC+0x278>
