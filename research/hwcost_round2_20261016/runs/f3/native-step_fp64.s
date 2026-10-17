000000000022db30 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64>:
  22db30:	push   %rbp
  22db31:	push   %r15
  22db33:	push   %r14
  22db35:	push   %r13
  22db37:	push   %r12
  22db39:	push   %rbx
  22db3a:	sub    $0x38,%rsp
  22db3e:	mov    %rcx,0x30(%rsp)
  22db43:	mov    %rdx,%r15
  22db46:	mov    %rsi,%rbp
  22db49:	mov    %rdi,%r12
  22db4c:	mov    0x10(%rsi),%r13
  22db50:	test   %r13,%r13
  22db53:	je     22dba0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x70>
  22db55:	mov    0x40(%rbp),%rdi
  22db59:	mov    0x8(%rbp),%rax
  22db5d:	mov    0x38(%rbp),%rcx
  22db61:	mov    0xf0(%r12),%r8
  22db69:	lea    -0x1(%r13),%rdx
  22db6d:	cmp    %rdx,%rdi
  22db70:	cmovb  %rdi,%rdx
  22db74:	inc    %rdx
  22db77:	cmp    $0x11,%rdx
  22db7b:	setb   %sil
  22db7f:	mov    %rcx,%r9
  22db82:	sub    %rax,%r9
  22db85:	cmp    $0x80,%r9
  22db8c:	setb   %r9b
  22db90:	or     %sil,%r9b
  22db93:	cmp    $0x1,%r9b
  22db97:	jne    22dbad <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x7d>
  22db99:	xor    %edx,%edx
  22db9b:	jmp    22dcc0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x190>
  22dba0:	mov    0x100(%r12),%rsi
  22dba8:	jmp    22dd81 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x251>
  22dbad:	mov    %edx,%esi
  22dbaf:	and    $0xf,%esi
  22dbb2:	mov    $0x10,%r9d
  22dbb8:	cmovne %rsi,%r9
  22dbbc:	sub    %r9,%rdx
  22dbbf:	vpbroadcastq %r8,%ymm0
  22dbc5:	xor    %esi,%esi
  22dbc7:	vpbroadcastq -0x21c220(%rip),%ymm1        # 119b0 <anon.b68ad25efc32a4f1d321fe93ba3d9fb4.405.llvm.17162556477609682337+0x10>
  22dbd0:	vmovdqu (%rax,%rsi,8),%ymm2
  22dbd5:	vpsrlq $0x1,%ymm2,%ymm3
  22dbda:	vpand  %ymm1,%ymm3,%ymm3
  22dbde:	kxnorw %k0,%k0,%k1
  22dbe2:	vpxor  %xmm4,%xmm4,%xmm4
  22dbe6:	vpgatherqq (%r15,%ymm3,8),%ymm4{%k1}
  22dbed:	vmovdqu 0x20(%rax,%rsi,8),%ymm3
  22dbf3:	vpsrlq $0x1,%ymm3,%ymm5
  22dbf8:	vpand  %ymm1,%ymm5,%ymm5
  22dbfc:	kxnorw %k0,%k0,%k1
  22dc00:	vpxor  %xmm6,%xmm6,%xmm6
  22dc04:	vpgatherqq (%r15,%ymm5,8),%ymm6{%k1}
  22dc0b:	vmovdqu 0x40(%rax,%rsi,8),%ymm5
  22dc11:	vpsrlq $0x1,%ymm5,%ymm7
  22dc16:	vpand  %ymm1,%ymm7,%ymm7
  22dc1a:	kxnorw %k0,%k0,%k1
  22dc1e:	vpxor  %xmm8,%xmm8,%xmm8
  22dc23:	vpgatherqq (%r15,%ymm7,8),%ymm8{%k1}
  22dc2a:	vmovdqu 0x60(%rax,%rsi,8),%ymm7
  22dc30:	vpsrlq $0x1,%ymm7,%ymm9
  22dc35:	vpand  %ymm1,%ymm9,%ymm9
  22dc39:	kxnorw %k0,%k0,%k1
  22dc3d:	vpxor  %xmm10,%xmm10,%xmm10
  22dc42:	vpgatherqq (%r15,%ymm9,8),%ymm10{%k1}
  22dc49:	vpcmpltuq %ymm2,%ymm4,%k1
  22dc50:	vpcmpltuq %ymm3,%ymm6,%k2
  22dc57:	vpcmpltuq %ymm5,%ymm8,%k3
  22dc5e:	vpcmpltuq %ymm7,%ymm10,%k4
  22dc65:	vpsubq %ymm2,%ymm4,%ymm2
  22dc69:	vpsubq %ymm3,%ymm6,%ymm3
  22dc6d:	vpsubq %ymm5,%ymm8,%ymm4
  22dc71:	vpsubq %ymm7,%ymm10,%ymm5
  22dc75:	vpaddq %ymm0,%ymm2,%ymm2{%k1}
  22dc7b:	vpaddq %ymm0,%ymm3,%ymm3{%k2}
  22dc81:	vpaddq %ymm0,%ymm4,%ymm4{%k3}
  22dc87:	vpaddq %ymm0,%ymm5,%ymm5{%k4}
  22dc8d:	vmovdqu %ymm2,(%rcx,%rsi,8)
  22dc92:	vmovdqu %ymm3,0x20(%rcx,%rsi,8)
  22dc98:	vmovdqu %ymm4,0x40(%rcx,%rsi,8)
  22dc9e:	vmovdqu %ymm5,0x60(%rcx,%rsi,8)
  22dca4:	add    $0x10,%rsi
  22dca8:	cmp    %rsi,%rdx
  22dcab:	jne    22dbd0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0xa0>
  22dcb1:	data16 data16 data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  22dcc0:	cmp    %rdx,%rdi
  22dcc3:	je     22e009 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x4d9>
  22dcc9:	mov    (%rax,%rdx,8),%rsi
  22dccd:	mov    %esi,%r9d
  22dcd0:	and    $0x3e,%r9d
  22dcd4:	mov    (%r15,%r9,4),%r9
  22dcd8:	sub    %rsi,%r9
  22dcdb:	mov    $0x0,%esi
  22dce0:	cmovb  %r8,%rsi
  22dce4:	add    %r9,%rsi
  22dce7:	mov    %rsi,(%rcx,%rdx,8)
  22dceb:	inc    %rdx
  22dcee:	cmp    %rdx,%r13
  22dcf1:	jne    22dcc0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x190>
  22dcf3:	mov    0xf8(%r12),%r9
  22dcfb:	mov    0x100(%r12),%rsi
  22dd03:	mov    0x50(%rbp),%r10
  22dd07:	mov    0x58(%rbp),%rax
  22dd0b:	mov    0x60(%rbp),%r11
  22dd0f:	xor    %ebx,%ebx
  22dd11:	jmp    22dd53 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x223>
  22dd13:	data16 data16 data16 cs nopw 0x0(%rax,%rax,1)
  22dd20:	mulx   %rsi,%r14,%rsi
  22dd25:	mov    %r9,%rdx
  22dd28:	imul   %r14,%rdx
  22dd2c:	mulx   %r8,%rdx,%rdx
  22dd31:	cmp    $0x1,%r14
  22dd35:	sbb    $0xffffffffffffffff,%rsi
  22dd39:	add    %rdx,%rsi
  22dd3c:	cmp    %r8,%rsi
  22dd3f:	mov    $0x0,%edx
  22dd44:	cmovae %r8,%rdx
  22dd48:	sub    %rdx,%rsi
  22dd4b:	inc    %rbx
  22dd4e:	cmp    %rbx,%r13
  22dd51:	je     22dd81 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x251>
  22dd53:	cmp    %rbx,%rax
  22dd56:	je     22e01c <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x4ec>
  22dd5c:	mov    %rsi,(%r10,%rbx,8)
  22dd60:	cmp    %rbx,%rdi
  22dd63:	je     22e032 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x502>
  22dd69:	mov    (%rcx,%rbx,8),%rdx
  22dd6d:	test   %rdx,%rdx
  22dd70:	jne    22dd20 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x1f0>
  22dd72:	inc    %r11
  22dd75:	mov    %r11,0x60(%rbp)
  22dd79:	inc    %rbx
  22dd7c:	cmp    %rbx,%r13
  22dd7f:	jne    22dd53 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x223>
  22dd81:	mov    %r12,%rdi
  22dd84:	vzeroupper
  22dd87:	call   *0x1873ab(%rip)        # 3b5138 <_DYNAMIC+0xcd0>
  22dd8d:	test   $0x1,%al
  22dd8f:	je     22dff6 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x4c6>
  22dd95:	mov    0x38(%rbp),%rax
  22dd99:	mov    %rax,0x10(%rsp)
  22dd9e:	mov    0x40(%rbp),%rcx
  22dda2:	mov    0xf0(%r12),%r10
  22ddaa:	mov    0xf8(%r12),%r11
  22ddb2:	lea    -0x1(%r13),%rax
  22ddb6:	mov    %r13,%rdi
  22ddb9:	mov    %rcx,0x8(%rsp)
  22ddbe:	cmp    %rcx,%rax
  22ddc1:	jae    22de78 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x348>
  22ddc7:	mov    0x58(%rbp),%rax
  22ddcb:	mov    0x50(%rbp),%rbx
  22ddcf:	xor    %r14d,%r14d
  22ddd2:	mov    %r13,%rdi
  22ddd5:	data16 cs nopw 0x0(%rax,%rax,1)
  22dde0:	test   %rdi,%rdi
  22dde3:	je     22de81 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x351>
  22dde9:	mov    0x10(%rsp),%rcx
  22ddee:	mov    -0x8(%rcx,%rdi,8),%rcx
  22ddf3:	dec    %rdi
  22ddf6:	test   %rcx,%rcx
  22ddf9:	je     22dde0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x2b0>
  22ddfb:	cmp    %rax,%rdi
  22ddfe:	jae    22e084 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x554>
  22de04:	mov    %rdx,%r8
  22de07:	mulx   (%rbx,%rdi,8),%rsi,%r9
  22de0d:	mov    %r11,%rdx
  22de10:	imul   %rsi,%rdx
  22de14:	mulx   %r10,%rdx,%rdx
  22de19:	cmp    $0x1,%rsi
  22de1d:	sbb    $0xffffffffffffffff,%r9
  22de21:	add    %rdx,%r9
  22de24:	cmp    %r10,%r9
  22de27:	mov    %r10,%rdx
  22de2a:	cmovb  %r14,%rdx
  22de2e:	sub    %rdx,%r9
  22de31:	mov    0x10(%rsp),%rdx
  22de36:	mov    %r9,(%rdx,%rdi,8)
  22de3a:	mov    %rcx,%rdx
  22de3d:	mulx   %r8,%rcx,%r8
  22de42:	mov    %r11,%rdx
  22de45:	imul   %rcx,%rdx
  22de49:	mulx   %r10,%rsi,%rsi
  22de4e:	mov    %r8,%rdx
  22de51:	cmp    $0x1,%rcx
  22de55:	sbb    $0xffffffffffffffff,%rdx
  22de59:	add    %rsi,%rdx
  22de5c:	cmp    %r10,%rdx
  22de5f:	mov    %r10,%rcx
  22de62:	cmovb  %r14,%rcx
  22de66:	sub    %rcx,%rdx
  22de69:	lea    -0x1(%rdi),%rcx
  22de6d:	cmp    0x8(%rsp),%rcx
  22de72:	jb     22dde0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x2b0>
  22de78:	test   %rdi,%rdi
  22de7b:	jne    22e06f <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x53f>
  22de81:	mov    %rbp,0x18(%rsp)
  22de86:	test   %r13,%r13
  22de89:	je     22dffa <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x4ca>
  22de8f:	mov    0x18(%rsp),%rax
  22de94:	mov    0x8(%rax),%rcx
  22de98:	mov    0x20(%rax),%r8
  22de9c:	mov    0x28(%rax),%rdx
  22dea0:	mov    %rdx,0x28(%rsp)
  22dea5:	mov    0x68(%rax),%rax
  22dea9:	mov    %rax,0x20(%rsp)
  22deae:	xor    %ebp,%ebp
  22deb0:	xor    %ebx,%ebx
  22deb2:	jmp    22decc <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x39c>
  22deb4:	data16 data16 cs nopw 0x0(%rax,%rax,1)
  22dec0:	inc    %rbx
  22dec3:	cmp    %rbx,%r13
  22dec6:	je     22dffa <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x4ca>
  22decc:	cmp    %rbx,0x28(%rsp)
  22ded1:	je     22e045 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x515>
  22ded7:	cmp    %rbx,0x8(%rsp)
  22dedc:	je     22e05a <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x52a>
  22dee2:	mov    (%rcx,%rbx,8),%r9
  22dee6:	mov    %r9d,%esi
  22dee9:	and    $0x3e,%esi
  22deec:	mov    0x100(%r15,%rsi,4),%rax
  22def4:	mov    (%r8,%rbx,8),%r14
  22def8:	sub    %r14,%rax
  22defb:	mov    $0x0,%edx
  22df00:	cmovb  %r10,%rdx
  22df04:	add    %rax,%rdx
  22df07:	mov    0x10(%rsp),%rax
  22df0c:	mulx   (%rax,%rbx,8),%r12,%rax
  22df12:	mov    %r11,%rdx
  22df15:	imul   %r12,%rdx
  22df19:	mulx   %r10,%rdx,%rdx
  22df1e:	cmp    $0x1,%r12
  22df22:	sbb    $0xffffffffffffffff,%rax
  22df26:	add    %rdx,%rax
  22df29:	cmp    %r10,%rax
  22df2c:	mov    %r10,%rdx
  22df2f:	cmovb  %rbp,%rdx
  22df33:	sub    %rdx,%rax
  22df36:	mov    %rax,%rdx
  22df39:	mulx   %rax,%r12,%rdi
  22df3e:	mov    %r11,%rdx
  22df41:	imul   %r12,%rdx
  22df45:	cmp    $0x1,%r12
  22df49:	mulx   %r10,%rdx,%rdx
  22df4e:	sbb    $0xffffffffffffffff,%rdi
  22df52:	add    %rdx,%rdi
  22df55:	cmp    %r10,%rdi
  22df58:	mov    %r10,%rdx
  22df5b:	cmovb  %rbp,%rdx
  22df5f:	sub    %rdx,%rdi
  22df62:	sub    %r9,%rdi
  22df65:	mov    $0x0,%edx
  22df6a:	cmovb  %r10,%rdx
  22df6e:	add    %rdi,%rdx
  22df71:	sub    (%r15,%rsi,4),%rdx
  22df75:	mov    $0x0,%esi
  22df7a:	cmovb  %r10,%rsi
  22df7e:	add    %rdx,%rsi
  22df81:	sub    %rsi,%r9
  22df84:	mov    $0x0,%edx
  22df89:	cmovb  %r10,%rdx
  22df8d:	add    %r9,%rdx
  22df90:	mulx   %rax,%rax,%rdi
  22df95:	mov    %r11,%rdx
  22df98:	imul   %rax,%rdx
  22df9c:	cmp    $0x1,%rax
  22dfa0:	mulx   %r10,%rax,%rax
  22dfa5:	sbb    $0xffffffffffffffff,%rdi
  22dfa9:	add    %rax,%rdi
  22dfac:	cmp    %r10,%rdi
  22dfaf:	mov    %r10,%rax
  22dfb2:	cmovb  %rbp,%rax
  22dfb6:	sub    %rax,%rdi
  22dfb9:	sub    %r14,%rdi
  22dfbc:	mov    $0x0,%eax
  22dfc1:	cmovb  %r10,%rax
  22dfc5:	add    %rdi,%rax
  22dfc8:	mov    %rsi,(%rcx,%rbx,8)
  22dfcc:	mov    %rax,(%r8,%rbx,8)
  22dfd0:	test   %rsi,0x30(%rsp)
  22dfd5:	jne    22dec0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x390>
  22dfdb:	mov    0x20(%rsp),%rdx
  22dfe0:	inc    %rdx
  22dfe3:	mov    0x18(%rsp),%rax
  22dfe8:	mov    %rdx,0x20(%rsp)
  22dfed:	mov    %rdx,0x68(%rax)
  22dff1:	jmp    22dec0 <crypto_lib::cryptanalysis::hwcost::kernel::step_fp64+0x390>
  22dff6:	incq   0x60(%rbp)
  22dffa:	add    $0x38,%rsp
  22dffe:	pop    %rbx
  22dfff:	pop    %r12
  22e001:	pop    %r13
  22e003:	pop    %r14
  22e005:	pop    %r15
  22e007:	pop    %rbp
  22e008:	ret
  22e009:	lea    0x17fbf8(%rip),%rdx        # 3adc08 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x3a0>
  22e010:	mov    %rdi,%rsi
  22e013:	vzeroupper
  22e016:	call   *0x18694c(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e01c:	lea    0x17fbb5(%rip),%rdx        # 3adbd8 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x370>
  22e023:	mov    %rax,%rdi
  22e026:	mov    %rax,%rsi
  22e029:	vzeroupper
  22e02c:	call   *0x186936(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e032:	lea    0x17fbb7(%rip),%rdx        # 3adbf0 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x388>
  22e039:	mov    %rdi,%rsi
  22e03c:	vzeroupper
  22e03f:	call   *0x186923(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e045:	lea    0x17fb14(%rip),%rdx        # 3adb60 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x2f8>
  22e04c:	mov    0x28(%rsp),%rdi
  22e051:	mov    %rdi,%rsi
  22e054:	call   *0x18690e(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e05a:	lea    0x17fb17(%rip),%rdx        # 3adb78 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x310>
  22e061:	mov    0x8(%rsp),%rdi
  22e066:	mov    %rdi,%rsi
  22e069:	call   *0x1868f9(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e06f:	dec    %rdi
  22e072:	lea    0x17fb2f(%rip),%rdx        # 3adba8 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x340>
  22e079:	mov    0x8(%rsp),%rsi
  22e07e:	call   *0x1868e4(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e084:	lea    0x17fb35(%rip),%rdx        # 3adbc0 <anon.db6fdc46c529212cfce610ad9ecec8c4.796.llvm.7120734616223197918+0x358>
  22e08b:	mov    %rax,%rsi
  22e08e:	call   *0x1868d4(%rip)        # 3b4968 <_DYNAMIC+0x500>
  22e094:	int3
  22e095:	int3
  22e096:	int3
  22e097:	int3
  22e098:	int3
  22e099:	int3
  22e09a:	int3
  22e09b:	int3
  22e09c:	int3
  22e09d:	int3
  22e09e:	int3
  22e09f:	int3

