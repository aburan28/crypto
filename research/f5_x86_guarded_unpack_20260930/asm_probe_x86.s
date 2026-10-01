	.file	"f5_unpack_asm.edc564a347f51d4f-cgu.0"
	.section	.text.checked_decode,"ax",@progbits
	.globl	checked_decode
	.p2align	4
	.type	checked_decode,@function
checked_decode:
.Lfunc_begin0:
	.cfi_startproc
	.cfi_personality 155, DW.ref.rust_eh_personality
	.cfi_lsda 27, .Lexception0
	pushq	%r15
	.cfi_def_cfa_offset 16
	pushq	%r14
	.cfi_def_cfa_offset 24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	.cfi_offset %r15, -16
	testq	%rsi, %rsi
	je	.LBB0_1
	leaq	(%rdi,%rsi,8), %r9
	xorl	%r11d, %r11d
	xorl	%r10d, %r10d
	jmp	.LBB0_4
	.p2align	4
.LBB0_5:
	movq	%r11, %rax
.LBB0_9:
	addq	$8, %rdi
	incq	%r10
	movq	%rax, %r11
	cmpq	%r9, %rdi
	je	.LBB0_10
.LBB0_4:
	movq	(%rdi), %rbx
	testq	%rbx, %rbx
	je	.LBB0_5
	movq	%r10, %r14
	shlq	$6, %r14
	.p2align	4
.LBB0_7:
	rep		bsfq	%rbx, %rsi
	orq	%r14, %rsi
	cmpq	%rcx, %rsi
	jae	.LBB0_11
	leaq	-1(%rbx), %r15
	leaq	1(%r11), %rax
	movq	(%rdx,%rsi,8), %rsi
	movq	%rsi, (%r8,%r11,8)
	movq	%rax, %r11
	andq	%r15, %rbx
	jne	.LBB0_7
	jmp	.LBB0_9
.LBB0_1:
	xorl	%eax, %eax
.LBB0_10:
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	retq
.LBB0_11:
	.cfi_def_cfa_offset 32
.Ltmp0:
	leaq	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.1(%rip), %rdx
	movq	%rsi, %rdi
	movq	%rcx, %rsi
	callq	*_RNvNtCsc36rpYXAlPq_4core9panicking18panic_bounds_check@GOTPCREL(%rip)
.Ltmp1:
	ud2
.LBB0_2:
.Ltmp2:
	callq	*_RNvNtCsc36rpYXAlPq_4core9panicking19panic_cannot_unwind@GOTPCREL(%rip)
.Lfunc_end0:
	.size	checked_decode, .Lfunc_end0-checked_decode
	.cfi_endproc
	.section	.gcc_except_table.checked_decode,"a",@progbits
	.p2align	2, 0x0
GCC_except_table0:
.Lexception0:
	.byte	255
	.byte	155
	.uleb128 .Lttbase0-.Lttbaseref0
.Lttbaseref0:
	.byte	1
	.uleb128 .Lcst_end0-.Lcst_begin0
.Lcst_begin0:
	.uleb128 .Ltmp0-.Lfunc_begin0
	.uleb128 .Ltmp1-.Ltmp0
	.uleb128 .Ltmp2-.Lfunc_begin0
	.byte	1
.Lcst_end0:
	.byte	127
	.byte	0
	.p2align	2, 0x0
.Lttbase0:
	.byte	0
	.p2align	2, 0x0

	.section	.text.guarded_decode,"ax",@progbits
	.globl	guarded_decode
	.p2align	4
	.type	guarded_decode,@function
guarded_decode:
.Lfunc_begin1:
	.cfi_startproc
	.cfi_personality 155, DW.ref.rust_eh_personality
	.cfi_lsda 27, .Lexception1
	pushq	%r15
	.cfi_def_cfa_offset 16
	pushq	%r14
	.cfi_def_cfa_offset 24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	.cfi_offset %r15, -16
	movq	%rcx, %r9
	movq	%rcx, %r10
	shrq	$6, %r10
	xorl	%eax, %eax
	andq	$63, %rcx
	setne	%al
	addq	%r10, %rax
	testq	%rsi, %rsi
	je	.LBB1_3
	testq	%rcx, %rcx
	je	.LBB1_3
	cmpq	%rax, %rsi
	jne	.LBB1_3
	movq	-8(%rdi,%rsi,8), %rax
	shrq	%cl, %rax
	testq	%rax, %rax
	jne	.LBB1_8
	jmp	.LBB1_6
.LBB1_3:
	cmpq	%rax, %rsi
	jbe	.LBB1_4
.LBB1_8:
	leaq	(%rdi,%rsi,8), %rsi
	xorl	%r11d, %r11d
	xorl	%r10d, %r10d
	jmp	.LBB1_9
	.p2align	4
.LBB1_10:
	movq	%r11, %rax
.LBB1_14:
	addq	$8, %rdi
	incq	%r10
	movq	%rax, %r11
	cmpq	%rsi, %rdi
	je	.LBB1_15
.LBB1_9:
	movq	(%rdi), %rbx
	testq	%rbx, %rbx
	je	.LBB1_10
	movq	%r10, %r14
	shlq	$6, %r14
	.p2align	4
.LBB1_12:
	rep		bsfq	%rbx, %rcx
	orq	%r14, %rcx
	cmpq	%r9, %rcx
	jae	.LBB1_16
	leaq	-1(%rbx), %r15
	leaq	1(%r11), %rax
	movq	(%rdx,%rcx,8), %rcx
	movq	%rcx, (%r8,%r11,8)
	movq	%rax, %r11
	andq	%r15, %rbx
	jne	.LBB1_12
	jmp	.LBB1_14
.LBB1_4:
	testq	%rsi, %rsi
	je	.LBB1_5
.LBB1_6:
	leaq	(%rdi,%rsi,8), %rcx
	xorl	%eax, %eax
	xorl	%esi, %esi
	jmp	.LBB1_19
	.p2align	4
.LBB1_18:
	addq	$8, %rdi
	incq	%rsi
	cmpq	%rcx, %rdi
	je	.LBB1_15
.LBB1_19:
	movq	(%rdi), %r10
	testq	%r10, %r10
	je	.LBB1_18
	movq	%rsi, %r9
	shlq	$6, %r9
	.p2align	4
.LBB1_21:
	movq	%r10, %r11
	rep		bsfq	%r10, %rbx
	orq	%r9, %rbx
	decq	%r10
	movq	(%rdx,%rbx,8), %rbx
	movq	%rbx, (%r8,%rax,8)
	incq	%rax
	andq	%r11, %r10
	jne	.LBB1_21
	jmp	.LBB1_18
.LBB1_5:
	xorl	%eax, %eax
.LBB1_15:
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	retq
.LBB1_16:
	.cfi_def_cfa_offset 32
.Ltmp3:
	leaq	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.1(%rip), %rdx
	movq	%rcx, %rdi
	movq	%r9, %rsi
	callq	*_RNvNtCsc36rpYXAlPq_4core9panicking18panic_bounds_check@GOTPCREL(%rip)
.Ltmp4:
	ud2
.LBB1_22:
.Ltmp5:
	callq	*_RNvNtCsc36rpYXAlPq_4core9panicking19panic_cannot_unwind@GOTPCREL(%rip)
.Lfunc_end1:
	.size	guarded_decode, .Lfunc_end1-guarded_decode
	.cfi_endproc
	.section	.gcc_except_table.guarded_decode,"a",@progbits
	.p2align	2, 0x0
GCC_except_table1:
.Lexception1:
	.byte	255
	.byte	155
	.uleb128 .Lttbase1-.Lttbaseref1
.Lttbaseref1:
	.byte	1
	.uleb128 .Lcst_end1-.Lcst_begin1
.Lcst_begin1:
	.uleb128 .Ltmp3-.Lfunc_begin1
	.uleb128 .Ltmp4-.Ltmp3
	.uleb128 .Ltmp5-.Lfunc_begin1
	.byte	1
.Lcst_end1:
	.byte	127
	.byte	0
	.p2align	2, 0x0
.Lttbase1:
	.byte	0
	.p2align	2, 0x0

	.type	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.0,@object
	.section	.rodata.str1.1,"aMS",@progbits,1
.Lanon.51a5b7b42f3970c5a338e42e041aa36c.0:
	.asciz	"/private/tmp/f5-unpack-asm.rs"
	.size	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.0, 30

	.type	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.1,@object
	.section	.data.rel.ro..Lanon.51a5b7b42f3970c5a338e42e041aa36c.1,"aw",@progbits
	.p2align	3, 0x0
.Lanon.51a5b7b42f3970c5a338e42e041aa36c.1:
	.quad	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.0
	.asciz	"\035\000\000\000\000\000\000\000\024\000\000\0002\000\000"
	.size	.Lanon.51a5b7b42f3970c5a338e42e041aa36c.1, 24

	.hidden	DW.ref.rust_eh_personality
	.weak	DW.ref.rust_eh_personality
	.section	.data.DW.ref.rust_eh_personality,"awG",@progbits,DW.ref.rust_eh_personality,comdat
	.p2align	3, 0x0
	.type	DW.ref.rust_eh_personality,@object
	.size	DW.ref.rust_eh_personality, 8
DW.ref.rust_eh_personality:
	.quad	rust_eh_personality
	.ident	"rustc version 1.98.0 (88d9e12ae 2026-08-18)"
	.section	".note.GNU-stack","",@progbits
