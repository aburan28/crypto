	.build_version macos, 11, 0
	.section	__TEXT,__text,regular,pure_instructions
	.p2align	2
__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E:
	.cfi_startproc
	sub	sp, sp, #160
	.cfi_def_cfa_offset 160
	stp	x28, x27, [sp, #64]
	stp	x26, x25, [sp, #80]
	stp	x24, x23, [sp, #96]
	stp	x22, x21, [sp, #112]
	stp	x20, x19, [sp, #128]
	stp	x29, x30, [sp, #144]
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	mov	x27, x0
	ldp	x16, x15, [x1]
	ldp	x13, x11, [x1, #16]
	ldr	x10, [x1, #32]
	ldp	x12, x14, [x2]
	ldp	x17, x1, [x2, #16]
	ldr	x2, [x2, #32]
	umulh	x3, x12, x16
	mul	x9, x12, x16
	umulh	x8, x2, x15
	mul	x4, x2, x15
	umulh	x5, x1, x13
	mul	x6, x1, x13
	umulh	x7, x17, x11
	mul	x19, x17, x11
	umulh	x20, x14, x10
	mul	x21, x14, x10
	adds	x19, x19, x21
	adc	x7, x7, x20
	adds	x6, x19, x6
	adc	x5, x7, x5
	adds	x4, x6, x4
	adc	x5, x5, x8
	mov	w8, #19
	umulh	x6, x4, x8
	madd	x5, x5, x8, x6
	mul	x4, x4, x8
	adds	x9, x4, x9
	str	x9, [sp, #56]
	adc	x3, x5, x3
	umulh	x4, x14, x16
	mul	x5, x14, x16
	umulh	x6, x12, x15
	mul	x7, x12, x15
	adds	x5, x5, x7
	adc	x4, x4, x6
	umulh	x6, x2, x13
	mul	x7, x2, x13
	umulh	x19, x1, x11
	mul	x20, x1, x11
	umulh	x21, x17, x10
	mul	x22, x17, x10
	adds	x20, x20, x22
	adc	x19, x19, x21
	adds	x7, x20, x7
	adc	x6, x19, x6
	umulh	x19, x7, x8
	madd	x6, x6, x8, x19
	mul	x7, x7, x8
	adds	x5, x5, x7
	adc	x6, x4, x6
	umulh	x9, x17, x16
	str	x9, [sp, #16]
	mul	x19, x17, x16
	umulh	x20, x14, x15
	mul	x21, x14, x15
	umulh	x22, x12, x13
	mul	x23, x12, x13
	umulh	x4, x2, x11
	mul	x24, x2, x11
	umulh	x25, x1, x10
	mul	x26, x1, x10
	adds	x24, x24, x26
	adc	x4, x4, x25
	umulh	x25, x24, x8
	madd	x25, x4, x8, x25
	umulh	x9, x1, x16
	str	x9, [sp, #24]
	mul	x7, x1, x16
	umulh	x28, x17, x15
	mul	x30, x17, x15
	mul	x9, x10, x2
	umulh	x4, x10, x2
	umulh	x26, x9, x8
	madd	x26, x4, x8, x26
	umulh	x0, x2, x16
	mul	x16, x2, x16
	stp	x16, x0, [sp, #40]
	umulh	x16, x1, x15
	str	x16, [sp, #32]
	mul	x15, x1, x15
	str	x15, [sp, #8]
	umulh	x4, x17, x13
	mul	x2, x17, x13
	umulh	x16, x14, x13
	mul	x15, x14, x13
	umulh	x1, x14, x11
	mul	x0, x14, x11
	umulh	x13, x12, x11
	mul	x11, x12, x11
	umulh	x17, x12, x10
	mul	x12, x12, x10
	ldr	x14, [sp, #56]
	extr	x10, x3, x14, #51
	adds	x10, x5, x10
	lsr	x3, x3, #51
	adc	x3, x6, x3
	adds	x5, x21, x23
	adc	x6, x20, x22
	adds	x5, x5, x19
	ldr	x19, [sp, #16]
	adc	x6, x6, x19
	mul	x19, x24, x8
	adds	x5, x5, x19
	adc	x6, x6, x25
	extr	x19, x3, x10, #51
	adds	x5, x5, x19
	lsr	x3, x3, #51
	adc	x3, x6, x3
	adds	x11, x15, x11
	adc	x13, x16, x13
	adds	x11, x11, x30
	adc	x13, x13, x28
	adds	x11, x11, x7
	ldr	x15, [sp, #24]
	adc	x13, x13, x15
	mul	x9, x9, x8
	adds	x9, x11, x9
	adc	x11, x13, x26
	extr	x13, x3, x5, #51
	adds	x9, x9, x13
	lsr	x13, x3, #51
	adc	x11, x11, x13
	adds	x12, x0, x12
	adc	x13, x1, x17
	adds	x12, x12, x2
	adc	x13, x13, x4
	ldr	x15, [sp, #8]
	adds	x12, x12, x15
	ldr	x15, [sp, #32]
	adc	x13, x13, x15
	ldr	x15, [sp, #40]
	adds	x12, x12, x15
	and	x14, x14, #0x7ffffffffffff
	ldr	x15, [sp, #48]
	adc	x13, x13, x15
	lsr	x15, x11, #51
	extr	x11, x11, x9, #51
	adds	x11, x12, x11
	adc	x12, x13, x15
	extr	x13, x12, x11, #51
	mul	x15, x13, x8
	adds	x14, x15, x14
	lsr	x12, x12, #51
	umulh	x13, x13, x8
	umull	x12, w12, w8
	adc	x12, x13, x12
	and	x10, x10, #0x7ffffffffffff
	and	x13, x5, #0x7ffffffffffff
	and	x9, x9, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	extr	x12, x12, x14, #51
	and	x14, x14, #0x7ffffffffffff
	add	x10, x12, x10
	and	x12, x10, #0x7ffffffffffff
	add	x10, x13, x10, lsr #51
	and	x13, x10, #0x7ffffffffffff
	add	x9, x9, x10, lsr #51
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	lsr	x11, x9, #51
	and	x9, x9, #0x7ffffffffffff
	umaddl	x11, w11, w8, x14
	and	x14, x11, #0x7ffffffffffff
	add	x11, x12, x11, lsr #51
	add	x12, x13, x11, lsr #51
	add	x10, x10, x12, lsr #51
	add	x9, x9, x10, lsr #51
	lsr	x13, x9, #51
	umaddl	x8, w13, w8, x14
	and	x11, x11, #0x7ffffffffffff
	stp	x8, x11, [x27]
	and	x8, x12, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	stp	x8, x10, [x27, #16]
	and	x8, x9, #0x7ffffffffffff
	str	x8, [x27, #32]
	ldp	x29, x30, [sp, #144]
	ldp	x20, x19, [sp, #128]
	ldp	x22, x21, [sp, #112]
	ldp	x24, x23, [sp, #96]
	ldp	x26, x25, [sp, #80]
	ldp	x28, x27, [sp, #64]
	add	sp, sp, #160
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	.cfi_restore w21
	.cfi_restore w22
	.cfi_restore w23
	.cfi_restore w24
	.cfi_restore w25
	.cfi_restore w26
	.cfi_restore w27
	.cfi_restore w28
	ret
	.cfi_endproc

	.p2align	2
__ZN16ecc_x25519_check3ecc6x255196x2551917h51fdf6a36a3c95f0E:
	.cfi_startproc
	stp	x28, x27, [sp, #-96]!
	.cfi_def_cfa_offset 96
	stp	x26, x25, [sp, #16]
	stp	x24, x23, [sp, #32]
	stp	x22, x21, [sp, #48]
	stp	x20, x19, [sp, #64]
	stp	x29, x30, [sp, #80]
	add	x29, sp, #80
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	sub	sp, sp, #864
	str	x0, [sp]
	mov	x20, #0
	mov	w8, #0
	ldp	q0, q1, [x1]
	sub	x9, x29, #144
	stp	q0, q1, [sp, #144]
	ldrb	w10, [sp, #144]
	and	w10, w10, #0xf8
	strb	w10, [sp, #144]
	ldrb	w10, [sp, #175]
	mov	w11, #64
	bfxil	w11, w10, #0, #6
	strb	w11, [sp, #175]
	ldrb	w10, [x2, #31]
	and	w10, w10, #0x7f
	ldr	q0, [x2]
	stur	q0, [x29, #-144]
	ldur	q0, [x2, #15]
	stur	q0, [x9, #15]
	sturb	w10, [x29, #-113]
	ldur	x10, [x29, #-144]
	ldur	x11, [x29, #-120]
	and	x12, x10, #0x7ffffffffffff
	ldur	x10, [x9, #6]
	ubfx	x24, x10, #3, #51
	ldur	x10, [x9, #12]
	ubfx	x28, x10, #6, #51
	ldur	x9, [x9, #19]
	ubfx	x25, x9, #1, #51
	ubfx	x21, x11, #12, #51
	stp	x12, x24, [sp, #184]
Lloh0:
	adrp	x9, l_anon.1519422ea68d41d5036578e62b7be619.0@PAGE
Lloh1:
	add	x9, x9, l_anon.1519422ea68d41d5036578e62b7be619.0@PAGEOFF
	ldp	q0, q1, [x9]
	stp	q0, q1, [sp, #224]
	movi.2d	v0, #0000000000000000
	mov	w23, #1
	stp	x28, x25, [sp, #200]
	mov	w9, #254
	str	x21, [sp, #216]
	str	xzr, [sp, #256]
	str	xzr, [sp, #304]
	mov	x27, #0
	mov	x22, #0
	mov	x19, #0
	str	xzr, [sp, #8]
	stp	q0, q0, [sp, #272]
LBB1_1:
	str	x9, [sp, #136]
	lsr	x9, x9, #3
	add	x10, sp, #144
	ldrb	w9, [x10, x9]
	ldr	x10, [sp, #136]
	and	x10, x10, #0x7
	lsr	x9, x9, x10
	and	x9, x9, #0x1
	str	x9, [sp, #128]
	eor	w8, w9, w8
	ldr	x16, [sp, #8]
	sub	x14, x16, w8, uxtb
	ldp	x8, x9, [sp, #224]
	eor	x10, x8, x12
	and	x10, x10, x14
	eor	x0, x10, x8
	eor	x8, x10, x12
	str	x8, [sp, #120]
	eor	x8, x9, x24
	and	x8, x8, x14
	eor	x1, x8, x9
	eor	x8, x8, x24
	str	x8, [sp, #104]
	ldp	x8, x9, [sp, #240]
	eor	x10, x8, x28
	and	x10, x10, x14
	eor	x4, x10, x8
	eor	x8, x10, x28
	str	x8, [sp, #88]
	eor	x8, x9, x25
	and	x8, x8, x14
	eor	x3, x8, x9
	eor	x8, x8, x25
	str	x8, [sp, #96]
	ldr	x8, [sp, #256]
	eor	x9, x8, x21
	and	x9, x9, x14
	eor	x2, x9, x8
	eor	x8, x9, x21
	str	x8, [sp, #112]
	ldp	x8, x9, [sp, #272]
	eor	x10, x8, x23
	and	x10, x10, x14
	eor	x5, x10, x8
	eor	x21, x10, x23
	eor	x8, x9, x27
	and	x8, x8, x14
	eor	x6, x8, x9
	eor	x8, x8, x27
	str	x8, [sp, #72]
	ldp	x8, x9, [sp, #288]
	eor	x10, x8, x22
	and	x10, x10, x14
	eor	x24, x10, x8
	eor	x8, x10, x22
	str	x8, [sp, #56]
	eor	x8, x9, x19
	and	x8, x8, x14
	eor	x15, x8, x9
	eor	x8, x8, x19
	str	x8, [sp, #64]
	ldr	x8, [sp, #304]
	eor	x9, x8, x20
	and	x9, x9, x14
	eor	x7, x9, x8
	eor	x8, x9, x20
	str	x8, [sp, #80]
	add	x8, x5, x0
	add	x9, x6, x1
	and	x10, x8, #0x7ffffffffffff
	adds	x8, x9, x8, lsr #51
	cset	w9, hs
	and	x19, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x24, x4
	adds	x8, x8, x9
	cset	w9, hs
	and	x20, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x15, x3
	adds	x8, x8, x9
	cset	w9, hs
	and	x22, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x7, x2
	adds	x8, x8, x9
	cset	w9, hs
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	mov	w14, #19
	umaddl	x9, w9, w14, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x19, x9, lsr #51
	and	x19, x9, #0x7ffffffffffff
	add	x9, x20, x9, lsr #51
	and	x20, x9, #0x7ffffffffffff
	add	x9, x22, x9, lsr #51
	and	x22, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w14, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x19, x9, lsr #51
	add	x25, x20, x9, lsr #51
	add	x22, x22, x25, lsr #51
	add	x8, x8, x22, lsr #51
	lsr	x19, x8, #51
	umaddl	x20, w19, w14, x10
	and	x19, x9, #0x7ffffffffffff
	and	x25, x25, #0x7ffffffffffff
	and	x26, x22, #0x7ffffffffffff
	and	x8, x8, #0x7ffffffffffff
	umulh	x9, x20, x20
	mul	x10, x20, x20
	umulh	x22, x8, x19
	mul	x27, x8, x19
	mul	x28, x26, x25
	adds	x27, x27, x28
	umulh	x28, x26, x25
	adc	x22, x22, x28
	mov	w17, #38
	umulh	x28, x27, x17
	madd	x28, x22, x17, x28
	mul	x22, x27, x17
	adds	x22, x22, x10
	adc	x27, x28, x9
	umulh	x9, x20, x26
	mul	x10, x20, x26
	mul	x28, x25, x19
	adds	x10, x10, x28
	umulh	x28, x25, x19
	adc	x9, x9, x28
	mul	x28, x8, x8
	umulh	x30, x8, x8
	umulh	x11, x28, x14
	madd	x11, x30, x14, x11
	extr	x9, x9, x10, #63
	mul	x28, x28, x14
	adds	x28, x28, x10, lsl #1
	adc	x30, x9, x11
	mul	x9, x20, x8
	mul	x10, x26, x19
	adds	x9, x9, x10
	umulh	x10, x20, x8
	umulh	x11, x26, x19
	adc	x10, x10, x11
	extr	x11, x10, x9, #63
	mul	x10, x26, x26
	umulh	x12, x26, x26
	umulh	x13, x10, x14
	madd	x13, x12, x14, x13
	mul	x12, x25, x25
	adds	x12, x12, x9, lsl #1
	umulh	x9, x25, x25
	adc	x9, x11, x9
	stp	x26, x8, [sp, #424]
	mul	x8, x8, x17
	mul	x10, x10, x14
	mul	x11, x8, x25
	adds	x10, x11, x10
	umulh	x11, x8, x25
	adc	x11, x11, x13
	lsl	x13, x20, #1
	mul	x14, x13, x19
	adds	x10, x10, x14
	umulh	x14, x8, x26
	mul	x8, x8, x26
	stp	x19, x25, [sp, #408]
	umulh	x26, x13, x25
	mul	x25, x13, x25
	umulh	x13, x13, x19
	extr	x23, x27, x22, #51
	adc	x11, x11, x13
	adds	x10, x10, x23
	lsr	x13, x27, #51
	adc	x11, x11, x13
	mul	x13, x19, x19
	adds	x8, x8, x13
	str	x20, [sp, #400]
	umulh	x13, x19, x19
	adc	x13, x14, x13
	adds	x8, x8, x25
	stp	x0, x1, [sp, #224]
	stp	x4, x3, [sp, #240]
	stp	x5, x6, [sp, #272]
	stp	x24, x15, [sp, #288]
	and	x14, x22, #0x7ffffffffffff
	adc	x13, x13, x26
	and	x19, x10, #0x7ffffffffffff
	extr	x10, x11, x10, #51
	lsr	x11, x11, #51
	adds	x8, x8, x10
	adc	x10, x13, x11
	and	x11, x8, #0x7ffffffffffff
	extr	x8, x10, x8, #51
	lsr	x10, x10, #51
	adds	x8, x28, x8
	adc	x10, x30, x10
	and	x13, x8, #0x7ffffffffffff
	extr	x8, x10, x8, #51
	lsr	x10, x10, #51
	adds	x8, x12, x8
	adc	x9, x9, x10
	lsr	x10, x9, #51
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	mov	w20, #19
	umulh	x12, x9, x20
	umull	x10, w10, w20
	mul	x9, x9, x20
	mov	w26, #19
	adds	x9, x9, x14
	adc	x10, x12, x10
	and	x12, x9, #0x7ffffffffffff
	extr	x9, x10, x9, #51
	add	x9, x9, x19
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x13, x9, lsr #51
	and	x13, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x12
	and	x12, x9, #0x7ffffffffffff
	add	x19, x10, x9, lsr #51
	add	x20, x11, x19, lsr #51
	add	x22, x13, x20, lsr #51
	add	x25, x8, x22, lsr #51
	lsr	x8, x25, #51
	subs	x9, x0, x5
	sbc	x10, x16, x16
	mov	x11, #-38
	movk	x11, #15, lsl #48
	adds	x9, x9, x11
	cinc	x10, x10, hs
	subs	x11, x1, x6
	sbc	x13, x16, x16
	extr	x14, x10, x9, #51
	adds	x11, x11, x14
	lsr	x10, x10, #51
	adc	x10, x13, x10
	subs	x13, x4, x24
	sbc	x14, x16, x16
	mov	x27, #4503599627370494
	adds	x11, x11, x27
	cinc	x10, x10, hs
	extr	x0, x10, x11, #51
	adds	x13, x13, x0
	lsr	x10, x10, #51
	adc	x10, x14, x10
	subs	x14, x3, x15
	sbc	x15, x16, x16
	adds	x13, x13, x27
	cinc	x10, x10, hs
	extr	x0, x10, x13, #51
	adds	x14, x14, x0
	lsr	x10, x10, #51
	adc	x10, x15, x10
	subs	x15, x2, x7
	sbc	x0, x16, x16
	mov	x28, x16
	adds	x14, x14, x27
	cinc	x10, x10, hs
	extr	x1, x10, x14, #51
	umaddl	x16, w8, w26, x12
	adds	x8, x15, x1
	lsr	x10, x10, #51
	adc	x10, x0, x10
	and	x9, x9, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	and	x12, x13, #0x7ffffffffffff
	and	x13, x14, #0x7ffffffffffff
	adds	x8, x8, x27
	cinc	x10, x10, hs
	extr	x14, x10, x8, #51
	lsr	x10, x10, #51
	and	x8, x8, #0x7ffffffffffff
	umulh	x15, x14, x26
	umull	x10, w10, w26
	mul	x14, x14, x26
	adds	x9, x14, x9
	adc	x10, x15, x10
	and	x14, x9, #0x7ffffffffffff
	extr	x9, x10, x9, #51
	add	x9, x9, x11
	and	x10, x9, #0x7ffffffffffff
	add	x9, x12, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x13, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x14
	and	x13, x9, #0x7ffffffffffff
	add	x9, x10, x9, lsr #51
	add	x10, x11, x9, lsr #51
	add	x11, x12, x10, lsr #51
	add	x8, x8, x11, lsr #51
	lsr	x12, x8, #51
	umaddl	x1, w12, w26, x13
	and	x0, x9, #0x7ffffffffffff
	and	x3, x10, #0x7ffffffffffff
	and	x4, x11, #0x7ffffffffffff
	and	x11, x8, #0x7ffffffffffff
	stp	x3, x4, [sp, #496]
	str	x11, [sp, #512]
	umulh	x8, x1, x1
	mul	x9, x1, x1
	umulh	x10, x11, x0
	mul	x12, x11, x0
	umulh	x13, x4, x3
	mul	x14, x4, x3
	adds	x12, x12, x14
	adc	x10, x10, x13
	umulh	x13, x12, x17
	madd	x10, x10, x17, x13
	mul	x12, x12, x17
	adds	x2, x12, x9
	adc	x8, x10, x8
	lsl	x9, x1, #1
	mul	x12, x4, x4
	umulh	x10, x4, x4
	umulh	x13, x12, x26
	madd	x13, x10, x26, x13
	umulh	x10, x1, x4
	mul	x14, x1, x4
	umulh	x15, x3, x0
	mul	x5, x3, x0
	adds	x14, x14, x5
	adc	x10, x10, x15
	mul	x15, x11, x11
	umulh	x5, x11, x11
	umulh	x6, x15, x26
	madd	x5, x5, x26, x6
	mul	x6, x11, x17
	extr	x7, x10, x14, #63
	mul	x10, x15, x26
	adds	x10, x10, x14, lsl #1
	adc	x5, x7, x5
	umulh	x14, x1, x11
	mul	x11, x1, x11
	mul	x15, x4, x0
	adds	x11, x11, x15
	umulh	x15, x4, x0
	adc	x14, x14, x15
	extr	x14, x14, x11, #63
	mul	x15, x3, x3
	adds	x11, x15, x11, lsl #1
	umulh	x15, x3, x3
	adc	x14, x14, x15
	mul	x15, x6, x3
	mul	x12, x12, x26
	adds	x12, x15, x12
	umulh	x15, x6, x4
	mul	x4, x6, x4
	umulh	x6, x6, x3
	adc	x13, x6, x13
	mul	x6, x9, x0
	adds	x12, x12, x6
	umulh	x6, x9, x3
	mul	x3, x9, x3
	umulh	x9, x9, x0
	adc	x9, x13, x9
	extr	x13, x8, x2, #51
	adds	x12, x12, x13
	lsr	x8, x8, #51
	adc	x8, x9, x8
	mul	x9, x0, x0
	adds	x9, x4, x9
	stp	x1, x0, [sp, #480]
	umulh	x13, x0, x0
	adc	x13, x15, x13
	and	x17, x19, #0x7ffffffffffff
	and	x1, x20, #0x7ffffffffffff
	and	x4, x22, #0x7ffffffffffff
	and	x7, x25, #0x7ffffffffffff
	and	x15, x2, #0x7ffffffffffff
	adds	x9, x9, x3
	and	x0, x12, #0x7ffffffffffff
	extr	x12, x8, x12, #51
	lsr	x8, x8, #51
	adc	x13, x13, x6
	adds	x9, x9, x12
	adc	x8, x13, x8
	and	x12, x9, #0x7ffffffffffff
	extr	x9, x8, x9, #51
	lsr	x8, x8, #51
	adds	x9, x10, x9
	adc	x8, x5, x8
	and	x10, x9, #0x7ffffffffffff
	extr	x9, x8, x9, #51
	lsr	x8, x8, #51
	adds	x9, x11, x9
	adc	x8, x14, x8
	lsr	x11, x8, #51
	extr	x8, x8, x9, #51
	and	x9, x9, #0x7ffffffffffff
	umulh	x13, x8, x26
	umull	x11, w11, w26
	mul	x8, x8, x26
	adds	x8, x8, x15
	adc	x11, x13, x11
	and	x13, x8, #0x7ffffffffffff
	extr	x8, x11, x8, #51
	add	x8, x8, x0
	and	x11, x8, #0x7ffffffffffff
	add	x8, x12, x8, lsr #51
	and	x12, x8, #0x7ffffffffffff
	add	x8, x10, x8, lsr #51
	and	x10, x8, #0x7ffffffffffff
	add	x8, x9, x8, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x13
	and	x13, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	add	x11, x12, x9, lsr #51
	add	x10, x10, x11, lsr #51
	add	x8, x8, x10, lsr #51
	lsr	x12, x8, #51
	umaddl	x12, w12, w26, x13
	and	x9, x9, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	and	x8, x8, #0x7ffffffffffff
	str	x12, [sp, #520]
	str	x9, [sp, #528]
	str	x11, [sp, #536]
	str	x10, [sp, #544]
	stp	x17, x16, [sp, #40]
	subs	x12, x16, x12
	sbc	x13, x28, x28
	mov	x14, #-38
	movk	x14, #15, lsl #48
	adds	x12, x12, x14
	cinc	x13, x13, hs
	subs	x9, x17, x9
	sbc	x14, x28, x28
	extr	x15, x13, x12, #51
	adds	x9, x9, x15
	lsr	x13, x13, #51
	adc	x13, x14, x13
	subs	x11, x1, x11
	mov	x5, x1
	stp	x4, x1, [sp, #24]
	sbc	x14, x28, x28
	adds	x9, x9, x27
	cinc	x13, x13, hs
	extr	x15, x13, x9, #51
	adds	x11, x11, x15
	lsr	x13, x13, #51
	adc	x13, x14, x13
	subs	x10, x4, x10
	mov	x6, x4
	sbc	x14, x28, x28
	adds	x11, x11, x27
	cinc	x13, x13, hs
	extr	x15, x13, x11, #51
	str	x8, [sp, #552]
	adds	x10, x10, x15
	lsr	x13, x13, #51
	adc	x13, x14, x13
	subs	x8, x7, x8
	str	x7, [sp, #16]
	sbc	x14, x28, x28
	adds	x10, x10, x27
	cinc	x13, x13, hs
	extr	x15, x13, x10, #51
	adds	x8, x8, x15
	lsr	x13, x13, #51
	adc	x13, x14, x13
	and	x12, x12, #0x7ffffffffffff
	and	x9, x9, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	adds	x8, x8, x27
	cinc	x13, x13, hs
	lsr	x14, x13, #51
	extr	x13, x13, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umulh	x15, x13, x26
	umull	x14, w14, w26
	mul	x13, x13, x26
	adds	x12, x13, x12
	adc	x13, x15, x14
	extr	x13, x13, x12, #51
	and	x12, x12, #0x7ffffffffffff
	add	x0, x13, x9
	add	x1, x11, x0, lsr #51
	add	x2, x10, x1, lsr #51
	add	x3, x8, x2, lsr #51
	lsr	x8, x3, #51
	umaddl	x4, w8, w26, x12
	ldr	x14, [sp, #120]
	mov	x23, x21
	add	x8, x21, x14
	ldp	x20, x19, [sp, #96]
	ldp	x30, x25, [sp, #64]
	add	x9, x25, x19
	and	x10, x8, #0x7ffffffffffff
	adds	x8, x9, x8, lsr #51
	cset	w9, hs
	and	x11, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	ldp	x24, x22, [sp, #80]
	ldr	x21, [sp, #56]
	add	x9, x21, x22
	adds	x8, x8, x9
	cset	w9, hs
	and	x12, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x30, x20
	adds	x8, x8, x9
	cset	w9, hs
	and	x13, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	ldr	x15, [sp, #112]
	add	x9, x24, x15
	adds	x8, x8, x9
	cset	w9, hs
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x12, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x13, x9, lsr #51
	and	x13, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	add	x11, x12, x9, lsr #51
	add	x12, x13, x11, lsr #51
	add	x8, x8, x12, lsr #51
	lsr	x13, x8, #51
	umaddl	x10, w13, w26, x10
	and	x9, x9, #0x7ffffffffffff
	str	x10, [sp, #600]
	str	x9, [sp, #608]
	and	x9, x11, #0x7ffffffffffff
	and	x10, x12, #0x7ffffffffffff
	str	x9, [sp, #616]
	str	x10, [sp, #624]
	and	x8, x8, #0x7ffffffffffff
	str	x8, [sp, #632]
	subs	x8, x14, x23
	sbc	x9, x28, x28
	mov	x10, #-38
	movk	x10, #15, lsl #48
	adds	x8, x8, x10
	cinc	x9, x9, hs
	subs	x10, x19, x25
	sbc	x11, x28, x28
	extr	x12, x9, x8, #51
	adds	x10, x10, x12
	lsr	x9, x9, #51
	adc	x9, x11, x9
	subs	x11, x22, x21
	sbc	x12, x28, x28
	adds	x10, x10, x27
	cinc	x9, x9, hs
	extr	x13, x9, x10, #51
	adds	x11, x11, x13
	lsr	x9, x9, #51
	adc	x9, x12, x9
	subs	x12, x20, x30
	sbc	x13, x28, x28
	adds	x11, x11, x27
	cinc	x9, x9, hs
	extr	x14, x9, x11, #51
	adds	x12, x12, x14
	lsr	x9, x9, #51
	adc	x9, x13, x9
	subs	x13, x15, x24
	sbc	x14, x28, x28
	adds	x12, x12, x27
	cinc	x9, x9, hs
	extr	x15, x9, x12, #51
	adds	x13, x13, x15
	lsr	x9, x9, #51
	adc	x9, x14, x9
	and	x8, x8, #0x7ffffffffffff
	adds	x13, x13, x27
	cinc	x9, x9, hs
	extr	x14, x9, x13, #51
	mul	x15, x14, x26
	adds	x8, x15, x8
	lsr	x9, x9, #51
	umulh	x14, x14, x26
	umull	x9, w9, w26
	adc	x9, x14, x9
	and	x10, x10, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	and	x12, x12, #0x7ffffffffffff
	and	x13, x13, #0x7ffffffffffff
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	add	x9, x9, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x12, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x13, x9, lsr #51
	lsr	x13, x9, #51
	and	x9, x9, #0x7ffffffffffff
	umaddl	x8, w13, w26, x8
	and	x13, x8, #0x7ffffffffffff
	add	x8, x10, x8, lsr #51
	add	x10, x11, x8, lsr #51
	add	x11, x12, x10, lsr #51
	add	x9, x9, x11, lsr #51
	lsr	x12, x9, #51
	umaddl	x12, w12, w26, x13
	and	x8, x8, #0x7ffffffffffff
	str	x12, [sp, #640]
	str	x8, [sp, #648]
	and	x8, x10, #0x7ffffffffffff
	and	x10, x11, #0x7ffffffffffff
	str	x8, [sp, #656]
	str	x10, [sp, #664]
	and	x8, x0, #0x7ffffffffffff
	and	x10, x1, #0x7ffffffffffff
	and	x11, x2, #0x7ffffffffffff
	and	x12, x3, #0x7ffffffffffff
	add	x8, x8, x4, lsr #51
	and	x13, x4, #0x7ffffffffffff
	add	x10, x10, x8, lsr #51
	add	x11, x11, x10, lsr #51
	add	x12, x12, x11, lsr #51
	and	x9, x9, #0x7ffffffffffff
	str	x9, [sp, #672]
	lsr	x9, x12, #51
	umaddl	x9, w9, w26, x13
	stp	x16, x17, [sp, #440]
	stp	x5, x6, [sp, #456]
	str	x7, [sp, #472]
	and	x13, x8, #0x7ffffffffffff
	stp	x9, x13, [sp, #88]
	and	x10, x10, #0x7ffffffffffff
	and	x8, x11, #0x7ffffffffffff
	stp	x10, x8, [sp, #104]
	str	x9, [sp, #560]
	str	x13, [sp, #568]
	and	x9, x12, #0x7ffffffffffff
	str	x9, [sp, #120]
	str	x10, [sp, #576]
	str	x8, [sp, #584]
	str	x9, [sp, #592]
	add	x0, sp, #680
	add	x1, sp, #640
	add	x2, sp, #400
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	sub	x0, x29, #224
	add	x1, sp, #600
	add	x2, sp, #480
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x0, [sp, #680]
	ldr	x17, [sp, #688]
	ldp	x3, x1, [x29, #-224]
	add	x8, x3, x0
	add	x9, x1, x17
	ldr	x16, [sp, #696]
	ldr	x14, [sp, #704]
	ldp	x2, x15, [x29, #-208]
	ldr	x12, [sp, #712]
	ldur	x13, [x29, #-192]
	and	x10, x8, #0x7ffffffffffff
	adds	x8, x9, x8, lsr #51
	cset	w9, hs
	and	x11, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x2, x16
	adds	x8, x8, x9
	cset	w9, hs
	and	x4, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x15, x14
	adds	x8, x8, x9
	cset	w9, hs
	and	x5, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x9, x13, x12
	adds	x8, x8, x9
	cset	w9, hs
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x4, x9, lsr #51
	and	x4, x9, #0x7ffffffffffff
	add	x9, x5, x9, lsr #51
	and	x5, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x10
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	add	x11, x4, x9, lsr #51
	add	x4, x5, x11, lsr #51
	add	x5, x8, x4, lsr #51
	lsr	x8, x5, #51
	umaddl	x6, w8, w26, x10
	and	x8, x9, #0x7ffffffffffff
	and	x10, x11, #0x7ffffffffffff
	and	x11, x4, #0x7ffffffffffff
	and	x5, x5, #0x7ffffffffffff
	umulh	x4, x6, x6
	mul	x9, x6, x6
	umulh	x7, x5, x8
	mul	x19, x5, x8
	umulh	x20, x11, x10
	mul	x22, x11, x10
	adds	x19, x19, x22
	adc	x7, x7, x20
	mov	w22, #38
	umulh	x20, x19, x22
	madd	x7, x7, x22, x20
	mul	x19, x19, x22
	mov	w30, #38
	adds	x9, x19, x9
	adc	x4, x7, x4
	mul	x7, x11, x11
	umulh	x19, x11, x11
	umulh	x20, x7, x26
	madd	x19, x19, x26, x20
	umulh	x20, x6, x11
	mul	x22, x6, x11
	umulh	x23, x10, x8
	mul	x24, x10, x8
	adds	x22, x22, x24
	adc	x20, x20, x23
	mul	x23, x5, x5
	umulh	x24, x5, x5
	umulh	x25, x23, x26
	madd	x24, x24, x26, x25
	lsl	x25, x6, #1
	extr	x20, x20, x22, #63
	mul	x23, x23, x26
	adds	x22, x23, x22, lsl #1
	mul	x23, x5, x30
	adc	x20, x20, x24
	umulh	x24, x6, x5
	mul	x5, x6, x5
	mul	x6, x11, x8
	adds	x5, x5, x6
	umulh	x6, x11, x8
	adc	x6, x24, x6
	extr	x6, x6, x5, #63
	mul	x24, x10, x10
	adds	x5, x24, x5, lsl #1
	umulh	x24, x10, x10
	adc	x6, x6, x24
	mul	x24, x23, x10
	mul	x7, x7, x26
	adds	x7, x24, x7
	umulh	x24, x23, x11
	mul	x11, x23, x11
	umulh	x23, x23, x10
	adc	x19, x23, x19
	mul	x23, x25, x8
	adds	x7, x7, x23
	umulh	x23, x25, x10
	mul	x10, x25, x10
	umulh	x25, x25, x8
	adc	x19, x19, x25
	extr	x25, x4, x9, #51
	adds	x7, x7, x25
	lsr	x4, x4, #51
	adc	x4, x19, x4
	mul	x19, x8, x8
	adds	x11, x11, x19
	umulh	x8, x8, x8
	adc	x8, x24, x8
	adds	x10, x11, x10
	adc	x8, x8, x23
	extr	x11, x4, x7, #51
	adds	x10, x10, x11
	and	x9, x9, #0x7ffffffffffff
	and	x11, x7, #0x7ffffffffffff
	lsr	x4, x4, #51
	adc	x8, x8, x4
	and	x4, x10, #0x7ffffffffffff
	extr	x10, x8, x10, #51
	lsr	x8, x8, #51
	adds	x10, x22, x10
	adc	x8, x20, x8
	and	x7, x10, #0x7ffffffffffff
	extr	x10, x8, x10, #51
	lsr	x8, x8, #51
	adds	x10, x5, x10
	adc	x8, x6, x8
	lsr	x5, x8, #51
	extr	x8, x8, x10, #51
	and	x10, x10, #0x7ffffffffffff
	umulh	x6, x8, x26
	umull	x5, w5, w26
	mul	x8, x8, x26
	adds	x8, x8, x9
	adc	x9, x6, x5
	and	x5, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	add	x8, x8, x11
	and	x9, x8, #0x7ffffffffffff
	add	x8, x4, x8, lsr #51
	and	x11, x8, #0x7ffffffffffff
	add	x8, x7, x8, lsr #51
	and	x4, x8, #0x7ffffffffffff
	add	x8, x10, x8, lsr #51
	lsr	x10, x8, #51
	and	x6, x8, #0x7ffffffffffff
	umaddl	x8, w10, w26, x5
	and	x5, x8, #0x7ffffffffffff
	add	x8, x9, x8, lsr #51
	add	x9, x11, x8, lsr #51
	add	x10, x4, x9, lsr #51
	add	x11, x6, x10, lsr #51
	lsr	x4, x11, #51
	subs	x0, x0, x3
	sbc	x3, x28, x28
	mov	x6, #-38
	movk	x6, #15, lsl #48
	adds	x0, x0, x6
	cinc	x3, x3, hs
	subs	x17, x17, x1
	sbc	x1, x28, x28
	extr	x6, x3, x0, #51
	adds	x17, x17, x6
	lsr	x3, x3, #51
	adc	x1, x1, x3
	subs	x16, x16, x2
	sbc	x2, x28, x28
	adds	x17, x17, x27
	cinc	x1, x1, hs
	extr	x3, x1, x17, #51
	adds	x16, x16, x3
	lsr	x1, x1, #51
	adc	x1, x2, x1
	subs	x14, x14, x15
	sbc	x15, x28, x28
	adds	x16, x16, x27
	cinc	x1, x1, hs
	extr	x2, x1, x16, #51
	adds	x14, x14, x2
	lsr	x1, x1, #51
	adc	x15, x15, x1
	subs	x12, x12, x13
	sbc	x13, x28, x28
	adds	x14, x14, x27
	cinc	x15, x15, hs
	extr	x1, x15, x14, #51
	umaddl	x2, w4, w26, x5
	str	x2, [sp, #80]
	adds	x12, x12, x1
	lsr	x15, x15, #51
	adc	x13, x13, x15
	and	x15, x0, #0x7ffffffffffff
	and	x17, x17, #0x7ffffffffffff
	and	x16, x16, #0x7ffffffffffff
	and	x14, x14, #0x7ffffffffffff
	adds	x12, x12, x27
	cinc	x13, x13, hs
	extr	x0, x13, x12, #51
	lsr	x13, x13, #51
	and	x12, x12, #0x7ffffffffffff
	umulh	x1, x0, x26
	umull	x13, w13, w26
	mul	x0, x0, x26
	adds	x15, x0, x15
	adc	x13, x1, x13
	and	x0, x15, #0x7ffffffffffff
	extr	x13, x13, x15, #51
	add	x13, x13, x17
	and	x15, x13, #0x7ffffffffffff
	add	x13, x16, x13, lsr #51
	and	x16, x13, #0x7ffffffffffff
	add	x13, x14, x13, lsr #51
	and	x14, x13, #0x7ffffffffffff
	add	x12, x12, x13, lsr #51
	lsr	x13, x12, #51
	and	x12, x12, #0x7ffffffffffff
	umaddl	x13, w13, w26, x0
	and	x17, x13, #0x7ffffffffffff
	add	x13, x15, x13, lsr #51
	add	x15, x16, x13, lsr #51
	add	x16, x14, x15, lsr #51
	add	x0, x12, x16, lsr #51
	lsr	x12, x0, #51
	umaddl	x1, w12, w26, x17
	and	x12, x13, #0x7ffffffffffff
	and	x14, x15, #0x7ffffffffffff
	and	x15, x16, #0x7ffffffffffff
	and	x0, x0, #0x7ffffffffffff
	umulh	x16, x1, x1
	mul	x13, x1, x1
	umulh	x17, x0, x12
	mul	x2, x0, x12
	umulh	x3, x15, x14
	mul	x4, x15, x14
	adds	x2, x2, x4
	adc	x17, x17, x3
	umulh	x3, x2, x30
	madd	x17, x17, x30, x3
	mul	x2, x2, x30
	adds	x13, x2, x13
	adc	x16, x17, x16
	lsl	x17, x1, #1
	mul	x2, x0, x30
	umulh	x3, x1, x15
	mul	x4, x1, x15
	umulh	x5, x14, x12
	mul	x6, x14, x12
	adds	x4, x4, x6
	adc	x3, x3, x5
	mul	x5, x0, x0
	umulh	x6, x0, x0
	umulh	x7, x5, x26
	madd	x6, x6, x26, x7
	mul	x7, x15, x15
	extr	x3, x3, x4, #63
	mul	x5, x5, x26
	adds	x4, x5, x4, lsl #1
	umulh	x5, x15, x15
	adc	x3, x3, x6
	umulh	x6, x1, x0
	mul	x0, x1, x0
	mul	x1, x15, x12
	adds	x0, x0, x1
	umulh	x1, x15, x12
	adc	x1, x6, x1
	umulh	x6, x7, x26
	extr	x1, x1, x0, #63
	madd	x5, x5, x26, x6
	mul	x6, x14, x14
	adds	x0, x6, x0, lsl #1
	umulh	x6, x14, x14
	adc	x1, x1, x6
	mul	x6, x2, x14
	mul	x7, x7, x26
	adds	x6, x6, x7
	umulh	x7, x2, x15
	mul	x15, x2, x15
	umulh	x2, x2, x14
	adc	x2, x2, x5
	mul	x5, x17, x12
	adds	x5, x6, x5
	umulh	x6, x17, x14
	mul	x14, x17, x14
	umulh	x17, x17, x12
	adc	x17, x2, x17
	extr	x2, x16, x13, #51
	adds	x2, x5, x2
	lsr	x16, x16, #51
	adc	x16, x17, x16
	mul	x17, x12, x12
	adds	x15, x15, x17
	umulh	x12, x12, x12
	adc	x12, x7, x12
	adds	x14, x15, x14
	extr	x15, x16, x2, #51
	adc	x12, x12, x6
	adds	x14, x14, x15
	lsr	x15, x16, #51
	adc	x12, x12, x15
	extr	x15, x12, x14, #51
	adds	x15, x4, x15
	and	x13, x13, #0x7ffffffffffff
	lsr	x12, x12, #51
	adc	x12, x3, x12
	extr	x16, x12, x15, #51
	lsr	x12, x12, #51
	adds	x16, x0, x16
	adc	x12, x1, x12
	extr	x17, x12, x16, #51
	mul	x0, x17, x26
	adds	x13, x0, x13
	lsr	x12, x12, #51
	umulh	x17, x17, x26
	umull	x12, w12, w26
	adc	x12, x17, x12
	and	x17, x2, #0x7ffffffffffff
	and	x14, x14, #0x7ffffffffffff
	and	x15, x15, #0x7ffffffffffff
	and	x16, x16, #0x7ffffffffffff
	extr	x12, x12, x13, #51
	and	x13, x13, #0x7ffffffffffff
	add	x12, x12, x17
	and	x17, x12, #0x7ffffffffffff
	add	x12, x14, x12, lsr #51
	and	x14, x12, #0x7ffffffffffff
	add	x12, x15, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x16, x12, lsr #51
	lsr	x16, x12, #51
	and	x12, x12, #0x7ffffffffffff
	umaddl	x13, w16, w26, x13
	and	x16, x13, #0x7ffffffffffff
	add	x13, x17, x13, lsr #51
	add	x14, x14, x13, lsr #51
	add	x15, x15, x14, lsr #51
	add	x12, x12, x15, lsr #51
	lsr	x17, x12, #51
	umaddl	x16, w17, w26, x16
	and	x13, x13, #0x7ffffffffffff
	stp	x16, x13, [x29, #-144]
	and	x13, x14, #0x7ffffffffffff
	and	x14, x15, #0x7ffffffffffff
	stp	x13, x14, [x29, #-128]
	and	x12, x12, #0x7ffffffffffff
	stur	x12, [x29, #-112]
	and	x24, x8, #0x7ffffffffffff
	and	x28, x9, #0x7ffffffffffff
	and	x25, x10, #0x7ffffffffffff
	and	x21, x11, #0x7ffffffffffff
	sub	x0, x29, #184
	add	x1, sp, #184
	sub	x2, x29, #144
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldp	x23, x27, [x29, #-184]
	ldp	x22, x19, [x29, #-168]
	ldur	x20, [x29, #-152]
	add	x0, sp, #224
	add	x1, sp, #440
	add	x2, sp, #520
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	mov	w17, #56129
	movk	w17, #1, lsl #16
	ldp	x9, x11, [sp, #88]
	umulh	x8, x9, x17
	mul	x9, x9, x17
	umulh	x10, x11, x17
	mul	x11, x11, x17
	ldp	x13, x15, [sp, #104]
	umulh	x12, x13, x17
	mul	x13, x13, x17
	umulh	x14, x15, x17
	mul	x15, x15, x17
	ldr	x0, [sp, #120]
	umulh	x16, x0, x17
	mul	x17, x0, x17
	and	x0, x9, #0x7ffffffffffff
	extr	x8, x8, x9, #51
	adds	x8, x8, x11
	cinc	x9, x10, hs
	and	x10, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	adds	x8, x8, x13
	cinc	x9, x12, hs
	and	x11, x8, #0x7ffffffffffff
	extr	x8, x9, x8, #51
	adds	x8, x8, x15
	cinc	x9, x14, hs
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	adds	x9, x9, x17
	cinc	x12, x16, hs
	extr	x12, x12, x9, #51
	and	x9, x9, #0x7ffffffffffff
	madd	x12, x12, x26, x0
	and	x13, x12, #0x7ffffffffffff
	add	x10, x10, x12, lsr #51
	and	x12, x10, #0x7ffffffffffff
	add	x10, x11, x10, lsr #51
	and	x11, x10, #0x7ffffffffffff
	add	x8, x8, x10, lsr #51
	and	x10, x8, #0x7ffffffffffff
	add	x8, x9, x8, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x13
	and	x13, x9, #0x7ffffffffffff
	add	x9, x12, x9, lsr #51
	add	x11, x11, x9, lsr #51
	add	x10, x10, x11, lsr #51
	add	x8, x8, x10, lsr #51
	lsr	x12, x8, #51
	umaddl	x12, w12, w26, x13
	and	x9, x9, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	and	x8, x8, #0x7ffffffffffff
	ldp	x13, x15, [sp, #40]
	add	x12, x12, x15
	add	x9, x9, x13
	ldp	x13, x15, [sp, #24]
	add	x11, x11, x15
	add	x10, x10, x13
	ldr	x13, [sp, #16]
	add	x8, x8, x13
	and	x13, x12, #0x7ffffffffffff
	add	x9, x9, x12, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x10, x9, lsr #51
	and	x10, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x13
	and	x13, x9, #0x7ffffffffffff
	add	x9, x12, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	and	x11, x9, #0x7ffffffffffff
	add	x9, x10, x9, lsr #51
	and	x10, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w26, x13
	and	x13, x9, #0x7ffffffffffff
	add	x9, x12, x9, lsr #51
	add	x11, x11, x9, lsr #51
	add	x10, x10, x11, lsr #51
	add	x8, x8, x10, lsr #51
	lsr	x12, x8, #51
	umaddl	x12, w12, w26, x13
	and	x9, x9, #0x7ffffffffffff
	and	x11, x11, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	stp	x12, x9, [x29, #-144]
	and	x8, x8, #0x7ffffffffffff
	stp	x11, x10, [x29, #-128]
	stur	x8, [x29, #-112]
	add	x0, sp, #272
	add	x1, sp, #560
	sub	x2, x29, #144
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x12, [sp, #80]
	ldp	x10, x9, [sp, #128]
	sub	x9, x9, #1
	mov	x8, x10
	cmn	x9, #1
	b.ne	LBB1_1
	fmov	d0, x12
	mov.d	v0[1], x24
	ldp	q1, q2, [sp, #224]
	fmov	d3, x28
	mov.d	v3[1], x25
	neg	x12, x10
	and	x8, x21, x12
	dup.2d	v4, x12
	bit.16b	v2, v3, v4
	bif.16b	v0, v1, v4
	stp	q0, q2, [sp, #224]
	ldr	x9, [sp, #256]
	sub	x13, x10, #1
	and	x9, x9, x13
	orr	x8, x8, x9
	str	x8, [sp, #256]
	ldp	x8, x9, [sp, #272]
	and	x8, x8, x13
	and	x10, x23, x12
	orr	x8, x10, x8
	and	x9, x9, x13
	and	x10, x27, x12
	orr	x9, x10, x9
	stp	x8, x9, [sp, #272]
	ldp	x10, x11, [sp, #288]
	and	x10, x10, x13
	and	x14, x22, x12
	orr	x10, x14, x10
	and	x11, x11, x13
	and	x14, x19, x12
	orr	x11, x14, x11
	stp	x10, x11, [sp, #288]
	ldr	x14, [sp, #304]
	and	x13, x14, x13
	and	x12, x20, x12
	orr	x12, x12, x13
	str	x12, [sp, #304]
	umulh	x13, x8, x8
	mul	x14, x8, x8
	umulh	x15, x12, x9
	mul	x16, x12, x9
	umulh	x17, x11, x10
	mul	x0, x11, x10
	adds	x16, x16, x0
	adc	x15, x15, x17
	mov	w19, #38
	umulh	x17, x16, x19
	madd	x15, x15, x19, x17
	mul	x16, x16, x19
	adds	x14, x16, x14
	adc	x13, x15, x13
	lsl	x15, x8, #1
	lsr	x16, x8, #63
	umulh	x17, x15, x9
	madd	x17, x16, x9, x17
	mul	x0, x15, x9
	mul	x1, x12, x19
	umulh	x2, x12, x19
	umulh	x3, x1, x10
	madd	x3, x2, x10, x3
	mul	x4, x1, x10
	mul	x5, x11, x11
	umulh	x6, x11, x11
	mov	w20, #19
	umulh	x7, x5, x20
	madd	x6, x6, x20, x7
	mul	x5, x5, x20
	umulh	x7, x15, x10
	madd	x16, x16, x10, x7
	mul	x15, x15, x10
	umulh	x7, x9, x9
	mul	x21, x9, x9
	adds	x15, x15, x21
	adc	x16, x16, x7
	umulh	x7, x1, x11
	madd	x2, x2, x11, x7
	mul	x1, x1, x11
	adds	x15, x15, x1
	adc	x16, x16, x2
	umulh	x1, x11, x8
	mul	x2, x11, x8
	umulh	x7, x10, x9
	mul	x21, x10, x9
	adds	x2, x2, x21
	adc	x1, x1, x7
	extr	x1, x1, x2, #63
	mul	x7, x12, x12
	umulh	x21, x12, x12
	umulh	x22, x7, x20
	madd	x21, x21, x20, x22
	mul	x7, x7, x20
	adds	x2, x7, x2, lsl #1
	adc	x1, x21, x1
	umulh	x7, x12, x8
	mul	x8, x12, x8
	umulh	x12, x11, x9
	mul	x9, x11, x9
	adds	x8, x8, x9
	adc	x9, x7, x12
	extr	x9, x9, x8, #63
	umulh	x11, x10, x10
	mul	x10, x10, x10
	adds	x8, x10, x8, lsl #1
	adc	x9, x9, x11
	and	x10, x14, #0x7ffffffffffff
	extr	x11, x13, x14, #51
	lsr	x12, x13, #51
	adds	x13, x5, x0
	adc	x14, x6, x17
	adds	x13, x13, x4
	adc	x14, x14, x3
	adds	x11, x13, x11
	adc	x12, x14, x12
	and	x13, x11, #0x7ffffffffffff
	lsr	x14, x12, #51
	extr	x11, x12, x11, #51
	adds	x11, x15, x11
	adc	x12, x16, x14
	and	x14, x11, #0x7ffffffffffff
	lsr	x15, x12, #51
	extr	x11, x12, x11, #51
	adds	x11, x2, x11
	adc	x12, x1, x15
	and	x15, x11, #0x7ffffffffffff
	lsr	x16, x12, #51
	extr	x11, x12, x11, #51
	adds	x8, x8, x11
	adc	x9, x9, x16
	lsr	x11, x9, #51
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umulh	x12, x9, x20
	umull	x11, w11, w20
	mul	x9, x9, x20
	adds	x9, x9, x10
	adc	x10, x12, x11
	extr	x10, x10, x9, #51
	and	x9, x9, #0x7ffffffffffff
	add	x10, x10, x13
	and	x11, x10, #0x7ffffffffffff
	add	x10, x14, x10, lsr #51
	and	x12, x10, #0x7ffffffffffff
	add	x10, x15, x10, lsr #51
	and	x13, x10, #0x7ffffffffffff
	add	x8, x8, x10, lsr #51
	lsr	x10, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w10, w20, x9
	and	x10, x9, #0x7ffffffffffff
	add	x9, x11, x9, lsr #51
	add	x11, x12, x9, lsr #51
	add	x12, x13, x11, lsr #51
	add	x13, x8, x12, lsr #51
	lsr	x8, x13, #51
	umaddl	x10, w8, w20, x10
	and	x9, x9, #0x7ffffffffffff
	and	x8, x11, #0x7ffffffffffff
	and	x11, x12, #0x7ffffffffffff
	and	x12, x13, #0x7ffffffffffff
	stp	x10, x9, [sp, #360]
	stp	x8, x11, [sp, #376]
	str	x12, [sp, #392]
	umulh	x13, x10, x10
	mul	x14, x10, x10
	umulh	x15, x12, x9
	mul	x16, x12, x9
	umulh	x17, x11, x8
	mul	x0, x11, x8
	adds	x16, x16, x0
	adc	x15, x15, x17
	umulh	x17, x16, x19
	madd	x15, x15, x19, x17
	mul	x16, x16, x19
	adds	x14, x16, x14
	adc	x13, x15, x13
	lsl	x15, x10, #1
	umulh	x16, x15, x9
	mul	x17, x15, x9
	mul	x0, x12, x19
	umulh	x1, x0, x8
	mul	x2, x0, x8
	mul	x3, x11, x11
	umulh	x4, x11, x11
	umulh	x5, x3, x20
	madd	x4, x4, x20, x5
	mul	x3, x3, x20
	umulh	x5, x15, x8
	mul	x15, x15, x8
	umulh	x6, x9, x9
	mul	x7, x9, x9
	umulh	x21, x0, x11
	mul	x0, x0, x11
	umulh	x22, x10, x11
	mul	x23, x10, x11
	umulh	x24, x8, x9
	mul	x25, x8, x9
	adds	x23, x23, x25
	adc	x22, x22, x24
	extr	x22, x22, x23, #63
	mul	x24, x12, x12
	umulh	x25, x12, x12
	umulh	x26, x24, x20
	madd	x25, x25, x20, x26
	mul	x24, x24, x20
	adds	x23, x24, x23, lsl #1
	adc	x22, x22, x25
	umulh	x24, x10, x12
	mul	x10, x10, x12
	umulh	x12, x11, x9
	mul	x9, x11, x9
	adds	x9, x10, x9
	adc	x10, x24, x12
	extr	x10, x10, x9, #63
	umulh	x11, x8, x8
	mul	x8, x8, x8
	adds	x8, x8, x9, lsl #1
	adc	x9, x10, x11
	and	x10, x14, #0x7ffffffffffff
	lsr	x11, x13, #51
	extr	x12, x13, x14, #51
	adds	x13, x2, x3
	adc	x14, x1, x4
	adds	x13, x13, x17
	adc	x14, x14, x16
	adds	x12, x13, x12
	adc	x11, x14, x11
	and	x13, x12, #0x7ffffffffffff
	lsr	x14, x11, #51
	extr	x11, x11, x12, #51
	adds	x12, x0, x7
	adc	x16, x21, x6
	adds	x12, x12, x15
	adc	x15, x16, x5
	adds	x11, x12, x11
	adc	x12, x15, x14
	and	x14, x11, #0x7ffffffffffff
	lsr	x15, x12, #51
	extr	x11, x12, x11, #51
	adds	x11, x23, x11
	adc	x12, x22, x15
	and	x15, x11, #0x7ffffffffffff
	lsr	x16, x12, #51
	extr	x11, x12, x11, #51
	adds	x8, x8, x11
	adc	x9, x9, x16
	lsr	x11, x9, #51
	extr	x9, x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umulh	x12, x9, x20
	umull	x11, w11, w20
	mul	x9, x9, x20
	adds	x9, x9, x10
	adc	x10, x12, x11
	and	x11, x9, #0x7ffffffffffff
	extr	x9, x10, x9, #51
	add	x9, x9, x13
	and	x10, x9, #0x7ffffffffffff
	add	x9, x14, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x15, x9, lsr #51
	and	x13, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w20, x11
	and	x11, x9, #0x7ffffffffffff
	add	x9, x10, x9, lsr #51
	add	x10, x12, x9, lsr #51
	add	x12, x13, x10, lsr #51
	add	x8, x8, x12, lsr #51
	lsr	x13, x8, #51
	umaddl	x11, w13, w20, x11
	and	x9, x9, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	and	x12, x12, #0x7ffffffffffff
	and	x8, x8, #0x7ffffffffffff
	umulh	x13, x11, x11
	mul	x14, x11, x11
	umulh	x15, x8, x9
	mul	x16, x8, x9
	umulh	x17, x12, x10
	mul	x0, x12, x10
	adds	x16, x16, x0
	adc	x15, x15, x17
	umulh	x17, x16, x19
	madd	x15, x15, x19, x17
	mul	x16, x16, x19
	adds	x14, x16, x14
	adc	x13, x15, x13
	lsl	x15, x11, #1
	umulh	x16, x15, x9
	mul	x17, x15, x9
	mul	x0, x8, x19
	umulh	x1, x0, x10
	mul	x2, x0, x10
	mul	x3, x12, x12
	umulh	x4, x12, x12
	umulh	x5, x3, x20
	madd	x4, x4, x20, x5
	mul	x3, x3, x20
	umulh	x5, x15, x10
	mul	x15, x15, x10
	umulh	x6, x9, x9
	mul	x7, x9, x9
	umulh	x21, x0, x12
	umulh	x22, x11, x12
	mul	x23, x11, x12
	umulh	x24, x10, x9
	mul	x25, x10, x9
	adds	x23, x23, x25
	adc	x22, x22, x24
	extr	x22, x22, x23, #63
	mul	x0, x0, x12
	mul	x24, x8, x8
	umulh	x25, x8, x8
	umulh	x26, x24, x20
	madd	x25, x25, x20, x26
	mul	x24, x24, x20
	adds	x23, x24, x23, lsl #1
	adc	x22, x22, x25
	umulh	x24, x11, x8
	mul	x8, x11, x8
	umulh	x11, x12, x9
	mul	x9, x12, x9
	adds	x8, x8, x9
	adc	x9, x24, x11
	extr	x9, x9, x8, #63
	umulh	x11, x10, x10
	mul	x10, x10, x10
	adds	x8, x10, x8, lsl #1
	adc	x9, x9, x11
	and	x10, x14, #0x7ffffffffffff
	lsr	x11, x13, #51
	extr	x12, x13, x14, #51
	adds	x13, x2, x3
	adc	x14, x1, x4
	adds	x13, x13, x17
	adc	x14, x14, x16
	adds	x12, x13, x12
	adc	x11, x14, x11
	and	x13, x12, #0x7ffffffffffff
	lsr	x14, x11, #51
	extr	x11, x11, x12, #51
	adds	x12, x0, x7
	adc	x16, x21, x6
	adds	x12, x12, x15
	adc	x15, x16, x5
	adds	x11, x12, x11
	adc	x12, x15, x14
	and	x14, x11, #0x7ffffffffffff
	lsr	x15, x12, #51
	extr	x11, x12, x11, #51
	adds	x11, x23, x11
	adc	x12, x22, x15
	and	x15, x11, #0x7ffffffffffff
	lsr	x16, x12, #51
	extr	x11, x12, x11, #51
	adds	x8, x8, x11
	adc	x9, x9, x16
	extr	x11, x9, x8, #51
	lsr	x9, x9, #51
	and	x8, x8, #0x7ffffffffffff
	umulh	x12, x11, x20
	umull	x9, w9, w20
	mul	x11, x11, x20
	adds	x10, x11, x10
	adc	x9, x12, x9
	and	x11, x10, #0x7ffffffffffff
	extr	x9, x9, x10, #51
	add	x9, x9, x13
	and	x10, x9, #0x7ffffffffffff
	add	x9, x14, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x15, x9, lsr #51
	and	x13, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w20, x11
	and	x11, x9, #0x7ffffffffffff
	add	x9, x10, x9, lsr #51
	add	x10, x12, x9, lsr #51
	add	x12, x13, x10, lsr #51
	add	x8, x8, x12, lsr #51
	lsr	x13, x8, #51
	umaddl	x11, w13, w20, x11
	and	x9, x9, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	and	x12, x12, #0x7ffffffffffff
	stp	x11, x9, [x29, #-144]
	and	x8, x8, #0x7ffffffffffff
	stp	x10, x12, [x29, #-128]
	stur	x8, [x29, #-112]
	add	x0, sp, #400
	sub	x1, x29, #144
	add	x2, sp, #272
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	add	x0, sp, #440
	add	x1, sp, #400
	add	x2, sp, #360
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldp	x10, x9, [sp, #440]
	ldp	x8, x11, [sp, #456]
	ldr	x12, [sp, #472]
	umulh	x13, x10, x10
	mul	x14, x10, x10
	umulh	x15, x12, x9
	mul	x16, x12, x9
	umulh	x17, x11, x8
	mul	x0, x11, x8
	adds	x16, x16, x0
	adc	x15, x15, x17
	umulh	x17, x16, x19
	madd	x15, x15, x19, x17
	mul	x16, x16, x19
	adds	x14, x16, x14
	adc	x13, x15, x13
	lsl	x15, x10, #1
	lsr	x16, x10, #63
	umulh	x17, x15, x9
	madd	x17, x16, x9, x17
	mul	x0, x15, x9
	mul	x1, x12, x19
	umulh	x2, x12, x19
	umulh	x3, x1, x8
	madd	x3, x2, x8, x3
	mul	x4, x1, x8
	mul	x5, x11, x11
	umulh	x6, x11, x11
	umulh	x7, x5, x20
	madd	x6, x6, x20, x7
	mul	x5, x5, x20
	umulh	x7, x15, x8
	madd	x16, x16, x8, x7
	mul	x15, x15, x8
	umulh	x7, x9, x9
	mul	x21, x9, x9
	adds	x15, x15, x21
	adc	x16, x16, x7
	umulh	x7, x1, x11
	madd	x2, x2, x11, x7
	mul	x1, x1, x11
	adds	x15, x15, x1
	adc	x16, x16, x2
	umulh	x1, x11, x10
	mul	x2, x11, x10
	umulh	x7, x8, x9
	mul	x21, x8, x9
	adds	x2, x2, x21
	adc	x1, x1, x7
	extr	x1, x1, x2, #63
	mul	x7, x12, x12
	umulh	x21, x12, x12
	umulh	x22, x7, x20
	madd	x21, x21, x20, x22
	mul	x7, x7, x20
	adds	x2, x7, x2, lsl #1
	adc	x1, x1, x21
	umulh	x7, x12, x10
	mul	x10, x12, x10
	umulh	x12, x11, x9
	mul	x9, x11, x9
	adds	x9, x10, x9
	adc	x10, x7, x12
	extr	x10, x10, x9, #63
	umulh	x11, x8, x8
	mul	x8, x8, x8
	adds	x8, x8, x9, lsl #1
	adc	x9, x10, x11
	and	x10, x14, #0x7ffffffffffff
	lsr	x11, x13, #51
	extr	x12, x13, x14, #51
	adds	x13, x5, x0
	adc	x14, x6, x17
	adds	x13, x13, x4
	adc	x14, x14, x3
	adds	x12, x13, x12
	adc	x11, x14, x11
	and	x13, x12, #0x7ffffffffffff
	lsr	x14, x11, #51
	extr	x11, x11, x12, #51
	adds	x11, x15, x11
	adc	x12, x16, x14
	and	x14, x11, #0x7ffffffffffff
	lsr	x15, x12, #51
	extr	x11, x12, x11, #51
	adds	x11, x2, x11
	adc	x12, x1, x15
	and	x15, x11, #0x7ffffffffffff
	lsr	x16, x12, #51
	extr	x11, x12, x11, #51
	adds	x8, x8, x11
	adc	x9, x9, x16
	extr	x11, x9, x8, #51
	lsr	x9, x9, #51
	and	x8, x8, #0x7ffffffffffff
	umulh	x12, x11, x20
	umull	x9, w9, w20
	mul	x11, x11, x20
	adds	x10, x11, x10
	adc	x9, x12, x9
	and	x11, x10, #0x7ffffffffffff
	extr	x9, x9, x10, #51
	add	x9, x9, x13
	and	x10, x9, #0x7ffffffffffff
	add	x9, x14, x9, lsr #51
	and	x12, x9, #0x7ffffffffffff
	add	x9, x15, x9, lsr #51
	and	x13, x9, #0x7ffffffffffff
	add	x8, x8, x9, lsr #51
	lsr	x9, x8, #51
	and	x8, x8, #0x7ffffffffffff
	umaddl	x9, w9, w20, x11
	and	x11, x9, #0x7ffffffffffff
	add	x9, x10, x9, lsr #51
	add	x10, x12, x9, lsr #51
	add	x12, x13, x10, lsr #51
	add	x8, x8, x12, lsr #51
	lsr	x13, x8, #51
	umaddl	x11, w13, w20, x11
	and	x9, x9, #0x7ffffffffffff
	and	x10, x10, #0x7ffffffffffff
	and	x12, x12, #0x7ffffffffffff
	stp	x11, x9, [x29, #-144]
	and	x8, x8, #0x7ffffffffffff
	stp	x10, x12, [x29, #-128]
	stur	x8, [x29, #-112]
	add	x0, sp, #480
	sub	x1, x29, #144
	add	x2, sp, #400
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldp	x10, x11, [sp, #480]
	ldp	x9, x12, [sp, #496]
	mov	w8, #5
	ldr	x13, [sp, #512]
LBB1_3:
	umulh	x14, x10, x10
	mul	x15, x10, x10
	umulh	x16, x13, x11
	mul	x17, x13, x11
	umulh	x0, x12, x9
	mul	x1, x12, x9
	adds	x17, x17, x1
	adc	x16, x16, x0
	umulh	x0, x17, x19
	madd	x16, x16, x19, x0
	mul	x17, x17, x19
	adds	x15, x17, x15
	adc	x14, x16, x14
	lsl	x16, x10, #1
	lsr	x17, x10, #63
	umulh	x0, x16, x11
	madd	x0, x17, x11, x0
	mul	x1, x16, x11
	mul	x2, x13, x19
	umulh	x3, x13, x19
	umulh	x4, x2, x9
	madd	x4, x3, x9, x4
	mul	x5, x2, x9
	mul	x6, x12, x12
	umulh	x7, x12, x12
	umulh	x21, x6, x20
	madd	x7, x7, x20, x21
	mul	x6, x6, x20
	umulh	x21, x16, x9
	madd	x17, x17, x9, x21
	mul	x16, x16, x9
	umulh	x21, x11, x11
	mul	x22, x11, x11
	umulh	x23, x2, x12
	madd	x3, x3, x12, x23
	mul	x2, x2, x12
	umulh	x23, x12, x10
	mul	x24, x12, x10
	umulh	x25, x9, x11
	mul	x26, x9, x11
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x13, x13
	umulh	x26, x13, x13
	umulh	x27, x25, x20
	madd	x26, x26, x20, x27
	mul	x25, x25, x20
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x13, x10
	mul	x10, x13, x10
	umulh	x13, x12, x11
	mul	x11, x12, x11
	adds	x10, x10, x11
	adc	x11, x25, x13
	extr	x11, x11, x10, #63
	umulh	x12, x9, x9
	mul	x9, x9, x9
	adds	x9, x9, x10, lsl #1
	adc	x10, x11, x12
	and	x11, x15, #0x7ffffffffffff
	lsr	x12, x14, #51
	extr	x13, x14, x15, #51
	adds	x14, x6, x1
	adc	x15, x7, x0
	adds	x14, x14, x5
	adc	x15, x15, x4
	adds	x13, x14, x13
	adc	x12, x15, x12
	and	x14, x13, #0x7ffffffffffff
	lsr	x15, x12, #51
	extr	x12, x12, x13, #51
	adds	x13, x16, x22
	adc	x16, x17, x21
	adds	x13, x13, x2
	adc	x16, x16, x3
	adds	x12, x13, x12
	adc	x13, x16, x15
	and	x15, x12, #0x7ffffffffffff
	lsr	x16, x13, #51
	extr	x12, x13, x12, #51
	adds	x12, x24, x12
	adc	x13, x23, x16
	and	x16, x12, #0x7ffffffffffff
	lsr	x17, x13, #51
	extr	x12, x13, x12, #51
	adds	x9, x9, x12
	adc	x10, x10, x17
	lsr	x12, x10, #51
	extr	x10, x10, x9, #51
	and	x9, x9, #0x7ffffffffffff
	umulh	x13, x10, x20
	umull	x12, w12, w20
	mul	x10, x10, x20
	adds	x10, x10, x11
	adc	x11, x13, x12
	and	x12, x10, #0x7ffffffffffff
	extr	x10, x11, x10, #51
	add	x10, x10, x14
	and	x11, x10, #0x7ffffffffffff
	add	x10, x15, x10, lsr #51
	and	x13, x10, #0x7ffffffffffff
	add	x10, x16, x10, lsr #51
	and	x14, x10, #0x7ffffffffffff
	add	x9, x9, x10, lsr #51
	lsr	x10, x9, #51
	and	x9, x9, #0x7ffffffffffff
	umaddl	x10, w10, w20, x12
	and	x12, x10, #0x7ffffffffffff
	add	x11, x11, x10, lsr #51
	add	x13, x13, x11, lsr #51
	add	x14, x14, x13, lsr #51
	add	x15, x9, x14, lsr #51
	lsr	x9, x15, #51
	umaddl	x10, w9, w20, x12
	and	x11, x11, #0x7ffffffffffff
	and	x9, x13, #0x7ffffffffffff
	and	x12, x14, #0x7ffffffffffff
	and	x13, x15, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_3
	stp	x10, x11, [x29, #-144]
	stp	x9, x12, [x29, #-128]
	stur	x13, [x29, #-112]
	add	x0, sp, #520
	sub	x1, x29, #144
	add	x2, sp, #480
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x12, [sp, #520]
	ldr	x13, [sp, #528]
	ldr	x11, [sp, #536]
	ldr	x14, [sp, #544]
	mov	w8, #10
	mov	w9, #38
	mov	w10, #19
	ldr	x15, [sp, #552]
LBB1_5:
	umulh	x16, x12, x12
	mul	x17, x12, x12
	umulh	x0, x15, x13
	mul	x1, x15, x13
	umulh	x2, x14, x11
	mul	x3, x14, x11
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x12, #1
	lsr	x1, x12, #63
	umulh	x2, x0, x13
	madd	x2, x1, x13, x2
	mul	x3, x0, x13
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x13, x13
	mul	x22, x13, x13
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x14, x12
	mul	x24, x14, x12
	umulh	x25, x11, x13
	mul	x26, x11, x13
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x15, x12
	mul	x12, x15, x12
	umulh	x15, x14, x13
	mul	x13, x14, x13
	adds	x12, x12, x13
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x19, x3
	adc	x17, x20, x2
	adds	x16, x16, x7
	adc	x17, x17, x6
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x0, x22
	adc	x0, x1, x21
	adds	x15, x15, x4
	adc	x0, x0, x5
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x13, x13, x12, lsr #51
	add	x15, x15, x13, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x12, w11, w10, x14
	and	x13, x13, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_5
	stp	x12, x13, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	add	x0, sp, #560
	sub	x1, x29, #144
	add	x2, sp, #520
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x12, [sp, #560]
	ldr	x13, [sp, #568]
	ldr	x11, [sp, #576]
	ldr	x14, [sp, #584]
	mov	w8, #20
	mov	w9, #38
	mov	w10, #19
	ldr	x15, [sp, #592]
LBB1_7:
	umulh	x16, x12, x12
	mul	x17, x12, x12
	umulh	x0, x15, x13
	mul	x1, x15, x13
	umulh	x2, x14, x11
	mul	x3, x14, x11
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x12, #1
	lsr	x1, x12, #63
	umulh	x2, x0, x13
	madd	x2, x1, x13, x2
	mul	x3, x0, x13
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x13, x13
	mul	x22, x13, x13
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x14, x12
	mul	x24, x14, x12
	umulh	x25, x11, x13
	mul	x26, x11, x13
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x15, x12
	mul	x12, x15, x12
	umulh	x15, x14, x13
	mul	x13, x14, x13
	adds	x12, x12, x13
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x19, x3
	adc	x17, x20, x2
	adds	x16, x16, x7
	adc	x17, x17, x6
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x0, x22
	adc	x0, x1, x21
	adds	x15, x15, x4
	adc	x0, x0, x5
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x13, x13, x12, lsr #51
	add	x15, x15, x13, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x12, w11, w10, x14
	and	x13, x13, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_7
	stp	x12, x13, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	add	x0, sp, #600
	sub	x1, x29, #144
	add	x2, sp, #560
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x12, [sp, #600]
	ldr	x13, [sp, #608]
	ldr	x11, [sp, #616]
	ldr	x14, [sp, #624]
	mov	w8, #10
	mov	w9, #38
	mov	w10, #19
	ldr	x15, [sp, #632]
LBB1_9:
	umulh	x16, x12, x12
	mul	x17, x12, x12
	umulh	x0, x15, x13
	mul	x1, x15, x13
	umulh	x2, x14, x11
	mul	x3, x14, x11
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x12, #1
	lsr	x1, x12, #63
	umulh	x2, x0, x13
	madd	x2, x1, x13, x2
	mul	x3, x0, x13
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x13, x13
	mul	x22, x13, x13
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x14, x12
	mul	x24, x14, x12
	umulh	x25, x11, x13
	mul	x26, x11, x13
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x15, x12
	mul	x12, x15, x12
	umulh	x15, x14, x13
	mul	x13, x14, x13
	adds	x12, x12, x13
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x19, x3
	adc	x17, x20, x2
	adds	x16, x16, x7
	adc	x17, x17, x6
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x0, x22
	adc	x0, x1, x21
	adds	x15, x15, x4
	adc	x0, x0, x5
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x13, x13, x12, lsr #51
	add	x15, x15, x13, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x12, w11, w10, x14
	and	x13, x13, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_9
	stp	x12, x13, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	add	x0, sp, #640
	sub	x1, x29, #144
	add	x2, sp, #520
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x13, [sp, #640]
	ldr	x12, [sp, #648]
	ldr	x11, [sp, #656]
	ldr	x14, [sp, #664]
	mov	w8, #50
	mov	w9, #38
	mov	w10, #19
	ldr	x15, [sp, #672]
LBB1_11:
	umulh	x16, x13, x13
	mul	x17, x13, x13
	umulh	x0, x12, x15
	mul	x1, x12, x15
	umulh	x2, x11, x14
	mul	x3, x11, x14
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x13, #1
	lsr	x1, x13, #63
	umulh	x2, x0, x12
	madd	x2, x1, x12, x2
	mul	x3, x0, x12
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x12, x12
	mul	x22, x12, x12
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x13, x14
	mul	x24, x13, x14
	umulh	x25, x12, x11
	mul	x26, x12, x11
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x13, x15
	mul	x13, x13, x15
	umulh	x15, x12, x14
	mul	x12, x12, x14
	adds	x12, x13, x12
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x7, x19
	adc	x17, x6, x20
	adds	x16, x16, x3
	adc	x17, x17, x2
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x22, x4
	adc	x2, x21, x5
	adds	x15, x15, x0
	adc	x0, x2, x1
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x12, x13, x12, lsr #51
	add	x15, x15, x12, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x13, w11, w10, x14
	and	x12, x12, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_11
	stp	x13, x12, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	add	x0, sp, #680
	sub	x1, x29, #144
	add	x2, sp, #640
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldr	x13, [sp, #680]
	ldr	x12, [sp, #688]
	ldr	x11, [sp, #696]
	ldr	x14, [sp, #704]
	mov	w8, #100
	mov	w9, #38
	mov	w10, #19
	ldr	x15, [sp, #712]
LBB1_13:
	umulh	x16, x13, x13
	mul	x17, x13, x13
	umulh	x0, x12, x15
	mul	x1, x12, x15
	umulh	x2, x11, x14
	mul	x3, x11, x14
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x13, #1
	lsr	x1, x13, #63
	umulh	x2, x0, x12
	madd	x2, x1, x12, x2
	mul	x3, x0, x12
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x12, x12
	mul	x22, x12, x12
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x13, x14
	mul	x24, x13, x14
	umulh	x25, x12, x11
	mul	x26, x12, x11
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x13, x15
	mul	x13, x13, x15
	umulh	x15, x12, x14
	mul	x12, x12, x14
	adds	x12, x13, x12
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x7, x19
	adc	x17, x6, x20
	adds	x16, x16, x3
	adc	x17, x17, x2
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x22, x4
	adc	x2, x21, x5
	adds	x15, x15, x0
	adc	x0, x2, x1
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x12, x13, x12, lsr #51
	add	x15, x15, x12, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x13, w11, w10, x14
	and	x12, x12, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_13
	stp	x13, x12, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	sub	x0, x29, #224
	sub	x1, x29, #144
	add	x2, sp, #680
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldp	x13, x12, [x29, #-224]
	ldp	x11, x14, [x29, #-208]
	mov	w8, #50
	mov	w9, #38
	mov	w10, #19
	ldur	x15, [x29, #-192]
LBB1_15:
	umulh	x16, x13, x13
	mul	x17, x13, x13
	umulh	x0, x12, x15
	mul	x1, x12, x15
	umulh	x2, x11, x14
	mul	x3, x11, x14
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x13, #1
	lsr	x1, x13, #63
	umulh	x2, x0, x12
	madd	x2, x1, x12, x2
	mul	x3, x0, x12
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x12, x12
	mul	x22, x12, x12
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x13, x14
	mul	x24, x13, x14
	umulh	x25, x12, x11
	mul	x26, x12, x11
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x13, x15
	mul	x13, x13, x15
	umulh	x15, x12, x14
	mul	x12, x12, x14
	adds	x12, x13, x12
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x7, x19
	adc	x17, x6, x20
	adds	x16, x16, x3
	adc	x17, x17, x2
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x22, x4
	adc	x2, x21, x5
	adds	x15, x15, x0
	adc	x0, x2, x1
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x12, x13, x12, lsr #51
	add	x15, x15, x12, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x13, w11, w10, x14
	and	x12, x12, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_15
	stp	x13, x12, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	sub	x0, x29, #184
	sub	x1, x29, #144
	add	x2, sp, #640
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	ldp	x13, x12, [x29, #-184]
	ldp	x11, x14, [x29, #-168]
	mov	w8, #5
	mov	w9, #38
	mov	w10, #19
	ldur	x15, [x29, #-152]
LBB1_17:
	umulh	x16, x13, x13
	mul	x17, x13, x13
	umulh	x0, x12, x15
	mul	x1, x12, x15
	umulh	x2, x11, x14
	mul	x3, x11, x14
	adds	x1, x1, x3
	adc	x0, x0, x2
	umulh	x2, x1, x9
	madd	x0, x0, x9, x2
	mul	x1, x1, x9
	adds	x17, x1, x17
	adc	x16, x0, x16
	lsl	x0, x13, #1
	lsr	x1, x13, #63
	umulh	x2, x0, x12
	madd	x2, x1, x12, x2
	mul	x3, x0, x12
	mul	x4, x15, x9
	umulh	x5, x15, x9
	umulh	x6, x4, x11
	madd	x6, x5, x11, x6
	mul	x7, x4, x11
	mul	x19, x14, x14
	umulh	x20, x14, x14
	umulh	x21, x19, x10
	madd	x20, x20, x10, x21
	mul	x19, x19, x10
	umulh	x21, x0, x11
	madd	x1, x1, x11, x21
	mul	x0, x0, x11
	umulh	x21, x12, x12
	mul	x22, x12, x12
	umulh	x23, x4, x14
	madd	x5, x5, x14, x23
	mul	x4, x4, x14
	umulh	x23, x13, x14
	mul	x24, x13, x14
	umulh	x25, x12, x11
	mul	x26, x12, x11
	adds	x24, x24, x26
	adc	x23, x23, x25
	extr	x23, x23, x24, #63
	mul	x25, x15, x15
	umulh	x26, x15, x15
	umulh	x27, x25, x10
	madd	x26, x26, x10, x27
	mul	x25, x25, x10
	adds	x24, x25, x24, lsl #1
	adc	x23, x23, x26
	umulh	x25, x13, x15
	mul	x13, x13, x15
	umulh	x15, x12, x14
	mul	x12, x12, x14
	adds	x12, x13, x12
	adc	x13, x25, x15
	extr	x13, x13, x12, #63
	umulh	x14, x11, x11
	mul	x11, x11, x11
	adds	x11, x11, x12, lsl #1
	adc	x12, x13, x14
	and	x13, x17, #0x7ffffffffffff
	lsr	x14, x16, #51
	extr	x15, x16, x17, #51
	adds	x16, x7, x19
	adc	x17, x6, x20
	adds	x16, x16, x3
	adc	x17, x17, x2
	adds	x15, x16, x15
	adc	x14, x17, x14
	and	x16, x15, #0x7ffffffffffff
	lsr	x17, x14, #51
	extr	x14, x14, x15, #51
	adds	x15, x22, x4
	adc	x2, x21, x5
	adds	x15, x15, x0
	adc	x0, x2, x1
	adds	x14, x15, x14
	adc	x15, x0, x17
	and	x17, x14, #0x7ffffffffffff
	lsr	x0, x15, #51
	extr	x14, x15, x14, #51
	adds	x14, x24, x14
	adc	x15, x23, x0
	and	x0, x14, #0x7ffffffffffff
	lsr	x1, x15, #51
	extr	x14, x15, x14, #51
	adds	x11, x11, x14
	adc	x12, x12, x1
	lsr	x14, x12, #51
	extr	x12, x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umulh	x15, x12, x10
	umull	x14, w14, w10
	mul	x12, x12, x10
	adds	x12, x12, x13
	adc	x13, x15, x14
	and	x14, x12, #0x7ffffffffffff
	extr	x12, x13, x12, #51
	add	x12, x12, x16
	and	x13, x12, #0x7ffffffffffff
	add	x12, x17, x12, lsr #51
	and	x15, x12, #0x7ffffffffffff
	add	x12, x0, x12, lsr #51
	and	x16, x12, #0x7ffffffffffff
	add	x11, x11, x12, lsr #51
	lsr	x12, x11, #51
	and	x11, x11, #0x7ffffffffffff
	umaddl	x12, w12, w10, x14
	and	x14, x12, #0x7ffffffffffff
	add	x12, x13, x12, lsr #51
	add	x15, x15, x12, lsr #51
	add	x16, x16, x15, lsr #51
	add	x17, x11, x16, lsr #51
	lsr	x11, x17, #51
	umaddl	x13, w11, w10, x14
	and	x12, x12, #0x7ffffffffffff
	and	x11, x15, #0x7ffffffffffff
	and	x14, x16, #0x7ffffffffffff
	and	x15, x17, #0x7ffffffffffff
	subs	x8, x8, #1
	b.ne	LBB1_17
	stp	x13, x12, [x29, #-144]
	stp	x11, x14, [x29, #-128]
	stur	x15, [x29, #-112]
	add	x0, sp, #320
	sub	x1, x29, #144
	add	x2, sp, #440
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	sub	x0, x29, #144
	add	x1, sp, #224
	add	x2, sp, #320
	bl	__ZN16ecc_x25519_check3ecc6f255192Fe3mul17h4dc1a5ad8a0c6ac5E
	mov	x10, #-19
	movk	x10, #7, lsl #48
	ldp	x8, x9, [x29, #-144]
	cmp	x8, x10
	cset	w11, lo
	mov	x10, #-2251799813685247
	mov	x12, #2251799813685247
	cmp	x9, x12
	add	x13, x9, x10
	sub	x14, x13, x11
	ccmp	x13, x11, #0, hs
	cset	w11, lo
	ldp	x13, x15, [x29, #-128]
	cmp	x13, x12
	add	x16, x13, x10
	sub	x17, x16, x11
	ccmp	x16, x11, #0, hs
	cset	w11, lo
	cmp	x15, x12
	add	x12, x15, x10
	sub	x16, x12, x11
	ccmp	x12, x11, #0, hs
	cset	w11, lo
	ldur	x12, [x29, #-112]
	mov	x0, #2251799813685246
	cmp	x12, x0
	add	x10, x12, x10
	sub	x0, x10, x11
	ccmp	x10, x11, #0, hi
	cset	w10, hs
	sbfx	x10, x10, #0, #1
	stur	x10, [x29, #-184]
	sub	x10, x29, #184
	; InlineAsm Start
	; InlineAsm End
	ldur	x10, [x29, #-184]
	and	x11, x0, x10
	bic	x12, x12, x10
	orr	x11, x11, x12
	lsl	x11, x11, #12
	and	x12, x16, x10
	bic	x15, x15, x10
	orr	x12, x12, x15
	orr	x11, x11, x12, lsr #39
	lsl	x12, x12, #25
	and	x15, x10, x17
	bic	x13, x13, x10
	orr	x13, x15, x13
	lsl	x15, x13, #38
	and	x14, x10, x14
	bic	x9, x9, x10
	orr	x9, x14, x9
	orr	x14, x15, x9, lsr #13
	mov	x15, #19
	movk	x15, #65528, lsl #48
	add	x15, x8, x15
	and	x15, x10, x15
	bic	x8, x8, x10
	orr	x8, x15, x8
	orr	x8, x8, x9, lsl #51
	ldr	x9, [sp]
	stp	x8, x14, [x9]
	orr	x8, x12, x13, lsr #26
	stp	x8, x11, [x9, #16]
	add	sp, sp, #864
	.cfi_def_cfa wsp, 96
	ldp	x29, x30, [sp, #80]
	ldp	x20, x19, [sp, #64]
	ldp	x22, x21, [sp, #48]
	ldp	x24, x23, [sp, #32]
	ldp	x26, x25, [sp, #16]
	ldp	x28, x27, [sp], #96
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	.cfi_restore w21
	.cfi_restore w22
	.cfi_restore w23
	.cfi_restore w24
	.cfi_restore w25
	.cfi_restore w26
	.cfi_restore w27
	.cfi_restore w28
	ret
	.loh AdrpAdd	Lloh0, Lloh1
	.cfi_endproc

	.p2align	2
__ZN16ecc_x25519_check3hex17h0b884dadc6f00174E:
Lfunc_begin0:
	.cfi_startproc
	.cfi_personality 155, _rust_eh_personality
	.cfi_lsda 16, Lexception0
	sub	sp, sp, #208
	.cfi_def_cfa_offset 208
	stp	x28, x27, [sp, #112]
	stp	x26, x25, [sp, #128]
	stp	x24, x23, [sp, #144]
	stp	x22, x21, [sp, #160]
	stp	x20, x19, [sp, #176]
	stp	x29, x30, [sp, #192]
	add	x29, sp, #192
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	.cfi_remember_state
	mov	x20, x1
	mov	x19, x0
	str	x1, [sp, #32]
	add	x8, sp, #32
Lloh2:
	adrp	x25, __ZN45_$LT$$RF$T$u20$as$u20$core..fmt..LowerHex$GT$3fmt17h80359584ca9d4511E@PAGE
Lloh3:
	add	x25, x25, __ZN45_$LT$$RF$T$u20$as$u20$core..fmt..LowerHex$GT$3fmt17h80359584ca9d4511E@PAGEOFF
	stp	x8, x25, [sp, #64]
Lloh4:
	adrp	x0, l_anon.1519422ea68d41d5036578e62b7be619.1@PAGE
Lloh5:
	add	x0, x0, l_anon.1519422ea68d41d5036578e62b7be619.1@PAGEOFF
	add	x8, sp, #8
	add	x1, sp, #64
	bl	__ZN5alloc3fmt6format12format_inner17hfd0ba00c746e1ec5E
	ldr	x8, [sp, #8]
	mov	x9, #-9223372036854775808
	cmp	x8, x9
	b.ne	LBB2_2
	mov	w8, #1
	stp	xzr, x8, [x19]
	str	xzr, [x19, #16]
	b	LBB2_11
LBB2_2:
	add	x26, x20, #1
	ldur	q0, [sp, #8]
	str	q0, [sp, #32]
	ldr	x21, [sp, #24]
	str	x21, [sp, #48]
	mov	w27, #31
	add	x28, sp, #88
Lloh6:
	adrp	x20, l_anon.1519422ea68d41d5036578e62b7be619.1@PAGE
Lloh7:
	add	x20, x20, l_anon.1519422ea68d41d5036578e62b7be619.1@PAGEOFF
	b	LBB2_4
LBB2_3:
	add	x26, x26, #1
	subs	x27, x27, #1
	b.eq	LBB2_10
LBB2_4:
	stp	x26, x28, [sp, #88]
	str	x25, [sp, #104]
Ltmp0:
	add	x8, sp, #64
	add	x1, sp, #96
	mov	x0, x20
	bl	__ZN5alloc3fmt6format12format_inner17hfd0ba00c746e1ec5E
Ltmp1:
	ldp	x22, x23, [sp, #64]
	ldr	x24, [sp, #80]
	ldr	x8, [sp, #32]
	sub	x8, x8, x21
	cmp	x24, x8
	b.hi	LBB2_8
LBB2_6:
	ldr	x8, [sp, #40]
	add	x0, x8, x21
	mov	x1, x23
	mov	x2, x24
	bl	_memcpy
	add	x21, x21, x24
	str	x21, [sp, #48]
	cbz	x22, LBB2_3
	mov	x0, x23
	mov	x1, x22
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
	b	LBB2_3
LBB2_8:
Ltmp3:
	add	x0, sp, #32
	mov	x1, x21
	mov	x2, x24
	bl	__ZN5alloc7raw_vec20RawVecInner$LT$A$GT$7reserve21do_reserve_and_handle17h879fecadec90dfa4E
Ltmp4:
	ldr	x21, [sp, #48]
	b	LBB2_6
LBB2_10:
	ldr	q0, [sp, #32]
	str	q0, [x19]
	ldr	x8, [sp, #48]
	str	x8, [x19, #16]
LBB2_11:
	.cfi_def_cfa wsp, 208
	ldp	x29, x30, [sp, #192]
	ldp	x20, x19, [sp, #176]
	ldp	x22, x21, [sp, #160]
	ldp	x24, x23, [sp, #144]
	ldp	x26, x25, [sp, #128]
	ldp	x28, x27, [sp, #112]
	add	sp, sp, #208
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	.cfi_restore w21
	.cfi_restore w22
	.cfi_restore w23
	.cfi_restore w24
	.cfi_restore w25
	.cfi_restore w26
	.cfi_restore w27
	.cfi_restore w28
	ret
LBB2_12:
	.cfi_restore_state
Ltmp5:
	mov	x19, x0
	cbnz	x22, LBB2_15
	ldr	x1, [sp, #32]
	cbnz	x1, LBB2_17
LBB2_14:
	mov	x0, x19
	bl	__Unwind_Resume
LBB2_15:
	mov	x0, x23
	mov	x1, x22
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
	ldr	x1, [sp, #32]
	cbz	x1, LBB2_14
	b	LBB2_17
LBB2_16:
Ltmp2:
	mov	x19, x0
	ldr	x1, [sp, #32]
	cbz	x1, LBB2_14
LBB2_17:
	ldr	x0, [sp, #40]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
	mov	x0, x19
	bl	__Unwind_Resume
	.loh AdrpAdd	Lloh4, Lloh5
	.loh AdrpAdd	Lloh2, Lloh3
	.loh AdrpAdd	Lloh6, Lloh7
Lfunc_end0:
	.cfi_endproc
	.section	__TEXT,__gcc_except_tab
	.p2align	2, 0x0
GCC_except_table2:
Lexception0:
	.byte	255
	.byte	255
	.byte	1
	.uleb128 Lcst_end0-Lcst_begin0
Lcst_begin0:
	.uleb128 Lfunc_begin0-Lfunc_begin0
	.uleb128 Ltmp0-Lfunc_begin0
	.byte	0
	.byte	0
	.uleb128 Ltmp0-Lfunc_begin0
	.uleb128 Ltmp1-Ltmp0
	.uleb128 Ltmp2-Lfunc_begin0
	.byte	0
	.uleb128 Ltmp1-Lfunc_begin0
	.uleb128 Ltmp3-Ltmp1
	.byte	0
	.byte	0
	.uleb128 Ltmp3-Lfunc_begin0
	.uleb128 Ltmp4-Ltmp3
	.uleb128 Ltmp5-Lfunc_begin0
	.byte	0
	.uleb128 Ltmp4-Lfunc_begin0
	.uleb128 Lfunc_end0-Ltmp4
	.byte	0
	.byte	0
Lcst_end0:
	.p2align	2, 0x0

	.section	__TEXT,__text,regular,pure_instructions
	.private_extern	__ZN16ecc_x25519_check4main17h62f82c4e5abb7c1dE
	.globl	__ZN16ecc_x25519_check4main17h62f82c4e5abb7c1dE
	.p2align	2
__ZN16ecc_x25519_check4main17h62f82c4e5abb7c1dE:
Lfunc_begin1:
	.cfi_startproc
	.cfi_personality 155, _rust_eh_personality
	.cfi_lsda 16, Lexception1
	sub	sp, sp, #304
	.cfi_def_cfa_offset 304
	stp	x26, x25, [sp, #224]
	stp	x24, x23, [sp, #240]
	stp	x22, x21, [sp, #256]
	stp	x20, x19, [sp, #272]
	stp	x29, x30, [sp, #288]
	add	x29, sp, #288
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_remember_state
Lloh8:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.2@PAGE
Lloh9:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.2@PAGEOFF
	add	x0, sp, #40
	mov	w2, #64
	bl	__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E
Lloh10:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.3@PAGE
Lloh11:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.3@PAGEOFF
	add	x0, sp, #144
	mov	w2, #64
	bl	__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E
	add	x0, sp, #8
	add	x1, sp, #40
	add	x2, sp, #144
	bl	__ZN16ecc_x25519_check3ecc6x255196x2551917h51fdf6a36a3c95f0E
Lloh12:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.4@PAGE
Lloh13:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.4@PAGEOFF
	sub	x0, x29, #112
	mov	w2, #64
	bl	__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E
	ldp	x8, x9, [sp, #8]
	ldp	x10, x11, [x29, #-112]
	ldp	x12, x13, [sp, #24]
	ldp	x14, x15, [x29, #-96]
	cmp	x8, x10
	ccmp	x9, x11, #0, eq
	ccmp	x12, x14, #0, eq
	ccmp	x13, x15, #0, eq
	b.ne	LBB3_15
Lloh14:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.5@PAGE
Lloh15:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.5@PAGEOFF
	add	x0, sp, #40
	mov	w2, #64
	bl	__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E
Lloh16:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.6@PAGE
Lloh17:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.6@PAGEOFF
	add	x0, sp, #144
	mov	w2, #64
	bl	__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E
	add	x0, sp, #8
	add	x1, sp, #40
	add	x2, sp, #144
	bl	__ZN16ecc_x25519_check3ecc6x255196x2551917h51fdf6a36a3c95f0E
Lloh18:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.7@PAGE
Lloh19:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.7@PAGEOFF
	sub	x0, x29, #112
	mov	w2, #64
	bl	__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E
	ldp	x8, x9, [sp, #8]
	ldp	x10, x11, [x29, #-112]
	ldp	x12, x13, [sp, #24]
	ldp	x14, x15, [x29, #-96]
	cmp	x8, x10
	ccmp	x9, x11, #0, eq
	ccmp	x12, x14, #0, eq
	ccmp	x13, x15, #0, eq
	b.ne	LBB3_15
	mov	x25, #53493
	movk	x25, #13342, lsl #16
	movk	x25, #9852, lsl #32
	movk	x25, #35745, lsl #48
	mov	w20, #100
	add	x21, sp, #72
	add	x22, sp, #96
	add	x23, sp, #120
Lloh20:
	adrp	x24, __ZN60_$LT$alloc..string..String$u20$as$u20$core..fmt..Display$GT$3fmt17hf54eb1e03ec5d969E@PAGE
Lloh21:
	add	x24, x24, __ZN60_$LT$alloc..string..String$u20$as$u20$core..fmt..Display$GT$3fmt17hf54eb1e03ec5d969E@PAGEOFF
Lloh22:
	adrp	x19, l_anon.1519422ea68d41d5036578e62b7be619.8@PAGE
Lloh23:
	add	x19, x19, l_anon.1519422ea68d41d5036578e62b7be619.8@PAGEOFF
	b	LBB3_4
LBB3_3:
	subs	w20, w20, #1
	b.eq	LBB3_14
LBB3_4:
	eor	x8, x25, x25, lsl #13
	eor	x8, x8, x8, lsr #7
	eor	x8, x8, x8, lsl #17
	eor	x9, x8, x8, lsl #13
	eor	x9, x9, x9, lsr #7
	eor	x9, x9, x9, lsl #17
	stp	x8, x9, [sp, #8]
	eor	x8, x9, x9, lsl #13
	eor	x8, x8, x8, lsr #7
	eor	x8, x8, x8, lsl #17
	eor	x9, x8, x8, lsl #13
	eor	x9, x9, x9, lsr #7
	eor	x9, x9, x9, lsl #17
	stp	x8, x9, [sp, #24]
	eor	x8, x9, x9, lsl #13
	eor	x8, x8, x8, lsr #7
	eor	x8, x8, x8, lsl #17
	eor	x9, x8, x8, lsl #13
	eor	x9, x9, x9, lsr #7
	eor	x9, x9, x9, lsl #17
	stp	x8, x9, [sp, #40]
	eor	x8, x9, x9, lsl #13
	eor	x8, x8, x8, lsr #7
	eor	x8, x8, x8, lsl #17
	eor	x9, x8, x8, lsl #13
	eor	x9, x9, x9, lsr #7
	eor	x25, x9, x9, lsl #17
	stp	x8, x25, [sp, #56]
	add	x0, sp, #72
	add	x1, sp, #8
	bl	__ZN16ecc_x25519_check3hex17h0b884dadc6f00174E
Ltmp6:
	add	x0, sp, #96
	add	x1, sp, #40
	bl	__ZN16ecc_x25519_check3hex17h0b884dadc6f00174E
Ltmp7:
Ltmp9:
	add	x0, sp, #144
	add	x1, sp, #8
	add	x2, sp, #40
	bl	__ZN16ecc_x25519_check3ecc6x255196x2551917h51fdf6a36a3c95f0E
Ltmp10:
Ltmp11:
	add	x0, sp, #120
	add	x1, sp, #144
	bl	__ZN16ecc_x25519_check3hex17h0b884dadc6f00174E
Ltmp12:
	stp	x21, x24, [x29, #-112]
	stp	x22, x24, [x29, #-96]
	stp	x23, x24, [x29, #-80]
Ltmp14:
	sub	x1, x29, #112
	mov	x0, x19
	bl	__ZN3std2io5stdio6_print17h82ffd9d3f3a616e2E
Ltmp15:
	ldr	x1, [sp, #120]
	cbz	x1, LBB3_10
	ldr	x0, [sp, #128]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
LBB3_10:
	ldr	x1, [sp, #96]
	cbz	x1, LBB3_12
	ldr	x0, [sp, #104]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
LBB3_12:
	ldr	x1, [sp, #72]
	cbz	x1, LBB3_3
	ldr	x0, [sp, #80]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
	b	LBB3_3
LBB3_14:
	.cfi_def_cfa wsp, 304
	ldp	x29, x30, [sp, #288]
	ldp	x20, x19, [sp, #272]
	ldp	x22, x21, [sp, #256]
	ldp	x24, x23, [sp, #240]
	ldp	x26, x25, [sp, #224]
	add	sp, sp, #304
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	.cfi_restore w21
	.cfi_restore w22
	.cfi_restore w23
	.cfi_restore w24
	.cfi_restore w25
	.cfi_restore w26
	ret
LBB3_15:
	.cfi_restore_state
	add	x0, sp, #8
	sub	x1, x29, #112
	bl	__ZN4core9panicking13assert_failed17hfecf9c4736e337dfE
LBB3_16:
Ltmp8:
	mov	x19, x0
	b	LBB3_22
LBB3_17:
Ltmp16:
	mov	x19, x0
	ldr	x1, [sp, #120]
	cbz	x1, LBB3_20
	ldr	x0, [sp, #128]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
	b	LBB3_20
LBB3_19:
Ltmp13:
	mov	x19, x0
LBB3_20:
	ldr	x1, [sp, #96]
	cbz	x1, LBB3_22
	ldr	x0, [sp, #104]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
LBB3_22:
	ldr	x1, [sp, #72]
	cbz	x1, LBB3_24
	ldr	x0, [sp, #80]
	mov	w2, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_dealloc
LBB3_24:
	mov	x0, x19
	bl	__Unwind_Resume
	.loh AdrpAdd	Lloh12, Lloh13
	.loh AdrpAdd	Lloh10, Lloh11
	.loh AdrpAdd	Lloh8, Lloh9
	.loh AdrpAdd	Lloh18, Lloh19
	.loh AdrpAdd	Lloh16, Lloh17
	.loh AdrpAdd	Lloh14, Lloh15
	.loh AdrpAdd	Lloh22, Lloh23
	.loh AdrpAdd	Lloh20, Lloh21
Lfunc_end1:
	.cfi_endproc
	.section	__TEXT,__gcc_except_tab
	.p2align	2, 0x0
GCC_except_table3:
Lexception1:
	.byte	255
	.byte	255
	.byte	1
	.uleb128 Lcst_end1-Lcst_begin1
Lcst_begin1:
	.uleb128 Lfunc_begin1-Lfunc_begin1
	.uleb128 Ltmp6-Lfunc_begin1
	.byte	0
	.byte	0
	.uleb128 Ltmp6-Lfunc_begin1
	.uleb128 Ltmp7-Ltmp6
	.uleb128 Ltmp8-Lfunc_begin1
	.byte	0
	.uleb128 Ltmp9-Lfunc_begin1
	.uleb128 Ltmp12-Ltmp9
	.uleb128 Ltmp13-Lfunc_begin1
	.byte	0
	.uleb128 Ltmp14-Lfunc_begin1
	.uleb128 Ltmp15-Ltmp14
	.uleb128 Ltmp16-Lfunc_begin1
	.byte	0
	.uleb128 Ltmp15-Lfunc_begin1
	.uleb128 Lfunc_end1-Ltmp15
	.byte	0
	.byte	0
Lcst_end1:
	.p2align	2, 0x0

	.section	__TEXT,__text,regular,pure_instructions
	.p2align	2
__ZN16ecc_x25519_check5hex3217hcdef10b643f8ca75E:
	.cfi_startproc
	sub	sp, sp, #64
	.cfi_def_cfa_offset 64
	stp	x29, x30, [sp, #48]
	add	x29, sp, #48
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_remember_state
	mov	x9, #0
	movi.2d	v0, #0000000000000000
	stp	q0, q0, [sp]
	mov	w10, #1
	mov	x11, sp
	mov	x12, #9223372036854775807
LBB4_1:
	lsl	x8, x9, #1
	add	x3, x8, #2
	cbz	x9, LBB4_6
	cmp	x9, x12
	b.eq	LBB4_15
	cmp	x8, x2
	b.hs	LBB4_5
	ldrsb	w13, [x1, x8]
	cmn	w13, #64
	b.ge	LBB4_6
	b	LBB4_16
LBB4_5:
	b.ne	LBB4_16
LBB4_6:
	cmp	x3, x2
	b.hs	LBB4_8
	ldrsb	w13, [x1, x3]
	cmn	w13, #65
	b.gt	LBB4_9
	b	LBB4_16
LBB4_8:
	b.ne	LBB4_16
LBB4_9:
	mov	w14, #0
	add	x13, x9, #1
	add	x15, x1, x8
	mov	x16, x15
	ldrb	w8, [x16], #1
	cmp	w8, #43
	cinc	x8, x10, ne
	csel	x15, x15, x16, ne
LBB4_10:
	ldrb	w16, [x15], #1
	sub	w17, w16, #65
	and	w17, w17, #0xffffffdf
	add	w17, w17, #10
	sub	w3, w16, #48
	cmp	w16, #57
	csel	w16, w17, w3, hi
	cmp	w16, #15
	b.hi	LBB4_14
	orr	w14, w16, w14, lsl #4
	subs	x8, x8, #1
	b.ne	LBB4_10
	strb	w14, [x11, x9]
	mov	x9, x13
	cmp	x13, #32
	b.ne	LBB4_1
	ldp	q0, q1, [sp]
	stp	q0, q1, [x0]
	.cfi_def_cfa wsp, 64
	ldp	x29, x30, [sp, #48]
	add	sp, sp, #64
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	ret
LBB4_14:
	.cfi_restore_state
	mov	w8, #1
	sturb	w8, [x29, #-1]
Lloh24:
	adrp	x0, l_anon.1519422ea68d41d5036578e62b7be619.16@PAGE
Lloh25:
	add	x0, x0, l_anon.1519422ea68d41d5036578e62b7be619.16@PAGEOFF
Lloh26:
	adrp	x3, l_anon.1519422ea68d41d5036578e62b7be619.15@PAGE
Lloh27:
	add	x3, x3, l_anon.1519422ea68d41d5036578e62b7be619.15@PAGEOFF
Lloh28:
	adrp	x4, l_anon.1519422ea68d41d5036578e62b7be619.12@PAGE
Lloh29:
	add	x4, x4, l_anon.1519422ea68d41d5036578e62b7be619.12@PAGEOFF
	sub	x2, x29, #1
	mov	w1, #43
	bl	__ZN4core6result13unwrap_failed17h68056ace71cfeb36E
LBB4_15:
	mov	x3, #0
	mov	x8, #-2
LBB4_16:
Lloh30:
	adrp	x4, l_anon.1519422ea68d41d5036578e62b7be619.11@PAGE
Lloh31:
	add	x4, x4, l_anon.1519422ea68d41d5036578e62b7be619.11@PAGEOFF
	mov	x0, x1
	mov	x1, x2
	mov	x2, x8
	bl	__ZN4core3str16slice_error_fail17hb778393aee5928e5E
	.loh AdrpAdd	Lloh28, Lloh29
	.loh AdrpAdd	Lloh26, Lloh27
	.loh AdrpAdd	Lloh24, Lloh25
	.loh AdrpAdd	Lloh30, Lloh31
	.cfi_endproc

	.private_extern	__ZN3std2rt10lang_start17h7fe016abdcd2ab37E
	.globl	__ZN3std2rt10lang_start17h7fe016abdcd2ab37E
	.p2align	2
__ZN3std2rt10lang_start17h7fe016abdcd2ab37E:
	.cfi_startproc
	sub	sp, sp, #32
	.cfi_def_cfa_offset 32
	stp	x29, x30, [sp, #16]
	add	x29, sp, #16
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	mov	x4, x3
	mov	x3, x2
	mov	x2, x1
	str	x0, [sp, #8]
Lloh32:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.13@PAGE
Lloh33:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.13@PAGEOFF
	add	x0, sp, #8
	bl	__ZN3std2rt19lang_start_internal17h2359c59178685e12E
	.cfi_def_cfa wsp, 32
	ldp	x29, x30, [sp, #16]
	add	sp, sp, #32
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	ret
	.loh AdrpAdd	Lloh32, Lloh33
	.cfi_endproc

	.p2align	2
__ZN3std2rt10lang_start28_$u7b$$u7b$closure$u7d$$u7d$17hf607e04c5d37f1a4E:
	.cfi_startproc
	stp	x29, x30, [sp, #-16]!
	.cfi_def_cfa_offset 16
	mov	x29, sp
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	ldr	x0, [x0]
	bl	__ZN3std3sys9backtrace28__rust_begin_short_backtrace17h5e48724e633e9ef4E
	mov	w0, #0
	.cfi_def_cfa wsp, 16
	ldp	x29, x30, [sp], #16
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	ret
	.cfi_endproc

	.p2align	2
__ZN3std3sys9backtrace28__rust_begin_short_backtrace17h5e48724e633e9ef4E:
	.cfi_startproc
	stp	x29, x30, [sp, #-16]!
	.cfi_def_cfa_offset 16
	mov	x29, sp
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	blr	x0
	; InlineAsm Start
	; InlineAsm End
	.cfi_def_cfa wsp, 16
	ldp	x29, x30, [sp], #16
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	ret
	.cfi_endproc

	.p2align	2
__ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17ha06d04fc56c10028E:
	.cfi_startproc
	ldr	x0, [x0]
	ldr	w8, [x1, #16]
	tbnz	w8, #25, LBB8_3
	tbnz	w8, #26, LBB8_4
	b	__ZN4core3fmt3num3imp51_$LT$impl$u20$core..fmt..Display$u20$for$u20$u8$GT$3fmt17hfc49c4e53f6424f0E
LBB8_3:
	b	__ZN4core3fmt3num52_$LT$impl$u20$core..fmt..LowerHex$u20$for$u20$u8$GT$3fmt17h524eee85f99e5c22E
LBB8_4:
	b	__ZN4core3fmt3num52_$LT$impl$u20$core..fmt..UpperHex$u20$for$u20$u8$GT$3fmt17hb02454fe5a059e44E
	.cfi_endproc

	.p2align	2
__ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hc84990cc29f17970E:
	.cfi_startproc
	sub	sp, sp, #64
	.cfi_def_cfa_offset 64
	stp	x20, x19, [sp, #32]
	stp	x29, x30, [sp, #48]
	add	x29, sp, #48
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	ldr	x20, [x0]
	add	x8, sp, #8
	mov	x0, x1
	bl	__ZN4core3fmt9Formatter10debug_list17hb7fb9b7fe9766868E
	str	x20, [sp, #24]
Lloh34:
	adrp	x19, l_anon.1519422ea68d41d5036578e62b7be619.14@PAGE
Lloh35:
	add	x19, x19, l_anon.1519422ea68d41d5036578e62b7be619.14@PAGEOFF
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #1
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #2
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #3
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #4
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #5
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #6
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #7
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #8
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #9
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #10
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #11
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #12
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #13
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #14
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #15
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #16
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #17
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #18
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #19
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #20
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #21
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #22
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #23
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #24
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #25
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #26
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #27
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #28
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #29
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #30
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x8, x20, #31
	str	x8, [sp, #24]
	add	x0, sp, #8
	add	x1, sp, #24
	mov	x2, x19
	bl	__ZN4core3fmt8builders9DebugList5entry17h5c1ce299239fd84bE
	add	x0, sp, #8
	bl	__ZN4core3fmt8builders9DebugList6finish17h22a4a78ffc880d42E
	.cfi_def_cfa wsp, 64
	ldp	x29, x30, [sp, #48]
	ldp	x20, x19, [sp, #32]
	add	sp, sp, #64
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	ret
	.loh AdrpAdd	Lloh34, Lloh35
	.cfi_endproc

	.p2align	2
__ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E:
	.cfi_startproc
	mov	x8, x1
	ldr	x9, [x0]
	ldrb	w9, [x9]
Lloh36:
	adrp	x10, l_switch.table._ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E@PAGE
Lloh37:
	add	x10, x10, l_switch.table._ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E@PAGEOFF
	ldr	x2, [x10, x9, lsl #3]
Lloh38:
	adrp	x10, l_switch.table._ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E.56@PAGE
Lloh39:
	add	x10, x10, l_switch.table._ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E.56@PAGEOFF
	ldr	x1, [x10, x9, lsl #3]
	mov	x0, x8
	b	__ZN4core3fmt9Formatter9write_str17h51bae54d1ac15513E
	.loh AdrpAdd	Lloh38, Lloh39
	.loh AdrpAdd	Lloh36, Lloh37
	.cfi_endproc

	.p2align	2
__ZN45_$LT$$RF$T$u20$as$u20$core..fmt..LowerHex$GT$3fmt17h80359584ca9d4511E:
	.cfi_startproc
	ldr	x0, [x0]
	b	__ZN4core3fmt3num52_$LT$impl$u20$core..fmt..LowerHex$u20$for$u20$u8$GT$3fmt17h524eee85f99e5c22E
	.cfi_endproc

	.p2align	2
__ZN4core3ops8function6FnOnce40call_once$u7b$$u7b$vtable.shim$u7d$$u7d$17h8a007958b3710bddE:
	.cfi_startproc
	stp	x29, x30, [sp, #-16]!
	.cfi_def_cfa_offset 16
	mov	x29, sp
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	ldr	x0, [x0]
	bl	__ZN3std3sys9backtrace28__rust_begin_short_backtrace17h5e48724e633e9ef4E
	mov	w0, #0
	.cfi_def_cfa wsp, 16
	ldp	x29, x30, [sp], #16
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	ret
	.cfi_endproc

	.p2align	2
__ZN4core9panicking13assert_failed17hfecf9c4736e337dfE:
	.cfi_startproc
	sub	sp, sp, #32
	stp	x29, x30, [sp, #16]
	add	x29, sp, #16
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	stp	x0, x1, [sp]
Lloh40:
	adrp	x2, l_anon.1519422ea68d41d5036578e62b7be619.17@PAGE
Lloh41:
	add	x2, x2, l_anon.1519422ea68d41d5036578e62b7be619.17@PAGEOFF
Lloh42:
	adrp	x7, l_anon.1519422ea68d41d5036578e62b7be619.10@PAGE
Lloh43:
	add	x7, x7, l_anon.1519422ea68d41d5036578e62b7be619.10@PAGEOFF
	mov	x1, sp
	add	x3, sp, #8
	mov	w0, #0
	mov	x4, x2
	mov	x5, #0
	bl	__ZN4core9panicking19assert_failed_inner17h291f2aa7a9cad7e1E
	.loh AdrpAdd	Lloh42, Lloh43
	.loh AdrpAdd	Lloh40, Lloh41
	.cfi_endproc

	.p2align	2
__ZN5alloc7raw_vec20RawVecInner$LT$A$GT$11finish_grow17h00f29fb6c44d8bb4E:
	.cfi_startproc
	stp	x20, x19, [sp, #-32]!
	.cfi_def_cfa_offset 32
	stp	x29, x30, [sp, #16]
	add	x29, sp, #16
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	mov	x20, x0
	tbz	x3, #63, LBB14_2
	mov	x19, #0
	mov	w8, #1
	mov	w9, #8
	b	LBB14_10
LBB14_2:
	mov	x19, x3
	cbz	x1, LBB14_4
	mov	x0, x2
	mov	w2, #1
	mov	x3, x19
	bl	__RNvCs5A4hkLM3k7S_7___rustc14___rust_realloc
	cbnz	x0, LBB14_8
	b	LBB14_6
LBB14_4:
	cbz	x19, LBB14_7
	bl	__RNvCs5A4hkLM3k7S_7___rustc35___rust_no_alloc_shim_is_unstable_v2
	mov	x0, x19
	mov	w1, #1
	bl	__RNvCs5A4hkLM3k7S_7___rustc12___rust_alloc
	cbnz	x0, LBB14_8
LBB14_6:
	mov	w8, #1
	str	x8, [x20, #8]
	b	LBB14_9
LBB14_7:
	mov	w0, #1
LBB14_8:
	mov	x8, #0
	str	x0, [x20, #8]
LBB14_9:
	mov	w9, #16
LBB14_10:
	str	x19, [x20, x9]
	str	x8, [x20]
	.cfi_def_cfa wsp, 32
	ldp	x29, x30, [sp, #16]
	ldp	x20, x19, [sp], #32
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	ret
	.cfi_endproc

	.p2align	2
__ZN5alloc7raw_vec20RawVecInner$LT$A$GT$7reserve21do_reserve_and_handle17h879fecadec90dfa4E:
	.cfi_startproc
	sub	sp, sp, #64
	.cfi_def_cfa_offset 64
	stp	x20, x19, [sp, #32]
	stp	x29, x30, [sp, #48]
	add	x29, sp, #48
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_remember_state
	adds	x8, x2, x1
	b.hs	LBB15_3
	mov	x19, x0
	ldp	x1, x2, [x0]
	lsl	x9, x1, #1
	cmp	x8, x9
	csel	x8, x8, x9, hi
	mov	w9, #8
	cmp	x8, #8
	csel	x20, x8, x9, hi
	add	x0, sp, #8
	mov	x3, x20
	bl	__ZN5alloc7raw_vec20RawVecInner$LT$A$GT$11finish_grow17h00f29fb6c44d8bb4E
	ldr	w8, [sp, #8]
	cmp	w8, #1
	b.eq	LBB15_4
	ldr	x8, [sp, #16]
	stp	x20, x8, [x19]
	.cfi_def_cfa wsp, 64
	ldp	x29, x30, [sp, #48]
	ldp	x20, x19, [sp, #32]
	add	sp, sp, #64
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	.cfi_restore w19
	.cfi_restore w20
	ret
LBB15_3:
	.cfi_restore_state
	mov	x0, #0
	bl	__ZN5alloc7raw_vec12handle_error17h8d88ebedd284f5bcE
LBB15_4:
	ldp	x0, x1, [sp, #16]
	bl	__ZN5alloc7raw_vec12handle_error17h8d88ebedd284f5bcE
	.cfi_endproc

	.p2align	2
__ZN60_$LT$alloc..string..String$u20$as$u20$core..fmt..Display$GT$3fmt17hf54eb1e03ec5d969E:
	.cfi_startproc
	mov	x2, x1
	ldp	x8, x1, [x0, #8]
	mov	x0, x8
	b	__ZN42_$LT$str$u20$as$u20$core..fmt..Display$GT$3fmt17h09f4205c44ecba39E
	.cfi_endproc

	.p2align	2
__ZN68_$LT$core..num..error..ParseIntError$u20$as$u20$core..fmt..Debug$GT$3fmt17h5c13b96b098ecb49E:
	.cfi_startproc
	sub	sp, sp, #32
	.cfi_def_cfa_offset 32
	stp	x29, x30, [sp, #16]
	add	x29, sp, #16
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	mov	x8, x1
	str	x0, [sp, #8]
Lloh44:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.24@PAGE
Lloh45:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.24@PAGEOFF
Lloh46:
	adrp	x3, l_anon.1519422ea68d41d5036578e62b7be619.25@PAGE
Lloh47:
	add	x3, x3, l_anon.1519422ea68d41d5036578e62b7be619.25@PAGEOFF
Lloh48:
	adrp	x6, l_anon.1519422ea68d41d5036578e62b7be619.23@PAGE
Lloh49:
	add	x6, x6, l_anon.1519422ea68d41d5036578e62b7be619.23@PAGEOFF
	add	x5, sp, #8
	mov	x0, x8
	mov	w2, #13
	mov	w4, #4
	bl	__ZN4core3fmt9Formatter26debug_struct_field1_finish17h21312d90f544c233E
	.cfi_def_cfa wsp, 32
	ldp	x29, x30, [sp, #16]
	add	sp, sp, #32
	.cfi_def_cfa_offset 0
	.cfi_restore w30
	.cfi_restore w29
	ret
	.loh AdrpAdd	Lloh48, Lloh49
	.loh AdrpAdd	Lloh46, Lloh47
	.loh AdrpAdd	Lloh44, Lloh45
	.cfi_endproc

	.globl	_main
	.p2align	2
_main:
	.cfi_startproc
	sub	sp, sp, #32
	stp	x29, x30, [sp, #16]
	add	x29, sp, #16
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	mov	x3, x1
	sxtw	x2, w0
Lloh50:
	adrp	x8, __ZN16ecc_x25519_check4main17h62f82c4e5abb7c1dE@PAGE
Lloh51:
	add	x8, x8, __ZN16ecc_x25519_check4main17h62f82c4e5abb7c1dE@PAGEOFF
	str	x8, [sp, #8]
Lloh52:
	adrp	x1, l_anon.1519422ea68d41d5036578e62b7be619.13@PAGE
Lloh53:
	add	x1, x1, l_anon.1519422ea68d41d5036578e62b7be619.13@PAGEOFF
	add	x0, sp, #8
	mov	w4, #0
	bl	__ZN3std2rt19lang_start_internal17h2359c59178685e12E
	ldp	x29, x30, [sp, #16]
	add	sp, sp, #32
	ret
	.loh AdrpAdd	Lloh52, Lloh53
	.loh AdrpAdd	Lloh50, Lloh51
	.cfi_endproc

	.section	__TEXT,__const
	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.0:
	.asciz	"\001\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000\000"

	.section	__TEXT,__literal8,8byte_literals
l_anon.1519422ea68d41d5036578e62b7be619.1:
	.asciz	"\303 \000\000i\002\000"

	.section	__TEXT,__const
l_anon.1519422ea68d41d5036578e62b7be619.2:
	.ascii	"a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"

l_anon.1519422ea68d41d5036578e62b7be619.3:
	.ascii	"e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c"

l_anon.1519422ea68d41d5036578e62b7be619.4:
	.ascii	"c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"

l_anon.1519422ea68d41d5036578e62b7be619.5:
	.ascii	"4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d"

l_anon.1519422ea68d41d5036578e62b7be619.6:
	.ascii	"e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493"

l_anon.1519422ea68d41d5036578e62b7be619.7:
	.ascii	"95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957"

	.section	__TEXT,__cstring,cstring_literals
l_anon.1519422ea68d41d5036578e62b7be619.8:
	.asciz	"\300\001 \300\001 \300\001\n"

l_anon.1519422ea68d41d5036578e62b7be619.9:
	.asciz	"research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.rs"

	.section	__DATA,__const
	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.10:
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.9
	.asciz	"F\000\000\000\000\000\000\000\035\000\000\000\t\000\000"

	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.11:
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.9
	.asciz	"F\000\000\000\000\000\000\000\n\000\000\000*\000\000"

	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.12:
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.9
	.asciz	"F\000\000\000\000\000\000\000\n\000\000\000B\000\000"

	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.13:
	.asciz	"\000\000\000\000\000\000\000\000\b\000\000\000\000\000\000\000\b\000\000\000\000\000\000"
	.quad	__ZN4core3ops8function6FnOnce40call_once$u7b$$u7b$vtable.shim$u7d$$u7d$17h8a007958b3710bddE
	.quad	__ZN3std2rt10lang_start28_$u7b$$u7b$closure$u7d$$u7d$17hf607e04c5d37f1a4E
	.quad	__ZN3std2rt10lang_start28_$u7b$$u7b$closure$u7d$$u7d$17hf607e04c5d37f1a4E

	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.14:
	.asciz	"\000\000\000\000\000\000\000\000\b\000\000\000\000\000\000\000\b\000\000\000\000\000\000"
	.quad	__ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17ha06d04fc56c10028E

	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.15:
	.asciz	"\000\000\000\000\000\000\000\000\001\000\000\000\000\000\000\000\001\000\000\000\000\000\000"
	.quad	__ZN68_$LT$core..num..error..ParseIntError$u20$as$u20$core..fmt..Debug$GT$3fmt17h5c13b96b098ecb49E

	.section	__TEXT,__const
l_anon.1519422ea68d41d5036578e62b7be619.16:
	.ascii	"called `Result::unwrap()` on an `Err` value"

	.section	__DATA,__const
	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.17:
	.asciz	"\000\000\000\000\000\000\000\000\b\000\000\000\000\000\000\000\b\000\000\000\000\000\000"
	.quad	__ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hc84990cc29f17970E

	.section	__TEXT,__const
l_anon.1519422ea68d41d5036578e62b7be619.18:
	.ascii	"Empty"

l_anon.1519422ea68d41d5036578e62b7be619.19:
	.ascii	"InvalidDigit"

l_anon.1519422ea68d41d5036578e62b7be619.20:
	.ascii	"PosOverflow"

l_anon.1519422ea68d41d5036578e62b7be619.21:
	.ascii	"NegOverflow"

	.section	__TEXT,__literal4,4byte_literals
l_anon.1519422ea68d41d5036578e62b7be619.22:
	.ascii	"Zero"

	.section	__DATA,__const
	.p2align	3, 0x0
l_anon.1519422ea68d41d5036578e62b7be619.23:
	.asciz	"\000\000\000\000\000\000\000\000\b\000\000\000\000\000\000\000\b\000\000\000\000\000\000"
	.quad	__ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E

	.section	__TEXT,__const
l_anon.1519422ea68d41d5036578e62b7be619.24:
	.ascii	"ParseIntError"

	.section	__TEXT,__literal4,4byte_literals
l_anon.1519422ea68d41d5036578e62b7be619.25:
	.ascii	"kind"

	.section	__TEXT,__const
	.p2align	3, 0x0
l_switch.table._ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E:
	.quad	5
	.quad	12
	.quad	11
	.quad	11
	.quad	4

	.section	__DATA,__const
	.p2align	3, 0x0
l_switch.table._ZN42_$LT$$RF$T$u20$as$u20$core..fmt..Debug$GT$3fmt17hcb5085fd0e84bbc8E.56:
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.18
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.19
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.20
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.21
	.quad	l_anon.1519422ea68d41d5036578e62b7be619.22

.subsections_via_symbols
