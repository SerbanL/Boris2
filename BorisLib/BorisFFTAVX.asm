.DATA

;useful constants
signmask   DQ 8000000000000000h		;sign mask for lower double of each pair: xorpd to change sign
		   DQ 0000000000000000h
		   DQ 8000000000000000h		
		   DQ 0000000000000000h

signmask_  DQ 0000000000000000h		;sign mask for upper double of each pair : xorpd to change sign
	       DQ 8000000000000000h
		   DQ 0000000000000000h		
	       DQ 8000000000000000h

_signmask_ DQ 8000000000000000h		;sign mask for both upper and lower doubles : xorpd to change sign
		   DQ 8000000000000000h
		   DQ 8000000000000000h		;sign mask for both upper and lower doubles : xorpd to change sign
		   DQ 8000000000000000h 	

zero	   DQ 0000000000000000h
		   DQ 0000000000000000h
		   DQ 0000000000000000h
		   DQ 0000000000000000h

;pointers, data must be 16-byte aligned (Re, Im)
F		DQ ?
F2		DQ ?

Kx		DQ ?
Ky		DQ ?

;twiddle factors, 16-byte aligned pairs (Re, Im)
cossin	DQ ?

;fft shuffling indices look up table
fftShuf DQ ?
fftShufInv DQ ?

Mx		DQ ?
My		DQ ?
Mz		DQ ?

Heffx	DQ ?
Heffy	DQ ?
Heffz	DQ ?

;dimensions
MAXFFT	DQ ?		

nx		DQ ?
ny		DQ ?
nz		DQ ?

nxnynz	DQ ?

N		DQ ?
M		DQ ?
K		DQ ?

powN	DQ ?
powM	DQ ?
powK	DQ ?

NMK		DQ ?

.CODE

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Call Calculations Functions;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {bool IsAVXEnabled(void);}

PUBLIC IsAVXEnabled

IsAVXEnabled PROC

	xor eax, eax
	cpuid
	cmp eax, 1 ; does CPUID support eax = 1?

	jb not_supported

	mov eax, 1

	cpuid
	and ecx, 018000000h ;check 27 bit (OS uses XSAVE/XRSTOR)
	cmp ecx, 018000000h ;and 28 (AVX supported by CPU)

	jne not_supported

	xor ecx, ecx ;XFEATURE_ENABLED_MASK/XCR0 register number = 0
	xgetbv ; XFEATURE_ENABLED_MASK register is in edx:eax
	and eax, 110b
	cmp eax, 110b ;check the AVX registers restore at context switch

	jne not_supported

	mov eax, 1

	ret

	not_supported:
	xor eax, eax

	ret

IsAVXEnabled ENDP


;;;;;;;;;;;;;;;;;;;;;;extern "C" {void FFTAVX_Radix4_DIT(ReIm *FX, int ldn, int N);}
;
; Do Radix-4 FFT using DIT on Input Fx, Output placed in FX
; Fx not affected, input expected to be contiguous in memory, arranged as (Re, Im) doubles, output likewise
; N FFT length ldn = log2(N)
; Fastcall convention : 
; rcx : FX
; rdx : ldn
; r8 : N

PUBLIC FFTAVX_Radix4_DIT
FFTAVX_Radix4_DIT PROC

	push rbx
	push rsi
	push rdi
	push r12
	push r13
	push r14
	push r15

;<RADIX 2 STEP?>
; rbx : loop counter
; rcx : FX
; rdx : ldn
; r8  : N
; r9  : ldn&1
; r10 : FX + i0 * 32
; r15  : MAXFFT

	mov rax, MAXFFT
	mov r15, rax							;need MAXFFT stored for later

	mov r9, rdx
	and r9, 1								;ldn&1
	jz FFTASM_Radix4_DIT_SkipRadix2

	mov r10, rcx							;FX + i0
	mov rbx, r8
	shr rbx, 1								;N / 2 : loop counter

FFTASM_Radix4_DIT_Radix2Loop:
	
	vmovupd ymm0, [r10]						;FX[i0]
	vmovupd ymm1, [r10 + 32]				;FX[i0 + 1]
	vaddpd  ymm2, ymm0, ymm1				;FX[i0] + FX[(i0 + 1)]
	vsubpd  ymm0, ymm0, ymm1				;FX[i0] - FX[(i0 + 1)]
	vmovupd [r10], ymm2						;store sum
	vmovupd [r10 + 32], ymm0				;store diff

	add r10, 64								;FX + i0 + 2

	dec rbx
	jnz FFTASM_Radix4_DIT_Radix2Loop

;</RADIX 2 STEP?>

	shr r15, 1								;MAXFFT / N2 since N2 needs to be 2 (set below). If we jump past this statement then N2 is 1 and this already holds the right value.

FFTASM_Radix4_DIT_SkipRadix2:

;<RADIX 4 PART>
; rax : N4 * 32
; rbx : N2 * 32
; rcx : FX
; rdx : outer loop counter
; rsi : j * 32
; rdi : r * 32
; r8  : N * 32
; r10 : i0 * 32 - used in inner loop only
; r11 : i1 * 32 - used in inner loop only
; r12 : i2 * 32 - used in inner loop only
; r13 : i3 * 32 - used in inner loop only
; r14 : cossin
; r15  : MAXFFT / N2

	mov r14, cossin

	shl r8, 5								;N * 32
	inc r9									;1 << p, where p = ldn&1 : if we jumped here because p was zero then 1 << p = 1. Otherwise p was 1 and now 1 << p = 2
	mov rbx, r9								;N2
	shl rbx, 5								;N2 * 32
	
	dec rdx
	sub rdx, r9
	shr rdx, 1								
	inc rdx									;(ldn - 1 - 1<<p) / 2 + 1: number of outer loop iterations

	vmovupd	ymm8, signmask					;will need this for multiplication by -1 on the lower doubles of each pair

FFTASM_Radix4_DIT_OuterLoop:

	mov rax, rbx							;N4 = N2
	shl rbx, 2								;N2*=4
	shr r15, 2								;MAXFFT / N2

	;treat case j = 0 separately due to no multiplications
	xor rdi, rdi							;r = 0 start index

	mov r10, rcx							;FX + i0, where i0 = 0 at start
	mov r11, r10
	add r11, rax							;FX + i1, where i1 = i0 + N4
	mov r12, r11
	add r12, rax							;FX + i2, where i2 = i1 + N4
	mov r13, r12
	add r13, rax							;FX + i3, where i3 = i2 + N4

FFTASM_Radix4_DIT_SpecialLoop:

	vmovupd ymm4, [r10]						;FX[i0]
	vmovupd ymm5, [r11]						;FX[i1]	
	vmovupd ymm6, [r12]						;FX[i2]	
	vmovupd ymm7, [r13]						;FX[i3]	

	vsubpd ymm0, ymm4, ymm5					;t1 = FX[i0] - FX[i1]
	vaddpd ymm4, ymm4, ymm5					;t0 = FX[i0] + FX[i1]
	
	vsubpd ymm5, ymm6, ymm7					;t3 = FX[i2] * exp1 - FX[i3]
	vaddpd ymm6, ymm6, ymm7					;t2 = FX[i2] * exp1 + FX[i3]
	vshufpd ymm5, ymm5, ymm5, 5				;Re3, Im3
	vxorpd ymm5, ymm5, ymm8					;Re3, -Im3 : t3 = i*t3 now

	vsubpd ymm7, ymm4, ymm6					;t0 - t2
	vaddpd ymm4, ymm4, ymm6					;t0 + t2
	
	vsubpd ymm6, ymm0, ymm5					;t1 - t3
	vaddpd ymm0, ymm0, ymm5					;t1 + t3

	vmovupd [r10], ymm4
	vmovupd [r11], ymm0
	vmovupd [r12], ymm7
	vmovupd [r13], ymm6	

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl FFTASM_Radix4_DIT_SpecialLoop

	mov rsi, 32								;j = 1 start index
	cmp rsi, rax
	jge FFTASM_Radix4_DIT_NextOuterLoop

FFTASM_Radix4_DIT_MiddleLoop:

	xor rdi, rdi							;r = 0 start index
	
	mov r13, r15
	imul r13, rsi							;(MAXFFT / N2) * j * 32
	mov r12, r14
	add r12, r13							;cossin + (MAXFFT / N2) * j
	
	vmovupd ymm1, [r12]						;exp  = cossin[(MAXFFT / N2) * j]
	vmovddup ymm11, ymm1					;cos1, cos1
	vunpckhpd ymm1, ymm1, ymm1				;sin1, sin1

	add r12, r13							;cossin + (MAXFFT / N2) * 2*j
	vmovupd ymm2, [r12]						;exp2 = cossin[(MAXFFT / N2) * 2*j]
	vmovddup ymm12, ymm2					;cos2, cos2
	vunpckhpd ymm2, ymm2, ymm2				;sin2, sin2

	add r12, r13							;cossin + (MAXFFT / N2) * 3*j
	vmovupd ymm3, [r12]						;exp3 = cossin[(MAXFFT / N2) * 3*j]
	vmovddup ymm13, ymm3					;cos3, cos3
	vunpckhpd ymm3, ymm3, ymm3				;sin3, sin3

	mov r10, rsi
	add r10, rcx							;FX + i0, where i0 = j at start
	mov r11, r10
	add r11, rax							;FX + i1, where i1 = i0 + N4
	mov r12, r11
	add r12, rax							;FX + i2, where i2 = i1 + N4
	mov r13, r12
	add r13, rax							;FX + i3, where i3 = i2 + N4

FFTASM_Radix4_DIT_InnerLoop:

	vmovupd ymm0, [r11]						;FX[i1]	
	vmulpd ymm5, ymm0, ymm12				;Im1*cos2, Re1*cos2
	vmulpd ymm0, ymm0, ymm2					;Im1*sin2, Re1*sin2
	vshufpd ymm0, ymm0, ymm0, 5				;Re1*sin2, Im1*sin2
	vaddsubpd ymm5, ymm5, ymm0				;FX[i1] * exp2

	vmovupd ymm0, [r12]						;FX[i2]	
	vmulpd ymm6, ymm0, ymm11				;Im2*cos1, Re2*cos1
	vmulpd ymm0, ymm0, ymm1					;Im2*sin1, Re2*sin1
	vshufpd ymm0, ymm0, ymm0, 5				;Re2*sin1, Im2*sin1
	vaddsubpd ymm6, ymm6, ymm0				;FX[i2] * exp1

	vmovupd ymm0, [r13]						;FX[i3]	
	vmulpd ymm7, ymm0, ymm13				;Im3*cos3, Re3*cos3
	vmulpd ymm0, ymm0, ymm3					;Im3*sin3, Re3*sin3
	vshufpd ymm0, ymm0, ymm0, 5				;Re3*sin3, Im3*sin3
	vaddsubpd ymm7, ymm7, ymm0				;FX[i3] * exp3

	vmovupd ymm4, [r10]						;FX[i0]

	vsubpd ymm0, ymm4, ymm5					;t1 = FX[i0] - FX[i1] * exp2
	vaddpd ymm4, ymm4, ymm5					;t0 = FX[i0] + FX[i1] * exp2
	
	vsubpd ymm5, ymm6, ymm7					;t3 = FX[i2] * exp1 - FX[i3] * exp3
	vaddpd ymm6, ymm6, ymm7					;t2 = FX[i2] * exp1 + FX[i3] * exp3
	vshufpd ymm5, ymm5, ymm5, 5				;Re3, Im3
	vxorpd ymm5, ymm5, ymm8					;Re3, -Im3 : t3 = i*t3 now

	vsubpd ymm7, ymm4, ymm6					;t0 - t2
	vaddpd ymm4, ymm4, ymm6					;t0 + t2
	
	vsubpd ymm6, ymm0, ymm5					;t1 - t3
	vaddpd ymm0, ymm0, ymm5					;t1 + t3

	vmovupd [r10], ymm4
	vmovupd [r11], ymm0
	vmovupd [r12], ymm7
	vmovupd [r13], ymm6	

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl FFTASM_Radix4_DIT_InnerLoop

	add rsi, 32
	cmp rsi, rax
	jl FFTASM_Radix4_DIT_MiddleLoop

FFTASM_Radix4_DIT_NextOuterLoop:

	dec rdx
	jnz FFTASM_Radix4_DIT_OuterLoop

;</RADIX 4 PART>

FFTASM_Radix4_DIT_return:

	;returning now
	pop r15
	pop r14
	pop r13
	pop r12
	pop rdi
	pop rsi
	pop rbx

	ret

FFTAVX_Radix4_DIT ENDP


;;;;;;;;;;;;;;;;;;;;;;extern "C" {void IFFTAVX_Radix4_DIF(ReIm *FX, int ldn, int N);}
;
; Do Radix-4 IFFT using DIF on Input Fx, Output placed in FX
; input expected to be contiguous in memory, arranged as (Re, Im) doubles, output likewise
; N FFT length ldn = log2(N)
; Fastcall convention : 
; rcx : FX
; rdx : ldn
; r8 : N

PUBLIC IFFTAVX_Radix4_DIF
IFFTAVX_Radix4_DIF PROC

	push rbx
	push rsi
	push rdi
	push r12
	push r13
	push r14
	push r15
		
;<RADIX 4 PART>
; rax : N4 * 32
; rbx : N2 * 32
; rcx : FX
; rdx : outer loop counter
; rsi : j * 32
; rdi : r * 32
; r8  : N * 32
; r9  : ldn
; r10 : i0 * 32 - used in inner loop only
; r11 : i1 * 32 - used in inner loop only
; r12 : i2 * 32 - used in inner loop only
; r13 : i3 * 32 - used in inner loop only
; r14 : cossin
; r15  : MAXFFT / N2

	mov r14, cossin

	mov r9, rdx								;save ldn for last part and to use rdx for div
	
	xor rdx, rdx
	mov rax, MAXFFT
	mov rsi, r8
	div esi									;MAXFFT / N
	mov r15, rax
	shr r15, 2								;MAXFFT / N2 = MAXFFT / (4*N)

	shl r8, 5								;N * 32
	mov rbx, r8								
	shl rbx, 2								;N2 = (4 * N) * 32
	mov rax, r8								;set N4 = N at start

	mov rdx, r9								;get ldn back
	dec rdx
	dec rdx
	shr rdx, 1														
	inc rdx									;(ldn - 2) / 2 + 1: number of outer loop iterations

	vmovupd	ymm8, signmask					;will need this for multiplication by -1 on the lower doubles
	vmovupd	ymm9, signmask_					;will need this for multiplication by -1 on the upper doubles

IFFTASM_Radix4_DIF_OuterLoop:

	shr rbx, 2								;N2 = N2 / 4
	shr rax, 2								;N4 = N2 / 4
	shl r15, 2								;MAXFFT / N2

	;treat j = 0 case separately due to no multiplications
	xor rdi, rdi							;r = 0 start index

	mov r10, rcx							;FX + i0, where i0 = 0 at start
	mov r11, r10
	add r11, rax							;FX + i1, where i1 = i0 + N4
	mov r12, r11
	add r12, rax							;FX + i2, where i2 = i1 + N4
	mov r13, r12
	add r13, rax							;FX + i3, where i3 = i2 + N4

IFFTASM_Radix4_DIF_SpecialLoop:

	vmovupd ymm4, [r10]						;FX[i0]
	vmovupd ymm5, [r11]						;FX[i1]	
	vmovupd ymm6, [r12]						;FX[i2]	
	vmovupd ymm7, [r13]						;FX[i3]

	vaddpd ymm0, ymm4, ymm6					;t0 = FX[i0] + FX[i2]
	vsubpd ymm4, ymm4, ymm6					;t1 = FX[i0] - FX[i2]

	vaddpd ymm6, ymm5, ymm7					;t2 = FX[i1] + FX[i3]
	vsubpd ymm5, ymm5, ymm7					;t3 = FX[i1] - FX[i3]
	vshufpd ymm5, ymm5, ymm5, 5				;Re3, Im3
	vxorpd ymm5, ymm5, ymm8					;Re3, -Im3 : t3 = i*t3 now

	vaddpd ymm7, ymm0, ymm6					;t0 + t2
	vsubpd ymm0, ymm0, ymm6					;t0 - t2	

	vaddpd ymm6, ymm4, ymm5					;t1 + t3
	vsubpd ymm4, ymm4, ymm5					;t1 - t3	

	vmovupd [r10], ymm7						;t0 + t2
	vmovupd [r11], ymm0						;t0 - t2
	vmovupd [r12], ymm4						;t1 - t3
	vmovupd [r13], ymm6						;t1 + t3

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl IFFTASM_Radix4_DIF_SpecialLoop	

	mov rsi, 32								;j = 1 start index
	cmp rsi, rax
	jge IFFTASM_Radix4_DIF_NextOuterLoop

IFFTASM_Radix4_DIF_MiddleLoop:

	xor rdi, rdi							;r = 0 start index
	
	mov r13, r15
	imul r13, rsi							;(MAXFFT / N2) * j * 32
	mov r12, r14
	add r12, r13							;cossin + (MAXFFT / N2) * j
	
	vmovupd ymm1, [r12]						;exp  = cossin[(MAXFFT / N2) * j]
	vxorpd ymm1, ymm1, ymm9					;conjugate exp : change sign of sin, the upper double
	vmovddup ymm11, ymm1					;cos1, cos1
	vunpckhpd ymm1, ymm1, ymm1				;sin1, sin1

	add r12, r13							;cossin + (MAXFFT / N2) * 2*j
	vmovupd ymm2, [r12]						;exp2 = cossin[(MAXFFT / N2) * 2*j]
	vxorpd ymm2, ymm2, ymm9					;conjugate exp : change sign of sin, the upper double
	vmovddup ymm12, ymm2					;cos2, cos2
	vunpckhpd ymm2, ymm2, ymm2				;sin2, sin2

	add r12, r13							;cossin + (MAXFFT / N2) * 3*j
	vmovupd ymm3, [r12]						;exp3 = cossin[(MAXFFT / N2) * 3*j]
	vxorpd ymm3, ymm3, ymm9					;conjugate exp : change sign of sin, the upper double
	vmovddup ymm13, ymm3					;cos3, cos3
	vunpckhpd ymm3, ymm3, ymm3				;sin3, sin3

	mov r10, rsi
	add r10, rcx							;FX + i0, where i0 = j at start
	mov r11, r10
	add r11, rax							;FX + i1, where i1 = i0 + N4
	mov r12, r11
	add r12, rax							;FX + i2, where i2 = i1 + N4
	mov r13, r12
	add r13, rax							;FX + i3, where i3 = i2 + N4

IFFTASM_Radix4_DIF_InnerLoop:

	vmovupd ymm4, [r10]						;FX[i0]
	vmovupd ymm5, [r11]						;FX[i1]	
	vmovupd ymm6, [r12]						;FX[i2]	
	vmovupd ymm7, [r13]						;FX[i3]	

	vaddpd ymm0, ymm4, ymm6					;t0 = FX[i0] + FX[i2]
	vsubpd ymm4, ymm4, ymm6					;t1 = FX[i0] - FX[i2]

	vaddpd ymm6, ymm5, ymm7					;t2 = FX[i1] + FX[i3]
	vsubpd ymm5, ymm5, ymm7					;t3 = FX[i1] - FX[i3]
	vshufpd ymm5, ymm5, ymm5, 5				;Re3, Im3
	vxorpd ymm5, ymm5, ymm8					;Re3, -Im3 : t3 = i*t3 now

	vaddpd ymm7, ymm0, ymm6					;t0 + t2
	vsubpd ymm0, ymm0, ymm6					;t0 - t2	

	vaddpd ymm6, ymm4, ymm5					;t1 + t3
	vsubpd ymm4, ymm4, ymm5					;t1 - t3	

	;form (t0-t2)*exp2
	vmulpd ymm5, ymm0, ymm12				;Im*cos2, Re*cos2
	vmulpd ymm0, ymm0, ymm2					;Im*sin2, Re*sin2
	vshufpd ymm0, ymm0, ymm0, 5				;Re*sin2, Im*sin2
	vaddsubpd ymm5, ymm5, ymm0				;(t0-t2)*exp2		

	;form (t1-t3)*exp
	vmulpd ymm0, ymm4, ymm1					;Im*sin, Re*sin
	vmulpd ymm4, ymm4, ymm11				;Im*cos, Re*cos
	vshufpd ymm0, ymm0, ymm0, 5				;Re*sin, Im*sin
	vaddsubpd ymm4, ymm4, ymm0				;(t1-t3)*exp

	;form (t1+t3)*exp3
	vmulpd ymm0, ymm6, ymm3					;Im*sin3, Re*sin3
	vmulpd ymm6, ymm6, ymm13				;Im*cos3, Re*cos3
	vshufpd ymm0, ymm0, ymm0, 5				;Re*sin3, Im*sin3
	vaddsubpd ymm6, ymm6, ymm0				;(t1+t3)*exp3

	vmovupd [r10], ymm7
	vmovupd [r11], ymm5
	vmovupd [r12], ymm4
	vmovupd [r13], ymm6	

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl IFFTASM_Radix4_DIF_InnerLoop

	add rsi, 32
	cmp rsi, rax
	jl IFFTASM_Radix4_DIF_MiddleLoop

IFFTASM_Radix4_DIF_NextOuterLoop:

	dec rdx
	jnz IFFTASM_Radix4_DIF_OuterLoop

;</RADIX 4 PART>

;<RADIX 2 STEP?>
; rcx : FX
; rbx : loop counter
; r8  : N * 32
; r9  : ldn
; r10 : i0 * 32

	and r9, 1								;ldn&1
	jz IFFTASM_Radix4_DIF_return

	mov r10, rcx							;FX + i0
	mov rbx, r8
	shr rbx, 6								;N / 2 : loop counter

IFFTASM_Radix4_DIF_Radix2Loop:
	
	vmovupd ymm0, [r10]						;FX[i0]
	vmovupd ymm1, [r10 + 32]				;FX[i0 + 1]
	vsubpd  ymm2, ymm0, ymm1				;FX[i0] - FX[(i0 + 1)]
	vaddpd  ymm0, ymm0, ymm1				;FX[i0] + FX[(i0 + 1)]
	vmovupd [r10], ymm0
	vmovupd [r10 + 32], ymm2	

	add r10, 64								;FX + i0 + 2

	dec rbx
	jnz IFFTASM_Radix4_DIF_Radix2Loop

;</RADIX 2 STEP?>

IFFTASM_Radix4_DIF_return:

	;returning now
	pop r15
	pop r14
	pop r13
	pop r12
	pop rdi
	pop rsi
	pop rbx

	ret

IFFTAVX_Radix4_DIF ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void CopyShuffleAVXZeroPad_Rows(int m);}
;
; rcx : m (row)

PUBLIC CopyShuffleAVXZeroPad_Rows

CopyShuffleAVXZeroPad_Rows PROC

	push rbx
	push rsi
	push rdi

;<LOOP>
; 
; rbx : F + m*N (output)
; rcx : (N/2) * 8
; rdx : loop counter
; rdi : Mx + m*nx + p (input)
; rsi : My + m*nx + p (input)
; r9 : fftShuffleInv + N - 1 + p
; r11 : (N/2) * 8 - 16

	mov rdx, N
	
	mov rbx, rdx
	imul rbx, rcx							;m*N
	shl rbx, 4								;m*N * 16
	add rbx, F								;F + m*N

	mov rsi, nx								
	imul rsi, rcx							;m*nx
	shl rsi, 3								;m*nx * 8

	mov rdi, Mx
	add rdi, rsi							;Mx + m*nx
	add rsi, My								;My + m*nx

	mov r9, rdx
	dec r9
	shl r9, 2								;(N - 1) * 4
	add r9, fftShufInv						;fftShuffleInv + N - 1

	mov rcx, rdx
	shl rcx, 2								;(N/2) * 8

	mov r11, rcx
	sub r11, 16								;(N/2) * 8 - 16

	shr rdx, 2								;(N / 2) / 2 : number of loop iterations

	vxorpd ymm0, ymm0, ymm0

CopyShuffleAVXZeroPad_Rows_Loop:

	;get data : 2 points from each Mx, My at offsets 0 and N/2 from current point
	vmovupd xmm1, [rdi]						;Mx b, a (Re)
	vmovupd xmm2, [rsi]						;My b, a (Im)

	add rdi, rcx							; + (N/2)  *8
	add rsi, rcx							; + (N/2)  *8
	vmovupd xmm3, [rdi]						;Mx b2, a2 (Re)
	vmovupd xmm4, [rsi]						;My b2, a2 (Im)

	;re-arrange data ready for writing

	vinsertf128 ymm1, ymm1, xmm3, 1			;Reb2, Rea2, Reb, Rea
	vinsertf128 ymm2, ymm2, xmm4, 1			;Imb2, Ima2, Imb, Ima

	vshufpd ymm6, ymm1, ymm2, 0				;Ima2, Rea2, Ima, Rea
	vshufpd ymm1, ymm1, ymm2, 15			;Imb2, Reb2, Imb, Reb

	;write arranged data
	
	;points a
	xor rax, rax
	mov eax, [r9]
	shl rax, 5								;setup index : fftShuffleInv[N - 1 + p] * 32
	add rax, rbx							;F + m*N + index*2
	vmovupd [rax], ymm6						;write points...
	vmovupd [rax + 32], ymm0				;and zero padding

	;points b
	add r9, 4								;next fftShuffle + N - 1 + p
	xor rax, rax
	mov eax, [r9]
	shl rax, 5								;setup index : fftShuffleInv[N - 1 + p] * 32
	add rax, rbx							;F + m*N + index*2
	vmovupd [rax], ymm1						;write points...
	vmovupd [rax + 32], ymm0				;and zero padding

	add r9, 4								;next fftShuffle + N - 1 + p
	sub rdi, r11							;next Mx + m*nx + 4*p
	sub rsi, r11							;next My + m*nx + 4*p

	dec rdx
	jnz CopyShuffleAVXZeroPad_Rows_Loop

;</LOOP>

	pop rdi
	pop rsi
	pop rbx

	ret

CopyShuffleAVXZeroPad_Rows ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void CopyShuffleAVXZeroPad_Cols(int n);}
;
; rcx : n (col)

PUBLIC CopyShuffleAVXZeroPad_Cols

CopyShuffleAVXZeroPad_Cols PROC

	push rbx
	push rdi

;<LOOP>
; 
; rbx : F2 + n*M (output)
; rcx : N (stride) * 16 * 2
; rdx : loop counter
; rdi : F + 2*n + p * stride (input)
; r9 : fftShuffleInv + M - 1 + p
;

	mov rdx, M
	
	mov rbx, rdx
	imul rbx, rcx							;n*M
	shl rbx, 4								;n*M * 16
	add rbx, F2								;F2 + n*M

	mov rdi, rcx
	shl rdi, 5								;2*n * 16
	add rdi, F								;F + 2*n

	mov r9, rdx
	dec r9
	shl r9, 2								;(M - 1) * 4
	add r9, fftShufInv						;fftShuffleInv + M - 1

	mov rcx, N
	shl rcx, 5								;N * 16 * 2

	shr rdx, 2								;N / 4 : number of loop iterations

	vxorpd ymm0, ymm0, ymm0

CopyShuffleAVXZeroPad_Cols_Loop:

	vmovupd ymm1, [rdi]						;Fx[p * stride + 1], Fx[p * stride]
	vmovupd ymm2, [rdi + 32]				;Fx[p * stride + 3], Fx[p * stride + 2]

	vextractf128 xmm3, ymm1, 1				;Fx[p * stride + 1]
	vinsertf128 ymm1, ymm1, xmm2, 1			;Fx[p * stride + 2], Fx[p * stride]
	vinsertf128 ymm2, ymm2, xmm3, 0			;Fx[p * stride + 3], Fx[p * stride + 1]

	xor rax, rax
	mov eax, [r9]
	shl rax, 5								;setup index : fftShuffleInv[N - 1 + p] * 32
	add rax, rbx							;F2 + n*M + index*2
	
	vmovupd [rax], ymm1						;write Fx[p * stride + 2], Fx[p * stride]
	vmovupd [rax + 32], ymm0				;zero pad FX[index*2 + 3], FX[index*2 + 2]

	add r9, 4								;next fftShuffle + N - 1 + p
	xor rax, rax
	mov eax, [r9]
	shl rax, 5								;setup index : fftShuffleInv[N - 1 + p] * 32
	add rax, rbx							;F2 + n*M + index*2
	
	vmovupd [rax], ymm2						;write Fx[p * stride + 3], Fx[p * stride + 1]
	vmovupd [rax + 32], ymm0				;zero pad FX[index*2 + 3], FX[index*2 + 2]

	add r9, 4								;next fftShuffle + N - 1 + p
	add rdi, rcx							;next F + 2*n + p * stride

	dec rdx
	jnz CopyShuffleAVXZeroPad_Cols_Loop

;</LOOP>

	pop rdi
	pop rbx

	ret

CopyShuffleAVXZeroPad_Cols ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void UnShuffleTruncateAVX_Rows(int m);}
;
; rcx : m (row)

PUBLIC UnShuffleTruncateAVX_Rows

UnShuffleTruncateAVX_Rows PROC

	push rbx
	push rdi
	push rsi

;<LOOP>
; 
; rbx : F2 + 2*m (output)
; rcx : M (stride) * 16
; rdx : loop counter
; rdi : F + m*N + 2*p (input)
; rsi : F + m*N + 2*p + N (input)
; r9 : fftShuffleInv + N - 1 + p
;

	mov rdx, N
	
	mov rdi, rdx
	shl rdi, 4								;N * 16
	mov rsi, rdi							;save N * 16
	imul rdi, rcx							;m*N * 16
	add rsi, rdi							;(m*N + N) * 16
	add rdi, F								;F + m*N
	add rsi, F								;F + m*N + N

	mov rsi, rdx
	shl rsi, 4								;N*16
	add rsi, rdi


	mov rbx, rcx
	shl rbx, 5								;2*m * 16
	add rbx, F2								;F2 + 2*m

	mov r9, rdx
	dec r9
	shl r9, 2								;(N - 1) * 4
	add r9, fftShufInv						;fftShuffleInv + N - 1

	mov rcx, M
	shl rcx, 4								;M * 16

	shr rdx, 2								;N / 4 : number of loop iterations

	vxorpd ymm0, ymm0, ymm0

UnShuffleTruncateAVX_Rows_Loop:

	vmovupd ymm1, [rdi]						;Fx[2 * p + 1], Fx[2 * p]
	vmovupd ymm2, [rsi]						;Fx[2 * p + 1 + N], Fx[2 * p + N]

	vextractf128 xmm3, ymm1, 1				;Fx[2 * p + 1]
	vinsertf128 ymm1, ymm1, xmm2, 1			;Fx[2 * p + N], Fx[2 * p]
	vinsertf128 ymm2, ymm2, xmm3, 0			;Fx[2 * p + 1 + N], Fx[2 * p + 1]

	xor rax, rax
	mov eax, [r9]							;index = fftShuffleInv[N - 1 + p]
	imul rax, rcx							;index * stride * 16
	add rax, rbx							;F2 + 2*m + index*stride
	
	vmovupd [rax], ymm1
	vmovupd [rax + 32], ymm2 

	add r9, 8								;next fftShuffleInv + N - 1 + p
	add rdi, 64								;next F + m*N + 2*p
	add rsi, 64								;next F + m*N + 2*p + N

	dec rdx
	jnz UnShuffleTruncateAVX_Rows_Loop

;</LOOP>

	pop rsi
	pop rdi
	pop rbx

	ret

UnShuffleTruncateAVX_Rows ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void UnShuffleTruncateAVX_Cols(int n);}
;
; rcx : n (col)

PUBLIC UnShuffleTruncateAVX_Cols

UnShuffleTruncateAVX_Cols PROC

	push rbx
	push rdi
	push rsi

;<LOOP>
; 
; rbx : F2 + n*M + 2*p (input)
; rcx : nx (stride) * 8
; rdx : loop counter
; rdi : Heffx + n (output)
; rsi : Heffy + n (output)
; r9 : fftShuffleInv + M - 1 + p
;

	mov rdx, M

	shl rcx, 3								;n*8
	
	mov rdi, Heffx							
	add rdi, rcx							;Heffx + n

	mov rsi, Heffy
	add rsi, rcx							;Heffy + n

	imul rcx, rdx							;n*M * 8
	shl rcx, 1								;n*M * 16

	mov rbx, F2
	add rbx, rcx							;F2 + n*M

	mov r9, rdx
	dec r9
	shl r9, 2								;(N - 1) * 4
	add r9, fftShufInv						;fftShuffleInv + M - 1

	mov rcx, nx
	shl rcx, 3								;nx * 8

	shr rdx, 1								;N / 2 : number of loop iterations

UnShuffleTruncateAVX_Cols_Loop:

	vmovupd ymm1, [rbx]						;Fx[2*p + 1], Fx[2*p]

	vextractf128 xmm2, ymm1, 1				;Im2p1, Re2p1
											;xmm1 has Im2p, Re2p

	vshufpd xmm3, xmm1, xmm2, 0				;Re2p1, Re2p
	vshufpd xmm4, xmm1, xmm2, 3				;Im2p1, Re2p1

	xor rax, rax
	mov eax, [r9]							;index = fftShuffleInv[M - 1 + p]
	imul rax, rcx							;index * stride * 8
	
	mov r10, rdi							
	add r10, rax							;Heffx + n + index*stride

	mov r11, rsi
	add r11, rax							;Heffy + n + index*stride
	
	vmovupd [r10], xmm3
	vmovupd [r11], xmm4 

	add r9, 8								;next fftShuffleInv + M - 1 + p
	add rbx, 64								;next F2 + n*M + 2*p

	dec rdx
	jnz UnShuffleTruncateAVX_Cols_Loop

;</LOOP>

	pop rsi
	pop rdi
	pop rbx

	ret

UnShuffleTruncateAVX_Cols ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Call Initialization Functions;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ZeroUpper(void);}

PUBLIC ZeroUpper
ZeroUpper PROC

	vzeroupper

	ret

ZeroUpper ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitAVXDemag_FFTSpaces(ReIm *F, ReIm *F2);}

PUBLIC InitAVXDemag_FFTSpaces
InitAVXDemag_FFTSpaces PROC

	mov F, rcx
	mov F2, rdx

	ret

InitAVXDemag_FFTSpaces ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitAVXDemag_Kernels(ReIm *Kx, ReIm *Ky);}

PUBLIC InitAVXDemag_Kernels
InitAVXDemag_Kernels PROC

	mov Kx, rcx
	mov Ky, rdx

	ret

InitAVXDemag_Kernels ENDP

;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitAVXMeshSpaces(double *Mx, double *My, double *Mz, double *Heffx, double *Heffy, double *Heffz);}

PUBLIC InitAVXMeshSpaces
InitAVXMeshSpaces PROC

	push rbp
	mov rbp, rsp

	mov Mx, rcx				;Mx
	mov My, rdx				;My
	mov Mz, r8				;Mz
	mov Heffx, r9			;Heffx				
	
	mov rax, [rbp + 48]		;Heffy
	mov Heffy, rax
	mov rax, [rbp + 56]		;Heffz
	mov Heffz, rax

	mov rsp, rbp
	pop rbp

	ret

InitAVXMeshSpaces ENDP


;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitAVXDemag_LUT(ReIm *cossin, int *fftShuffle, int *fftShuffleInv, int _MAXFFT);}

PUBLIC InitAVXDemag_LUT
InitAVXDemag_LUT PROC

	mov cossin, rcx
	mov fftShuf, rdx
	mov fftShufInv, r8
	mov MAXFFT, r9

	ret

InitAVXDemag_LUT ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitAVXDemag_sizes(int nx, int ny, int nz, int N, int M, int K);}

PUBLIC InitAVXDemag_sizes
InitAVXDemag_sizes PROC

	push rbp
	mov rbp, rsp

	mov nx, rcx				;nx
	mov ny, rdx				;ny
	mov rax, r8				;nz
	mov nz, rax				
	
	imul rcx, rdx
	imul rcx, rax			;nx * ny * nz
	mov nxnynz, rcx

	mov rax, r9				;N			
	mov N, rax
	mov rcx, rax			;N

	call GetPow2_
	mov powN, rax

	mov rax, [rbp + 48]		;M
	mov M, rax
	imul rcx, rax			;N*M

	call GetPow2_
	mov powM, rax

	mov rax, [rbp + 56]		;K
	mov K, rax
	imul rcx, rax			;N*M*K
	mov NMK, rcx

	call GetPow2_
	mov powK, rax

	mov rsp, rbp
	pop rbp

	ret

InitAVXDemag_sizes ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;Get log2 of number which is a positive power of 2. Input in rax, output also in rax.
;
GetPow2_ PROC

	push rcx

	xor rcx, rcx

GetPow2_Loop:

	inc rcx
	shr rax, 1
	jnz GetPow2_Loop
	
	dec rcx
	mov rax, rcx

	pop rcx

	ret

GetPow2_ ENDP


_TEXT   ENDS
END