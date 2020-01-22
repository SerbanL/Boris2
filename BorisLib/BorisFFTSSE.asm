.DATA

;useful constants
signmask   DQ 8000000000000000h		;sign mask for lower double : xorpd to change sign
		   DQ 0000000000000000h

signmask_  DQ 0000000000000000h		;sign mask for upper double : xorpd to change sign
	       DQ 8000000000000000h

_signmask_ DQ 8000000000000000h		;sign mask for both upper and lower doubles : xorpd to change sign
		   DQ 8000000000000000h 	

zero	   DQ 0000000000000000h
		   DQ 0000000000000000h

_xmm0	DQ ?
		DQ ?
_xmm1	DQ ?
		DQ ?
_xmm2	DQ ?
		DQ ?
_xmm3	DQ ?
		DQ ?
_xmm4	DQ ?
		DQ ?
_xmm5	DQ ?
		DQ ?
_xmm6	DQ ?
		DQ ?
_xmm7	DQ ?
		DQ ?
_xmm8	DQ ?
		DQ ?
_xmm9	DQ ?
		DQ ?
_xmm10	DQ ?
		DQ ?
_xmm11	DQ ?
		DQ ?
_xmm12	DQ ?
		DQ ?
_xmm13	DQ ?
		DQ ?
_xmm14	DQ ?
		DQ ?
_xmm15	DQ ?
		DQ ?

;pointers, data must be 16-byte aligned (Re, Im)
F		DQ ?
F2		DQ ?
F3		DQ ?
F4		DQ ?
F5		DQ ?
F6		DQ ?

Kx		DQ ?
Ky		DQ ?
Kz		DQ ?
Kxy		DQ ?
Kxz		DQ ?
Kyz		DQ ?

;twiddle factors, 16-byte aligned pairs (Re, Im)
cossin	DQ ?
sincos  DQ ?

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void CopyShuffleZeroPadASM_ReInput(double *Rex, double *Imx, ReIm *FX, int N);}
;
; Copy data in shuffled order to F space and zero pad, ready for DIT FFT
;
; rcx : Rex
; rdx : Imx
; r8 : FX
; r9 : N
;
PUBLIC CopyShuffleZeroPadASM_ReInput
CopyShuffleZeroPadASM_ReInput PROC

	push rsi

;<SHUFFLE>
; rcx : Rex (input)
; rdx : Imx (input)
; rsi : fftShufInv + N - 1 + p
; r8  : FX
; r9  : N/2 (loop counter)
;

	mov rsi, r9
	dec rsi
	shl rsi, 2								;(N - 1) * 4
	add rsi, fftShufInv						;fftShufInv + N - 1

	shr r9, 1								;N/2 - number of loop iterations

	xorpd xmm0, xmm0

CopyShuffleZeroPadASM_ReInput_Loop:

	xor rax, rax
	mov eax, [rsi]							;index = fftShuffleInv[N - 1 + p]

	shl rax, 4								;index * 16
	add rax, r8								;FX + index

	movlpd xmm1, QWORD PTR [rcx]			;Rex[p]
	movhpd xmm1, QWORD PTR [rdx]			;Imx[p]

	movupd [rax], xmm1						;FX[index]
	movupd [rax + 16], xmm0					;FX[index + 1] = ReIm(0, 0) : zero-padding

	add rsi, 4								;next fftShufInv + N - 1 + p
	add rcx, 8								;next Rex + p
	add rdx, 8								;next Imx + p

	dec r9
	jnz CopyShuffleZeroPadASM_ReInput_Loop

;</SHUFFLE>

	pop rsi

	ret

CopyShuffleZeroPadASM_ReInput ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void UnShuffleTruncateASM_ReOutput(ReIm *Fx, double *ReX, double *ImX, int N);}
;
; Copy data in shuffled order to F space and zero pad, ready for DIT FFT
;
; rcx : Fx
; rdx : ReX
; r8 : ImX
; r9 : N
;
PUBLIC UnShuffleTruncateASM_ReOutput
UnShuffleTruncateASM_ReOutput PROC

	push rsi

;<SHUFFLE>
; rcx : Fx (input)
; rdx : ReX (output)
; r8  : ImX (output)
; rsi : fftShufInv + N - 1 + p
; r9  : N/2 (loop counter)
;

	mov rsi, r9
	dec rsi
	shl rsi, 2								;(N - 1) * 4
	add rsi, fftShufInv						;fftShufInv + N - 1

	shr r9, 1								;N/2 - number of loop iterations

UnShuffleTruncateASM_ReOutput_Loop:

	xor rax, rax
	mov eax, [rsi]							;index = fftShuffleInv[N - 1 + p]
	shl rax, 3								;index*8

	movupd xmm1, [rcx]						;Fx[p] : Im, Re

	mov r10, rdx
	add r10, rax							;ReX + index
	movlpd QWORD PTR [r10], xmm1			;ReX[index] = Fx[p].Re

	mov r10, r8
	add r10, rax							;ImX + index
	movhpd QWORD PTR [r10], xmm1		    ;ImX[index] = Fx[p].Im

	add rsi, 8								;next fftShufInv + N - 1 + 2*p
	add rcx, 32								;next Fx + 2*p

	dec r9
	jnz UnShuffleTruncateASM_ReOutput_Loop

;</SHUFFLE>

	pop rsi

	ret

UnShuffleTruncateASM_ReOutput ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void CopyRealShuffleZeroPadASM(double *Rex, ReIm *FX, int N);}
;
; Copy data in shuffled order to F space and zero pad, ready for DIT FFT
;
; rcx : Rex
; rdx : FX
; r8 : N
;
PUBLIC CopyRealShuffleZeroPadASM
CopyRealShuffleZeroPadASM PROC

	push rsi

;<SHUFFLE>
; rcx : Rex (input)
; rdx : FX (output)
; rsi : fftShufInv + N - 1 + p
; r8  : N/2 (loop counter)
;

	mov rsi, r8
	dec rsi
	shl rsi, 2								;(N - 1) * 4
	add rsi, fftShufInv						;fftShufInv + N - 1

	shr r8, 1								;N/2 - number of loop iterations

	xorpd xmm0, xmm0

CopyRealShuffleZeroPadASM_Loop:

	xor rax, rax
	mov eax, [rsi]							;index = fftShuffleInv[N - 1 + p]

	shl rax, 4								;index * 16
	add rax, rdx							;FX + index

	movupd xmm1, [rcx]						;Rex[2*p + 1], Rex[2*p]

	movupd [rax], xmm1						;FX[index]
	movupd [rax + 16], xmm0					;FX[index + 1] = ReIm(0, 0) : zero-padding

	add rsi, 4								;next fftShufInv + N - 1 + p
	add rcx, 16								;next Rex + 2*p

	dec r8
	jnz CopyRealShuffleZeroPadASM_Loop

;</SHUFFLE>

	pop rsi

	ret

CopyRealShuffleZeroPadASM ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void CopyShuffleZeroPadASM(ReIm *Fx, ReIm *FX, int stride, int N);}
;
; Copy data in shuffled order to F space and zero pad, ready for DIT FFT
;
; rcx : Fx
; rdx : FX
; r8 : stride
; r9 : N
;
PUBLIC CopyShuffleZeroPadASM
CopyShuffleZeroPadASM PROC

	push rsi

;<SHUFFLE>
; rcx : Fx (input)
; rdx : FX (output)
; rsi : fftShufInv + N - 1 + p
; r8  : stride*16 in Fx
; r9  : N/2 (loop counter)
;

	mov rsi, r9
	dec rsi
	shl rsi, 2								;(N - 1) * 4
	add rsi, fftShufInv						;fftShufInv + N - 1

	shr r9, 1								;N/2 - number of loop iterations

	shl r8, 4								;stride * 16

	xorpd xmm0, xmm0

CopyShuffleZeroPadASM_Loop:

	xor rax, rax
	mov eax, [rsi]							;index = fftShuffleInv[N - 1 + p]

	shl rax, 4								;index * 16
	add rax, rdx							;FX + index

	movupd xmm1, [rcx]						;Fx[p * stride]

	movupd [rax], xmm1						;FX[index]
	movupd [rax + 16], xmm0					;FX[index + 1] = ReIm(0, 0) : zero-padding

	add rsi, 4								;next fftShufInv + N - 1 + p
	add rcx, r8								;next Fx + p

	dec r9
	jnz CopyShuffleZeroPadASM_Loop

;</SHUFFLE>

	pop rsi

	ret

CopyShuffleZeroPadASM ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void UnShuffleTruncateASM(ReIm *Fx, ReIm *FX, int stride, int N);}
;
; Copy data in shuffled order to F space and zero pad, ready for DIT FFT
;
; rcx : Fx
; rdx : FX
; r8 : stride
; r9 : N
;
PUBLIC UnShuffleTruncateASM
UnShuffleTruncateASM PROC

	push rsi

;<SHUFFLE>
; rcx : Fx (input)
; rdx : FX (output)
; rsi : fftShufInv + N - 1 + p
; r8  : stride*16 in FX
; r9  : N/2 (loop counter)
;

	mov rsi, r9
	dec rsi
	shl rsi, 2								;(N - 1) * 4
	add rsi, fftShufInv						;fftShufInv + N - 1

	shr r9, 1								;N/2 - number of loop iterations

	shl r8, 4								;stride * 16

UnShuffleTruncateASM_Loop:

	xor rax, rax
	mov eax, [rsi]							;index = fftShuffleInv[N - 1 + p]
	imul rax, r8							;index * stride * 16
	add rax, rdx							;FX + index*stride

	movupd xmm1, [rcx]						;Fx[p]

	movupd [rax], xmm1						;FX[index*stride]

	add rsi, 8								;next fftShufInv + N - 1 + 2*p
	add rcx, 32								;next Fx + 2*p

	dec r9
	jnz UnShuffleTruncateASM_Loop

;</SHUFFLE>

	pop rsi

	ret

UnShuffleTruncateASM ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void UnShuffleRealTruncateASM(ReIm *Fx, double *ReX, int N);}
;
; Copy data in shuffled order to F space and zero pad, ready for DIT FFT
;
; rcx : Fx
; rdx : ReX
; r8 : N
;
PUBLIC UnShuffleRealTruncateASM
UnShuffleRealTruncateASM PROC

	push rsi

;<SHUFFLE>
; rcx : Fx (input)
; rdx : ReX (output)
; rsi : fftShufInv + N - 1 + p
; r8  : N/2 (loop counter)
;

	mov rsi, r8
	dec rsi
	shl rsi, 2								;(N - 1) * 4
	add rsi, fftShufInv						;fftShufInv + N - 1

	shr r8, 1								;N/2 - number of loop iterations

UnShuffleRealTruncateASM_Loop:

	xor rax, rax
	mov eax, [rsi]							;index = fftShuffleInv[N - 1 + p]
	shl rax, 4								;2 * (index * 8)
	add rax, rdx							;FX + index

	movupd xmm1, [rcx]						;Fx[p].Im, Fx[p].Re
	
	movupd [rax], xmm1						;ReX[2*index] = Fx[p].Re, ReX[2*index + 1] = Fx[p].Im

	add rsi, 8								;next fftShufInv + N - 1 + 2*p
	add rcx, 32								;next Fx + 2*p

	dec r8
	jnz UnShuffleRealTruncateASM_Loop

;</SHUFFLE>

	pop rsi

	ret

UnShuffleRealTruncateASM ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void RealfromComplexFFTASM(ReIm *Fx, ReIm *FX, int N);}
;
; rcx : Fx (input)
; rdx : FX (output)
; r8  : N
;

PUBLIC RealfromComplexFFTASM

RealfromComplexFFTASM PROC

	push rsi

	mov rax, MAXFFT
	mov rsi, r8

	push rdx
	xor rdx, rdx
	div esi									;MAXFFT / N
	pop rdx

	shl rax, 3								;(MAXFFT / 2N) * 16

	mov r9, rcx								;Fx
	
	mov rcx, rax							;(MAXFFT / 2N) * 16
	mov rax, sincos							;sincos

	mov rsi, r8
	shl rsi, 4								;N * 16
	mov r10, r9
	add r10, rsi							;Fx + N

	dec r8									;N - 1 : loop counter
	
;<LOOP>
;
; rax : sincos + (MAXFFT / 2N) * n
; rcx : (MAXFFT / 2N) * 16
; rdx : FX + n
; r8 : loop counter
; r9 : Fx + n
; r10 : Fx + N - n
; 

	movupd xmm1, [r9]						;Fx[0].Im, Fx[0].Re
	xorpd xmm2, xmm2
	haddpd xmm1, xmm2						;0, Fx[0].Im + Fx[0].Re

	movupd [rdx], xmm1						;write FX[0]

	movupd xmm7, signmask_

	mov rsi, 2
	cvtsi2sd xmm6, esi						;2 in xmm6
	movddup xmm6, xmm6						;2 in both upper and lower

RealfromComplexFFTASM_Loop:

	add rax, rcx							;sincos + (MAXFFT / 2N) * n
	add r9, 16								;Fx + n
	sub r10, 16								;Fx + N - n
	add rdx, 16								;Fx + n

	movupd xmm1, [r9]						;Yn = Fx[n]
	movupd xmm2, [r10]						;YNmn = Fx[N - n]
	xorpd xmm2, xmm7						;~YNmn

	movupd xmm0, xmm1
	addpd xmm0, xmm2						;Yn + ~YNmn
	subpd xmm1, xmm2						;Yn - ~YNmn

	movupd xmm2, [rax]						;sincos[(MAXFFT/(2*N)) * n]
	movddup xmm3, xmm2						;Re, Re
	unpckhpd xmm2, xmm2						;Im, Im

	movupd xmm4, xmm1
	mulpd xmm4, xmm3						;Imdif * Re, Redif * Re
	mulpd xmm1, xmm2						;Imdif * Im, Redif * Im
	shufpd xmm1, xmm1, 1					;Redif * Im, Imdif * Im
	addsubpd xmm4, xmm1						;sincos[(MAXFFT/(2*N)) * n] * (Yn - ~YNmn)
	
	subpd xmm0, xmm4						;Yn + ~YNmn - sincos[(MAXFFT/(2*N)) * n] * (Yn - ~YNmn)

	divpd xmm0, xmm6						;divide both by 2

	movupd [rdx], xmm0						;write FX[n]

	dec r8
	jnz RealfromComplexFFTASM_Loop

	sub r10, 16								;Fx + N - n
	add rdx, 16								;Fx + n

	movupd xmm1, [r10]						;Fx[0].Im, Fx[0].Re
	xorpd xmm2, xmm2
	hsubpd xmm1, xmm2						;0, Fx[0].Re - Fx[0].Im

	movupd [rdx], xmm1						;write FX[N]

;</LOOP>

	pop rsi

	ret

RealfromComplexFFTASM ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void RealfromComplexIFFTASM(ReIm *Fx, ReIm *FX, int N);}
;
; rcx : Fx (input)
; rdx : FX (output)
; r8  : N
;

PUBLIC RealfromComplexIFFTASM

RealfromComplexIFFTASM PROC

	push rsi

	mov rax, MAXFFT
	mov rsi, r8

	push rdx
	xor rdx, rdx
	div esi									;MAXFFT / N
	pop rdx

	shl rax, 3								;(MAXFFT / 2N) * 16

	mov r9, rcx								;Fx
	
	mov rcx, rax							;(MAXFFT / 2N) * 16
	mov rax, sincos							;sincos

	mov rsi, r8
	shl rsi, 4								;N * 16
	mov r10, r9
	add r10, rsi							;Fx + N
	
;<LOOP>
;
; rax : sincos + (MAXFFT / 2N) * n
; rcx : (MAXFFT / 2N) * 16
; rdx : FX + n
; r8 : loop counter
; r9 : Fx + n
; r10 : Fx + N - n
; 

	movupd xmm7, signmask_

	mov rsi, 2
	cvtsi2sd xmm6, esi						;2 in xmm6
	movddup xmm6, xmm6						;2 in both upper and lower

RealfromComplexIFFTASM_Loop:

	movupd xmm1, [r9]						;Xn = Fx[n]
	movupd xmm2, [r10]						;XNmn = Fx[N - n]
	xorpd xmm2, xmm7						;~XNmn

	movupd xmm0, xmm1
	addpd xmm0, xmm2						;Xn + ~XNmn
	subpd xmm1, xmm2						;Xn - ~XNmn

	movupd xmm2, [rax]						;sincos[(MAXFFT/(2*N)) * n]
	xorpd xmm2, xmm7						;~sincos[(MAXFFT/(2*N)) * n]
	movddup xmm3, xmm2						;Re, Re
	unpckhpd xmm2, xmm2						;~Im, ~Im

	movupd xmm4, xmm1
	mulpd xmm4, xmm3						;Imdif * Re, Redif * Re
	mulpd xmm1, xmm2						;Imdif * ~Im, Redif * ~Im
	shufpd xmm1, xmm1, 1					;Redif * ~Im, Imdif * ~Im
	addsubpd xmm4, xmm1						;~sincos[(MAXFFT/(2*N)) * n] * (Xn - ~XNmn)
	
	subpd xmm0, xmm4						;Xn + ~XNmn - ~sincos[(MAXFFT/(2*N)) * n] * (Xn - ~XNmn)

	divpd xmm0, xmm6						;divide both by 2

	movupd [rdx], xmm0						;write FX[n]

	add rax, rcx							;sincos + (MAXFFT / 2N) * n
	add r9, 16								;Fx + n
	sub r10, 16								;Fx + N - n
	add rdx, 16								;Fx + n

	dec r8
	jnz RealfromComplexIFFTASM_Loop

;</LOOP>

	pop rsi

	ret

RealfromComplexIFFTASM ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ConvASM_2D(int startRow);} This needs to be redone
;
;
; rcx : (m) startRow
;

PUBLIC ConvASM_2D

ConvASM_2D PROC

	push rbx
	push rsi
	push rdi
	push r12
	push r13
	push r14
	push r15

;<KERNEL MULTIPLICATION>
; 
; rax : F2 + n + m*M             (input)
; rbx : F2 + (M - n) + (N - m)*M (input)
; rcx : F + n + m*M              (output)
; rdx : F + (M - n) + (N - m)*M  (output)
; r10 : loop counter
; r11 : Kx + n + m*M
; r12 : Kx + (M - n) + (N - m)*M
; r13 : Ky + n + m*M
; r14 : Ky + (M - n) + (N - m)*M
;

	mov r10, M
	dec r10									;M - 1 points per row to do

	mov rsi, rcx
	imul rsi, M
	shl rsi, 4								;startRow * M * 16
	
	mov rdi, N
	sub rdi, rcx
	imul rdi, M
	shl rdi, 4								;(N - startRow) * M * 16

	mov rcx, F
	add rcx, rsi							;F + m*M
	mov rax, F2
	add rax, rsi							;F2 + m*M
	mov r11, Kx
	add r11, rsi							;Kx + m*M
	mov r13, Ky
	add r13, rsi							;Ky + m*M

	mov rbx, M
	shl rbx, 4								;M * 16
	mov rdx, rbx
	mov r12, rbx
	mov r14, rbx
	
	add rdx, F								
	add rdx, rdi							;F + M + (N - startRow)*M
	add rbx, F2
	add rbx, rdi							;F2 + M + (N - startRow)*M
	add r12, Kx
	add r12, rdi							;Kx + M + (N - startRow)*M
	add r14, Ky
	add r14, rdi							;Ky + M + (N - startRow)*M

	movupd xmm8, _signmask_					;upper and lower double negation
	movupd xmm7, signmask_					;upper double negation : complex conjugation

ConvASM_2D_RowLoop:

	add rax, 16								;next n : F2 + n + m*M
	add r11, 16								;next n : Kx + n + m*M
	add r13, 16								;next n : Ky + n + m*M
	add rcx, 16								;next n in output : F + m + n*N

	sub rbx, 16								;next n : F2 + (M - n) + (N - m)*M
	sub r12, 16								;next n : Kx + (M - n) + (N - m)*M
	sub r14, 16								;next n : Ky + (M - n) + (N - m)*M
	sub rdx, 16								;next n in output : F + (N - m) + (M - n)*N

	movupd xmm1, [rax]						;F2[n + m*M]
	movupd xmm2, [rbx]						;F2[(M - n) + (N - m)*M]
	xorpd xmm2, xmm7						;~F2[(M - n) + (N - m)*M]

	movupd xmm3, xmm1
	addpd xmm1, xmm2						;sum
	subpd xmm3, xmm2						;dif

	movddup xmm2, xmm1						;ReS, ReS
	unpckhpd xmm1, xmm1						;ImS, ImS
	
	movddup xmm4, xmm3						;ReD, ReD
	unpckhpd xmm3, xmm3						;ImD, ImD

	movupd xmm5, [r11]						;Kx[n + m*M]
	movupd xmm6, [r13]						;Ky[n + m*M]

	movupd xmm0, xmm5
	mulpd xmm5, xmm2						;ImKx*ReS, ReKx*ReS
	mulpd xmm0, xmm1						;ImKx*ImS, ReKx*ImS
	shufpd xmm0, xmm0, 1					;ReKx*ImS, ImKx*ImS
	addsubpd xmm5, xmm0						;Kx * sum

	movupd xmm0, xmm6
	mulpd xmm6, xmm4						;ImKy*ReD, ReKy*ReD
	mulpd xmm0, xmm3						;ImKy*ImD, ReKy*ImD
	shufpd xmm0, xmm0, 1					;ReKy*ImD, ImKy*ImD
	addsubpd xmm6, xmm0						;Ky * dif

	subpd xmm5, xmm6						;Kx * sum - Ky * dif
	movupd [rcx], xmm5						;write F[m + n*N]

	movupd xmm5, [r12]						;Kx[(M - n) + (N - m)*M]
	movupd xmm6, [r14]						;Ky[(M - n) + (N - m)*M]

	xorpd xmm1, xmm8						;conjugate sum (change sign on imaginary parts)
	xorpd xmm3, xmm8						;conjugate dif (change sign on imaginary parts)

	movupd xmm0, xmm5
	mulpd xmm5, xmm2						;ImKx*ReS, ReKx*ReS
	mulpd xmm0, xmm1						;ImKx*ImS, ReKx*ImS
	shufpd xmm0, xmm0, 1					;ReKx*ImS, ImKx*ImS
	addsubpd xmm5, xmm0						;Kx * ~sum

	movupd xmm0, xmm6
	mulpd xmm6, xmm4						;ImKy*ReD, ReKy*ReD
	mulpd xmm0, xmm3						;ImKy*ImD, ReKy*ImD
	shufpd xmm0, xmm0, 1					;ReKy*ImD, ImKy*ImD
	addsubpd xmm6, xmm0						;Ky * ~dif

	addpd xmm5, xmm6						;Kx * ~sum + Ky * ~dif
	movupd [rdx], xmm5						;write F[(N - m) + (M - n)*N]

	dec r10
	jnz ConvASM_2D_RowLoop

;</KERNEL MULTIPLICATION>

	pop r15
	pop r14
	pop r13
	pop r12
	pop rdi
	pop rsi
	pop rbx

	ret

ConvASM_2D ENDP


;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ConvASM_3D(int startPlane);}
;
;
; rcx : (m) startPlane
;

PUBLIC ConvASM_3D

ConvASM_3D PROC

	push rbx
	push rsi
	push rdi
	push r12
	push r13
	push r14
	push r15

;<KERNEL MULTIPLICATION>
;
; rbx : loop counter 
;
; rcx : Kx + m*(N/2 + 1)*K + p
; rdx : Ky + m*(N/2 + 1)*K + p
; rsi : Kz + m*(N/2 + 1)*K + p
; rdi : Kxy + m*(N/2 + 1)*K + p
; r8  : Kxz + m*(N/2 + 1)*K + p
; r9  : Kyz + m*(N/2 + 1)*K + p
;
; r10 : F4 + m*(N/2 + 1)*K + p (input)
; r11 : F5 + m*(N/2 + 1)*K + p (input)
; r12 : F6 + m*(N/2 + 1)*K + p (input)
; r13 : F + m*(N/2 + 1)*K + p  (output)
; r14 : F2 + m*(N/2 + 1)*K + p (output)
; r15 : F3 + m*(N/2 + 1)*K + p (output)

	mov rax, N
	shr rax, 1
	inc rax
	imul rax, K								;(N/2 + 1)*K
	mov rbx, rax							;number of loop iterations
	imul rax, rcx							;(N/2 + 1)*K*m
	shl rax, 4								;(N/2 + 1)*K*m * 16

	mov rcx, Kx
	mov rdx, Ky
	mov rsi, Kz
	mov rdi, Kxy
	mov r8, Kxz
	mov r9, Kyz

	mov r10, F4
	mov r11, F5
	mov r12, F6
	mov r13, F
	mov r14, F2
	mov r15, F3

	add rcx, rax							;Kx + m*(N/2 + 1)*K
	add rdx, rax							;Ky + m*(N/2 + 1)*K
	add rsi, rax							;Kz + m*(N/2 + 1)*K
	add rdi, rax							;Kxy + m*(N/2 + 1)*K
	add r8, rax								;Kxz + m*(N/2 + 1)*K
	add r9, rax								;Kyz + m*(N/2 + 1)*K

	add r10, rax							;F4 + m*(N/2 + 1)*K
	add r11, rax							;F5 + m*(N/2 + 1)*K
	add r12, rax							;F6 + m*(N/2 + 1)*K
	add r13, rax							;F + m*(N/2 + 1)*K
	add r14, rax							;F2 + m*(N/2 + 1)*K
	add r15, rax							;F3 + m*(N/2 + 1)*K

ConvASM_3D_Loop:

		;for(int n = 0; n <= N/2; n++) {
		;	for(int k = 0; k < K; k++) {
		
		;		ReIm FMx, FMy, FMz;

		;		FMx = F4[k + n*K + m*(N/2 + 1)*K];
		;		FMy = F5[k + n*K + m*(N/2 + 1)*K];
		;		FMz = F6[k + n*K + m*(N/2 + 1)*K];

		;		F[k + n*K + m*(N/2 + 1)*K] = Kx[k + n*K + m*(N/2 + 1)*K] * FMx + Kxy[k + n*K + m*(N/2 + 1)*K] * FMy + Kxz[k + n*K + m*(N/2 + 1)*K] * FMz;
				
		;		F2[k + n*K + m*(N/2 + 1)*K] = Kxy[k + n*K + m*(N/2 + 1)*K] * FMx + Ky[k + n*K + m*(N/2 + 1)*K] * FMy + Kyz[k + n*K + m*(N/2 + 1)*K] * FMz;

		;		F3[k + n*K + m*(N/2 + 1)*K] = Kxz[k + n*K + m*(N/2 + 1)*K] * FMx + Kyz[k + n*K + m*(N/2 + 1)*K] * FMy + Kz[k + n*K + m*(N/2 + 1)*K] * FMz;
		;	}
 		;}

	;xmm1 : F
	;xmm2 : F2
	;xmm3 : F3

	;xmm4  : Fmx : Im, Im
	;xmm5  : Fmy : Im, Im
	;xmm6  : Fmz : Im, Im
	;xmm13 : Fmx : Re, Re
	;xmm14 : Fmy : Re, Re
	;xmm15 : Fmz : Re, Re

	movupd xmm4, [r10]						;Fmx
	movupd xmm5, [r11]						;Fmy
	movupd xmm6, [r12]						;Fmz

	movddup xmm13, xmm4						;ReFmx, ReFmx
	unpckhpd xmm4, xmm4						;ImFmx, ImFmx
	
	movddup xmm14, xmm5						;ReFmy, ReFmy
	unpckhpd xmm5, xmm5						;ImFmy, ImFmy
	
	movddup xmm15, xmm6						;ReFmz, ReFmz
	unpckhpd xmm6, xmm6						;ImFmz, ImFmz

	;F = Kx * Fmx
	movupd xmm7, [rcx]						;Kx

	movupd xmm1, xmm7						
	mulpd xmm1, xmm13						;ImKx * ReFmx, ReKx * ReFmx
	mulpd xmm7, xmm4						;ImKx * ImFmx, ReKx * ImFmx
	shufpd xmm7, xmm7, 1					;ReKx * ImFmx, ImKx * ImFmx
	addsubpd xmm1, xmm7						;F = Kx * Fmx

	;F2 = Ky * Fmy
	movupd xmm8, [rdx]						;Ky
	
	movupd xmm2, xmm8						
	mulpd xmm2, xmm14						;ImKy * ReFmy, ReKy * ReFmy
	mulpd xmm8, xmm5						;ImKy * ImFmy, ReKy * ImFmy
	shufpd xmm8, xmm8, 1					;ReKy * ImFmy, ImKy * ImFmy
	addsubpd xmm2, xmm8						;F2 = Ky * Fmy

	;F3 = Kz * Fmz
	movupd xmm9, [rsi]						;Kz

	movupd xmm3, xmm9						
	mulpd xmm3, xmm15						;ImKz * ReFmz, ReKz * ReFmz
	mulpd xmm9, xmm6						;ImKz * ImFmz, ReKz * ImFmz
	shufpd xmm9, xmm9, 1					;ReKz * ImFmz, ImKz * ImFmz
	addsubpd xmm3, xmm9						;F3 = Kz * Fmz
	

	;Kxy contributions
	movupd xmm10, [rdi]						;Kxy
	
	movupd xmm7, xmm10						
	mulpd xmm7, xmm14						;ImKxy * ReFmy, ReKxy * ReFmy
	movupd xmm8, xmm10
	mulpd xmm8, xmm5						;ImKxy * ImFmy, ReKxy * ImFmy
	shufpd xmm8, xmm8, 1					;ReKxy * ImFmy, ImKxy * ImFmy
	addsubpd xmm7, xmm8						;Kxy * Fmy
	addpd xmm1, xmm7						;contribution to F

	movupd xmm7, xmm10						
	mulpd xmm7, xmm13						;ImKxy * ReFmx, ReKxy * ReFmx
	mulpd xmm10, xmm4						;ImKxy * ImFmx, ReKxy * ImFmx
	shufpd xmm10, xmm10, 1					;ReKxy * ImFmx, ImKxy * ImFmx
	addsubpd xmm7, xmm10					;Kxy * Fmx
	addpd xmm2, xmm7						;contribution to F2


	;Kxz contributions
	movupd xmm11, [r8]						;Kxz
	
	movupd xmm7, xmm11						
	mulpd xmm7, xmm15						;ImKxz * ReFmz, ReKxz * ReFmz
	movupd xmm8, xmm11
	mulpd xmm8, xmm6						;ImKxz * ImFmz, ReKxz * ImFmz
	shufpd xmm8, xmm8, 1					;ReKxz * ImFmz, ImKxz * ImFmz
	addsubpd xmm7, xmm8						;Kxz * Fmz
	addpd xmm1, xmm7						;contribution to F

	movupd xmm7, xmm11						
	mulpd xmm7, xmm13						;ImKxz * ReFmx, ReKxz * ReFmx
	mulpd xmm11, xmm4						;ImKxz * ImFmx, ReKxz * ImFmx
	shufpd xmm11, xmm11, 1					;ReKxz * ImFmx, ImKxz * ImFmx
	addsubpd xmm7, xmm11					;Kxz * Fmx
	addpd xmm3, xmm7						;contribution to F3

	
	;Kyz contributions
	movupd xmm12, [r9]						;Kyz
	
	movupd xmm7, xmm12						
	mulpd xmm7, xmm15						;ImKyz * ReFmz, ReKyz * ReFmz
	movupd xmm8, xmm12
	mulpd xmm8, xmm6						;ImKyz * ImFmz, ReKyz * ImFmz
	shufpd xmm8, xmm8, 1					;ReKyz * ImFmz, ImKyz * ImFmz
	addsubpd xmm7, xmm8						;Kxz * Fmz
	addpd xmm2, xmm7						;contribution to F2

	movupd xmm7, xmm12						
	mulpd xmm7, xmm14						;ImKyz * ReFmy, ReKyz * ReFmy
	mulpd xmm12, xmm5						;ImKyz * ImFmy, ReKyz * ImFmy
	shufpd xmm12, xmm12, 1					;ReKyz * ImFmy, ImKyz * ImFmy
	addsubpd xmm7, xmm12					;Kyz * Fmy
	addpd xmm3, xmm7						;contribution to F3

	movupd [r13], xmm1						;write F
	movupd [r14], xmm2						;write F2
	movupd [r15], xmm3						;write F3

	add rcx, 16
	add rdx, 16
	add rsi, 16
	add rdi, 16
	add r8, 16
	add r9, 16
	add r10, 16
	add r11, 16
	add r12, 16
	add r13, 16
	add r14, 16
	add r15, 16

	dec rbx
	jnz ConvASM_3D_Loop

;</KERNEL MULTIPLICATION>

	pop r15
	pop r14
	pop r13
	pop r12
	pop rdi
	pop rsi
	pop rbx

	ret

ConvASM_3D ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void FFTASM_Radix4_DIT(ReIm *FX, int ldn, int N);}
;
; Do Radix-4 FFT using DIT on Input Fx, Output placed in FX
; Fx not affected, input expected to be contiguous in memory, arranged as (Re, Im) doubles, output likewise
; N FFT length ldn = log2(N)
; Fastcall convention : 
; rcx : FX
; rdx : ldn
; r8 : N

PUBLIC FFTASM_Radix4_DIT
FFTASM_Radix4_DIT PROC

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
; r10 : FX + i0 * 16
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
	
	movupd xmm0, [r10]						;FX[i0]
	movupd xmm1, [r10 + 16]					;FX[i0 + 1]
	movupd xmm2, xmm0
	addpd  xmm0, xmm1						;FX[i0] + FX[(i0 + 1)]
	subpd  xmm2, xmm1						;FX[i0] - FX[(i0 + 1)]
	movupd [r10], xmm0
	movupd [r10 + 16], xmm2	

	add r10, 32								;FX + i0 + 2

	dec rbx
	jnz FFTASM_Radix4_DIT_Radix2Loop

;</RADIX 2 STEP?>

	shr r15, 1								;MAXFFT / N2 since N2 needs to be 2 (set below). If we jump past this statement then N2 is 1 and this already holds the right value.

FFTASM_Radix4_DIT_SkipRadix2:

;<RADIX 4 PART>
; rax : N4 * 16
; rbx : N2 * 16
; rcx : FX
; rdx : outer loop counter
; rsi : j * 16
; rdi : r * 16
; r8  : N * 16
; r10 : i0 * 16 - used in inner loop only
; r11 : i1 * 16 - used in inner loop only
; r12 : i2 * 16 - used in inner loop only
; r13 : i3 * 16 - used in inner loop only
; r14 : cossin
; r15  : MAXFFT / N2

	mov r14, cossin

	shl r8, 4								;N * 16
	inc r9									;1 << p, where p = ldn&1 : if we jumped here because p was zero then 1 << p = 1. Otherwise p was 1 and now 1 << p = 2
	mov rbx, r9								;N2
	shl rbx, 4								;N2 * 16
	
	dec rdx
	sub rdx, r9
	shr rdx, 1								
	inc rdx									;(ldn - 1 - 1<<p) / 2 + 1: number of outer loop iterations

	movupd	xmm8, signmask					;will need this for multiplication by -1 on the lower double
	movupd	xmm9, signmask_					;will need this for multiplication by -1 on the upper double

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

	movupd xmm4, [r10]						;FX[i0]
	movupd xmm5, [r11]						;FX[i1]	
	movupd xmm6, [r12]						;FX[i2]	
	movupd xmm7, [r13]						;FX[i3]	

	movupd xmm0, xmm4
	addpd xmm4, xmm5						;t0 = FX[i0] + FX[i1]
	subpd xmm0, xmm5						;t1 = FX[i0] - FX[i1]
	
	movupd xmm5, xmm7
	subpd xmm5, xmm6						;t3 = FX[i3] - FX[i2] * exp1
	shufpd xmm5, xmm5, 1					;Re3, Im3
	xorpd xmm5, xmm8						;Re3, -Im3 : t3 = i*t3 now
	addpd xmm6, xmm7						;t2 = FX[i2] * exp1 + FX[i3]

	movupd xmm7, xmm4
	addpd xmm4, xmm6						;t0 + t2
	subpd xmm7, xmm6						;t0 - t2
	
	movupd xmm6, xmm0
	addpd xmm0, xmm5						;t1 + t3
	subpd xmm6, xmm5						;t1 - t3

	movupd [r10], xmm4
	movupd [r11], xmm0
	movupd [r12], xmm7
	movupd [r13], xmm6	

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl FFTASM_Radix4_DIT_SpecialLoop

	mov rsi, 16								;j = 1 start index
	cmp rsi, rax
	jge FFTASM_Radix4_DIT_NextOuterLoop

FFTASM_Radix4_DIT_MiddleLoop:

	xor rdi, rdi							;r = 0 start index
	
	mov r13, r15
	imul r13, rsi							;(MAXFFT / N2) * j * 16
	mov r12, r14
	add r12, r13							;cossin + (MAXFFT / N2) * j
	
	movupd xmm1, [r12]						;exp  = cossin[(MAXFFT / N2) * j]
	xorpd xmm1, xmm9						;conjugate exp : change sign of sin, the upper double
	movddup xmm11, xmm1						;cos1, cos1
	unpckhpd xmm1, xmm1						;sin1, sin1

	add r12, r13							;cossin + (MAXFFT / N2) * 2*j
	movupd xmm2, [r12]						;exp2 = cossin[(MAXFFT / N2) * 2*j]
	xorpd xmm2, xmm9						;conjugate exp : change sign of sin, the upper double
	movddup xmm12, xmm2						;cos2, cos2
	unpckhpd xmm2, xmm2						;sin2, sin2

	add r12, r13							;cossin + (MAXFFT / N2) * 3*j
	movupd xmm3, [r12]						;exp3 = cossin[(MAXFFT / N2) * 3*j]
	xorpd xmm3, xmm9						;conjugate exp : change sign of sin, the upper double
	movddup xmm13, xmm3						;cos3, cos3
	unpckhpd xmm3, xmm3						;sin3, sin3

	mov r10, rsi
	add r10, rcx							;FX + i0, where i0 = j at start
	mov r11, r10
	add r11, rax							;FX + i1, where i1 = i0 + N4
	mov r12, r11
	add r12, rax							;FX + i2, where i2 = i1 + N4
	mov r13, r12
	add r13, rax							;FX + i3, where i3 = i2 + N4

FFTASM_Radix4_DIT_InnerLoop:

	movupd xmm0, [r11]						;FX[i1]	
	movupd xmm5, xmm0
	mulpd xmm0, xmm2						;Im1*sin2, Re1*sin2
	mulpd xmm5, xmm12						;Im1*cos2, Re1*cos2
	shufpd xmm0, xmm0, 1					;Re1*sin2, Im1*sin2
	addsubpd xmm5, xmm0						;FX[i1] * exp2

	movupd xmm0, [r12]						;FX[i2]	
	movupd xmm6, xmm0
	mulpd xmm0, xmm1						;Im2*sin1, Re2*sin1
	mulpd xmm6, xmm11						;Im2*cos1, Re2*cos1
	shufpd xmm0, xmm0, 1					;Re2*sin1, Im2*sin1
	addsubpd xmm6, xmm0						;FX[i2] * exp1

	movupd xmm0, [r13]						;FX[i3]	
	movupd xmm7, xmm0
	mulpd xmm0, xmm3						;Im3*sin3, Re3*sin3
	mulpd xmm7, xmm13						;Im3*cos3, Re3*cos3
	shufpd xmm0, xmm0, 1					;Re3*sin3, Im3*sin3
	addsubpd xmm7, xmm0						;FX[i3] * exp3

	movupd xmm4, [r10]						;FX[i0]

	movupd xmm0, xmm4
	addpd xmm4, xmm5						;t0 = FX[i0] + FX[i1] * exp2
	subpd xmm0, xmm5						;t1 = FX[i0] - FX[i1] * exp2
	
	movupd xmm5, xmm7
	subpd xmm5, xmm6						;t3 = FX[i3] * exp3 - FX[i2] * exp1
	shufpd xmm5, xmm5, 1					;Re3, Im3
	xorpd xmm5, xmm8						;Re3, -Im3 : t3 = i*t3 now
	addpd xmm6, xmm7						;t2 = FX[i2] * exp1 + FX[i3] * exp3

	movupd xmm7, xmm4
	addpd xmm4, xmm6						;t0 + t2
	subpd xmm7, xmm6						;t0 - t2
	
	movupd xmm6, xmm0
	addpd xmm0, xmm5						;t1 + t3
	subpd xmm6, xmm5						;t1 - t3

	movupd [r10], xmm4
	movupd [r11], xmm0
	movupd [r12], xmm7
	movupd [r13], xmm6	

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl FFTASM_Radix4_DIT_InnerLoop

	add rsi, 16
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

FFTASM_Radix4_DIT ENDP


;;;;;;;;;;;;;;;;;;;;;;extern "C" {void IFFTASM_Radix4_DIF(ReIm *FX, int ldn, int N);}
;
; Do Radix-4 IFFT using DIF on Input Fx, Output placed in FX
; input expected to be contiguous in memory, arranged as (Re, Im) doubles, output likewise
; N FFT length ldn = log2(N)
; Fastcall convention : 
; rcx : FX
; rdx : ldn
; r8 : N

PUBLIC IFFTASM_Radix4_DIF
IFFTASM_Radix4_DIF PROC

	push rbx
	push rsi
	push rdi
	push r12
	push r13
	push r14
	push r15
		
;<RADIX 4 PART>
; rax : N4 * 16
; rbx : N2 * 16
; rcx : FX
; rdx : outer loop counter
; rsi : j * 16
; rdi : r * 16
; r8  : N * 16
; r9  : ldn
; r10 : i0 * 16 - used in inner loop only
; r11 : i1 * 16 - used in inner loop only
; r12 : i2 * 16 - used in inner loop only
; r13 : i3 * 16 - used in inner loop only
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

	shl r8, 4								;N * 16
	mov rbx, r8								
	shl rbx, 2								;N2 = (4 * N) * 16
	mov rax, r8								;set N4 = N at start

	mov rdx, r9								;get ldn back
	dec rdx
	dec rdx
	shr rdx, 1														
	inc rdx									;(ldn - 2) / 2 + 1: number of outer loop iterations

	movupd	xmm8, signmask					;will need this for multiplication by -1 on the lower double

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

	movupd xmm4, [r10]						;FX[i0]
	movupd xmm5, [r11]						;FX[i1]	
	movupd xmm6, [r12]						;FX[i2]	
	movupd xmm7, [r13]						;FX[i3]	

	movupd xmm0, xmm4
	addpd xmm0, xmm6						;t0 = FX[i0] + FX[i2]
	subpd xmm4, xmm6						;t1 = FX[i0] - FX[i2]

	movupd xmm6, xmm5
	addpd xmm6, xmm7						;t2 = FX[i1] + FX[i3]
	subpd xmm7, xmm5						;t3 = FX[i3] - FX[i1]
	shufpd xmm7, xmm7, 1					;Re3, Im3
	xorpd xmm7, xmm8						;Re3, -Im3 : t3 = i*t3 now

	movupd xmm5, xmm0
	addpd xmm5, xmm6						;t0 + t2
	subpd xmm0, xmm6						;t0 - t2	

	movupd xmm6, xmm4
	addpd xmm6, xmm7						;t1 + t3
	subpd xmm4, xmm7						;t1 - t3	

	movupd [r10], xmm5						;t0 + t2
	movupd [r11], xmm0						;t0 - t2
	movupd [r12], xmm4						;t1 - t3
	movupd [r13], xmm6						;t1 + t3

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl IFFTASM_Radix4_DIF_SpecialLoop	

	mov rsi, 16								;j = 1 start index
	cmp rsi, rax
	jge IFFTASM_Radix4_DIF_NextOuterLoop

IFFTASM_Radix4_DIF_MiddleLoop:

	xor rdi, rdi							;r = 0 start index
	
	mov r13, r15
	imul r13, rsi							;(MAXFFT / N2) * j * 16
	mov r12, r14
	add r12, r13							;cossin + (MAXFFT / N2) * j
	
	movupd xmm1, [r12]						;exp  = cossin[(MAXFFT / N2) * j]
	movddup xmm11, xmm1						;cos1, cos1
	unpckhpd xmm1, xmm1						;sin1, sin1

	add r12, r13							;cossin + (MAXFFT / N2) * 2*j
	movupd xmm2, [r12]						;exp2 = cossin[(MAXFFT / N2) * 2*j]
	movddup xmm12, xmm2						;cos2, cos2
	unpckhpd xmm2, xmm2						;sin2, sin2

	add r12, r13							;cossin + (MAXFFT / N2) * 3*j
	movupd xmm3, [r12]						;exp3 = cossin[(MAXFFT / N2) * 3*j]
	movddup xmm13, xmm3						;cos3, cos3
	unpckhpd xmm3, xmm3						;sin3, sin3

	mov r10, rsi
	add r10, rcx							;FX + i0, where i0 = j at start
	mov r11, r10
	add r11, rax							;FX + i1, where i1 = i0 + N4
	mov r12, r11
	add r12, rax							;FX + i2, where i2 = i1 + N4
	mov r13, r12
	add r13, rax							;FX + i3, where i3 = i2 + N4

IFFTASM_Radix4_DIF_InnerLoop:

	movupd xmm4, [r10]						;FX[i0]
	movupd xmm5, [r11]						;FX[i1]	
	movupd xmm6, [r12]						;FX[i2]	
	movupd xmm7, [r13]						;FX[i3]	

	movupd xmm0, xmm4
	addpd xmm0, xmm6						;t0 = FX[i0] + FX[i2]
	subpd xmm4, xmm6						;t1 = FX[i0] - FX[i2]

	movupd xmm6, xmm5
	addpd xmm6, xmm7						;t2 = FX[i1] + FX[i3]
	subpd xmm7, xmm5						;t3 = FX[i3] - FX[i1]
	shufpd xmm7, xmm7, 1					;Re3, Im3
	xorpd xmm7, xmm8						;Re3, -Im3 : t3 = i*t3 now

	movupd xmm5, xmm0
	addpd xmm5, xmm6						;t0 + t2
	subpd xmm0, xmm6						;t0 - t2	

	movupd xmm6, xmm4
	addpd xmm6, xmm7						;t1 + t3
	subpd xmm4, xmm7						;t1 - t3	

	;form (t0-t2)*exp2
	movupd xmm7, xmm0
	mulpd xmm0, xmm2						;Im*sin2, Re*sin2
	mulpd xmm7, xmm12						;Im*cos2, Re*cos2
	shufpd xmm0, xmm0, 1					;Re*sin2, Im*sin2
	addsubpd xmm7, xmm0						;(t0-t2)*exp2		

	;form (t1-t3)*exp
	movupd xmm0, xmm4
	mulpd xmm0, xmm1						;Im*sin, Re*sin
	mulpd xmm4, xmm11						;Im*cos, Re*cos
	shufpd xmm0, xmm0, 1					;Re*sin, Im*sin
	addsubpd xmm4, xmm0						;(t1-t3)*exp

	;form (t1+t3)*exp3
	movupd xmm0, xmm6
	mulpd xmm0, xmm3						;Im*sin3, Re*sin3
	mulpd xmm6, xmm13						;Im*cos3, Re*cos3
	shufpd xmm0, xmm0, 1					;Re*sin3, Im*sin3
	addsubpd xmm6, xmm0						;(t1+t3)*exp3

	movupd [r10], xmm5						;t0 + t2
	movupd [r11], xmm7						;(t0 - t2) * exp2
	movupd [r12], xmm4						;(t1 - t3) * exp
	movupd [r13], xmm6						;(t1 + t3) * exp3	

	add r10, rbx							;next i0 : add N2
	add r11, rbx							;next i1 : add N2
	add r12, rbx							;next i2 : add N2
	add r13, rbx							;next i3 : add N2

	add rdi, rbx
	cmp rdi, r8
	jl IFFTASM_Radix4_DIF_InnerLoop

	add rsi, 16
	cmp rsi, rax
	jl IFFTASM_Radix4_DIF_MiddleLoop

IFFTASM_Radix4_DIF_NextOuterLoop:

	dec rdx
	jnz IFFTASM_Radix4_DIF_OuterLoop

;</RADIX 4 PART>

;<RADIX 2 STEP?>
; rcx : FX
; rbx : loop counter
; r8  : N * 16
; r9  : ldn
; r10 : i0 * 16

	and r9, 1								;ldn&1
	jz IFFTASM_Radix4_DIF_return

	mov r10, rcx							;FX + i0
	mov rbx, r8
	shr rbx, 5								;N / 2 : loop counter

IFFTASM_Radix4_DIF_Radix2Loop:
	
	movupd xmm0, [r10]						;FX[i0]
	movupd xmm1, [r10 + 16]					;FX[i0 + 1]
	movupd xmm2, xmm0
	addpd  xmm0, xmm1						;FX[i0] + FX[(i0 + 1)]
	subpd  xmm2, xmm1						;FX[i0] - FX[(i0 + 1)]
	movupd [r10], xmm0
	movupd [r10 + 16], xmm2	

	add r10, 32								;FX + i0 + 2

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

IFFTASM_Radix4_DIF ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {double SetFields2DGetEnergy(int startRow, double *energyThreads);} This needs to be redone
;
; Used to finish effective field calculation and get energy value : 
; 1. Heffx, Heffy data needs division by NM since it has been IFFT-ed but not yet divided to finish the 2D IFFT
; 2. Heffz = -Mz since this is used for 2D demag field
; 3. energy value must be calculated
;
; This function needs to process data at given row (startRow)
;
; rcx : startRow
; rdx : energyThreads pointer

PUBLIC SetFields2DGetEnergy
SetFields2DGetEnergy PROC

	push rbx
	push r12
	push r13
	push r14

;<CALCULATION LOOP>
; rcx : Heffx + startRow * nx
; rdx : Heffy + startRow * nx
; r8  : Heffz + startRow * nx
; r9  : Mx + startRow * nx
; r10 : My + startRow * nx
; r11 : Mz + startRow * nx 
; r13 : loop counter
; r14 : energy threads pointer
; xmm1 : Heffx
; xmm2 : Heffy
; xmm3 : Heffz
; xmm4 : Mx
; xmm5 : My
; xmm6 : Mz
; xmm7 : N*M for division
; xmm8 : _signmask_

	mov r14, rdx							;save energyThreads

	mov rax, nx
	imul rcx, rax							;startRow * nx
	shl rcx, 3								;startRow * nx * 8								
	mov r12, rcx							;save start offset

	shr rax, 1								;set loop counter to nx / 2 : this is the number of pairs of points this function will process (do 2 at a time since SIMD)
	mov r13, rax							;use r13 as loop counter

	mov rcx, Heffx
	mov rdx, Heffy
	mov r8, Heffz
	mov r9, Mx
	mov r10, My
	mov r11, Mz

	add rcx, r12							;Heffx + startRow * nx
	add rdx, r12							;Heffy + startRow * nx
	add r8, r12								;Heffz + startRow * nx
	add r9, r12								;Mx + startRow * nx
	add r10, r12							;My + startRow * nx
	add r11, r12							;Mz + startRow * nx

	movupd xmm8, _signmask_

	mov rax, N
	imul rax, M								;N*M
	cvtsi2sd xmm7, eax
	movddup xmm7, xmm7						; N*M | N*M for ifft division

	xorpd xmm0, xmm0

SetFields2DGetEnergy_Loop:

	movupd xmm1, [rcx]						;get Heffx 2 points
	movupd xmm2, [rdx]						;get Heffy 2 points

	divpd xmm1, xmm7						;Heffx / NM
	divpd xmm2, xmm7						;Heffy / NM
	
	movupd xmm6, [r11]						;get Mz 2 points
	movupd xmm3, xmm6
	xorpd xmm3, xmm8						;Heffz = -Mz

	movupd xmm4, [r9]						;get Mx 2 points
	movupd xmm5, [r10]						;get My 2 points 

	movupd [rcx], xmm1						;write Heffx
	movupd [rdx], xmm2						;write Heffy
	movupd [r8], xmm3						;write Heffz

	;calculate energy

	mulpd xmm1, xmm4						;Mx*Heffx
	mulpd xmm2, xmm5						;My*Heffy
	mulpd xmm3, xmm6						;Mz*Heffz
	addpd xmm1, xmm2					
	addpd xmm1, xmm3						;Mx*Heffx + My*Heffy + Mz*Heffz
	addpd xmm0, xmm1						;energy += ... Note this is done in pairs of values, will need to add them at the end		
	
	add rcx, 16								;next point
	add rdx, 16								;next point
	add r8, 16								;next point
	add r9, 16								;next point
	add r10, 16								;next point
	add r11, 16								;next point

	dec r13
	jnz SetFields2DGetEnergy_Loop

	xorpd xmm1, xmm1
	haddpd xmm0, xmm1						

	movlpd xmm1, QWORD PTR [r14]
	addpd xmm0, xmm1
	movlpd QWORD PTR [r14], xmm0			;add to energy for this thread - can't get normal double return via xmm0 to work with OMP, don't understand why, so using this method instead.

;</CALCULATION LOOP>

	pop r14
	pop r13
	pop r12
	pop rbx
	
	ret

SetFields2DGetEnergy ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {double SetFields3DGetEnergy(int startRow, int startPlane, double *energyThreads);}
;
; Used to finish effective field calculation and get energy value : 
; 1. Heffx, Heffy, Heffz data needs division by NMK since it has been IFFT-ed but not yet divided to finish the 3D IFFT
; 2. energy value must be calculated
;
; This function needs to process data at given row in given plane (startRow, startPlane)
;
; rcx : startRow
; rdx : startPlane
; r8 : energyThreads pointer

PUBLIC SetFields3DGetEnergy
SetFields3DGetEnergy PROC

	push rbx
	push r12
	push r13
	push r14

;<CALCULATION LOOP>
; rcx : Heffx + startRow * nx + startPlane * nx * ny
; rdx : Heffy + startRow * nx + startPlane * nx * ny
; r8  : Heffz + startRow * nx + startPlane * nx * ny
; r9  : Mx + startRow * nx + startPlane * nx * ny
; r10 : My + startRow * nx + startPlane * nx * ny
; r11 : Mz + startRow * nx + startPlane * nx * ny
; r13 : loop counter
; r14 : energy threads pointer
; xmm1 : Heffx
; xmm2 : Heffy
; xmm3 : Heffz
; xmm4 : Mx
; xmm5 : My
; xmm6 : Mz
; xmm7 : N*M*K for division
; xmm8 : _signmask_

	mov r14, r8								;save energyThreads

	mov rax, nx
	imul rcx, rax							;startRow * nx
	imul rdx, rax							;startPlane * nx
	imul rdx, ny							;startPlane * nx * ny
	mov r12, rcx							
	add r12, rdx							;save start offset
	shl r12, 3								;x8 for double

	shr rax, 1								;set loop counter to nx / 2 : this is the number of pairs of points this function will process (do 2 at a time since SIMD)
	mov r13, rax							;use r13 as loop counter

	mov rcx, Heffx
	mov rdx, Heffy
	mov r8, Heffz
	mov r9, Mx
	mov r10, My
	mov r11, Mz

	add rcx, r12							;Heffx + startRow * nx + startPlane * nx * ny
	add rdx, r12							;Heffy + startRow * nx + startPlane * nx * ny
	add r8, r12								;Heffz + startRow * nx + startPlane * nx * ny
	add r9, r12								;Mx + startRow * nx + startPlane * nx * ny
	add r10, r12							;My + startRow * nx + startPlane * nx * ny
	add r11, r12							;Mz + startRow * nx + startPlane * nx * ny

	mov rax, N
	imul rax, M								;N*M
	imul rax, K								;N*M*K
	shr rax, 1								;(N/2)*M*K
	cvtsi2sd xmm7, eax
	movddup xmm7, xmm7						; (N/2)*M*K | (N/2)*M*K for ifft division

	xorpd xmm0, xmm0

SetFields3DGetEnergy_Loop:

	movupd xmm1, [rcx]						;get Heffx 2 points
	movupd xmm2, [rdx]						;get Heffy 2 points
	movupd xmm3, [r8]						;get Heffz 2 points

	divpd xmm1, xmm7						;Heffx / (N/2)MK
	divpd xmm2, xmm7						;Heffy / (N/2)MK
	divpd xmm3, xmm7						;Heffz / (N/2)MK

	movupd [rcx], xmm1						;write Heffx
	movupd [rdx], xmm2						;write Heffy
	movupd [r8], xmm3						;write Heffz
	
	movupd xmm4, [r9]						;get Mx 2 points
	movupd xmm5, [r10]						;get My 2 points 
	movupd xmm6, [r11]						;get Mz 2 points 

	;calculate energy

	mulpd xmm1, xmm4						;Mx*Heffx
	mulpd xmm2, xmm5						;My*Heffy
	mulpd xmm3, xmm6						;Mz*Heffz
	addpd xmm1, xmm2					
	addpd xmm1, xmm3						;Mx*Heffx + My*Heffy + Mz*Heffz
	addpd xmm0, xmm1						;energy += ... Note this is done in pairs of values, will need to add them at the end		
	
	add rcx, 16								;next point
	add rdx, 16								;next point
	add r8, 16								;next point
	add r9, 16								;next point
	add r10, 16								;next point
	add r11, 16								;next point

	dec r13
	jnz SetFields3DGetEnergy_Loop

	xorpd xmm1, xmm1
	haddpd xmm0, xmm1						

	movlpd xmm1, QWORD PTR [r14]
	addpd xmm0, xmm1
	movlpd QWORD PTR [r14], xmm0			;add to energy for this thread - can't get normal double return via xmm0 to work with OMP, don't understand why, so using this method instead.

;</CALCULATION LOOP>

	pop r14
	pop r13
	pop r12
	pop rbx
	
	ret

SetFields3DGetEnergy ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Call Initialization Functions;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_FFTSpaces(ReIm *F, ReIm *F2, ReIm *F3, ReIm *F4, ReIm *F5, ReIm *F6);}

PUBLIC InitASMDemag_FFTSpaces
InitASMDemag_FFTSpaces PROC

	push rbp
	mov rbp, rsp

	mov F, rcx
	mov F2, rdx
	mov F3, r8
	mov F4, r9

	mov rax, [rbp + 48]		;F5
	mov F5, rax

	mov rax, [rbp + 56]		;F6
	mov F6, rax

	pop rbp

	ret

InitASMDemag_FFTSpaces ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_Kernels(ReIm *Kx, ReIm *Ky, ReIm *Kz, ReIm *Kxy, ReIm *Kxz, ReIm *Kyz);}

PUBLIC InitASMDemag_Kernels
InitASMDemag_Kernels PROC

	push rbp
	mov rbp, rsp

	mov Kx, rcx
	mov Ky, rdx
	mov Kz, r8
	mov Kxy, r9

	mov rax, [rbp + 48]		;Kxz
	mov Kxz, rax

	mov rax, [rbp + 56]		;Kyz
	mov Kyz, rax

	pop rbp

	ret

InitASMDemag_Kernels ENDP

;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMMeshSpaces(double *Mx, double *My, double *Mz, double *Heffx, double *Heffy, double *Heffz);}

PUBLIC InitASMMeshSpaces
InitASMMeshSpaces PROC

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

InitASMMeshSpaces ENDP


;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_LUT(ReIm *cossin, ReIm *sincos, int *fftShuffle, int *fftShuffleInv, int _MAXFFT);}

PUBLIC InitASMDemag_LUT
InitASMDemag_LUT PROC

	push rbp
	mov rbp, rsp

	mov cossin, rcx
	mov sincos, rdx
	mov fftShuf, r8
	mov fftShufInv, r9
	
	mov rax, [rbp + 48]
	mov MAXFFT, rax

	pop rbp

	ret

InitASMDemag_LUT ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_sizes(int nx, int ny, int nz, int N, int M, int K);}

PUBLIC InitASMDemag_sizes
InitASMDemag_sizes PROC

	push rbp
	mov rbp, rsp

	mov nx, rcx				;nx
	mov ny, rdx				;ny
	mov nz, r8				;nz	
	
	imul rcx, rdx
	imul rcx, r8			;nx * ny * nz
	mov nxnynz, rcx

	mov N, r9				;N
	
	mov rcx, r9				;N
	
	mov rax, r9				;N
	call GetPow2
	mov powN, rax

	mov rax, [rbp + 48]		;M
	mov M, rax
	
	imul rcx, rax			;N*M

	call GetPow2
	mov powM, rax

	mov rax, [rbp + 56]		;K
	mov K, rax
	
	imul rcx, rax			;N*M*K
	
	mov NMK, rcx

	call GetPow2
	mov powK, rax

	mov rsp, rbp
	pop rbp

	ret

InitASMDemag_sizes ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void __SaveXMMRegs(void);}

PUBLIC __SaveXMMRegs
__SaveXMMRegs PROC

	movupd _xmm0, xmm0
	movupd _xmm1, xmm1
	movupd _xmm2, xmm2
	movupd _xmm3, xmm3
	movupd _xmm4, xmm4
	movupd _xmm5, xmm5
	movupd _xmm6, xmm6
	movupd _xmm7, xmm7
	movupd _xmm8, xmm8
	movupd _xmm9, xmm9
	movupd _xmm10, xmm10
	movupd _xmm11, xmm11
	movupd _xmm12, xmm12
	movupd _xmm13, xmm13
	movupd _xmm14, xmm14
	movupd _xmm15, xmm15

	ret

__SaveXMMRegs ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void __RestoreXMMRegs(void);}

PUBLIC __RestoreXMMRegs
__RestoreXMMRegs PROC

	movupd xmm0, _xmm0
	movupd xmm1, _xmm1
	movupd xmm2, _xmm2
	movupd xmm3, _xmm3
	movupd xmm4, _xmm4
	movupd xmm5, _xmm5
	movupd xmm6, _xmm6
	movupd xmm7, _xmm7
	movupd xmm8, _xmm8
	movupd xmm9, _xmm9
	movupd xmm10, _xmm10
	movupd xmm11, _xmm11
	movupd xmm12, _xmm12
	movupd xmm13, _xmm13
	movupd xmm14, _xmm14
	movupd xmm15, _xmm15

	ret

__RestoreXMMRegs ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;INTERNAL USE METHODS;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;Get log2 of number which is a positive power of 2. Input in rax, output also in rax.
;
GetPow2 PROC

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

GetPow2 ENDP

_TEXT   ENDS
END

