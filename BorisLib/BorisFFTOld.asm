; DEPRECATED
;
; SSE FFT methods for BORIS
; Author : Serban Lepadatu, SLepadatu@uclan.ac.uk
;
; Use x64 architecture and fastcall convention for external call routines.
;
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------
; MASM use routines :
; 1. FFT_1D_Radix4R_ZP_ASM_N :  
;							r8 input pointer, r9 output pointer, r10 calculation space pointer, rsi = cossin pointer, rax = 0, rbx = 16, rdx = N/2 (N must be even power of 2), r13 = Nzp * 8
;	 FFT_1D_Radix4R_ZP_Odd_ASM_S : Called by FFT_1D_Radix4R_ZP_ASM_N
;
; 2. FFT_1D_Radix4R_ZP_Odd_ASM_N : 
;							r8 input pointer, r9 output pointer, r10 calculation space pointer, rsi = cossin pointer, rax = 0, rbx = 16, rdx = N/2 (N must be odd power of 2), r13 = Nzp * 8
;	
;	_FFT_1D_Radix4R_ZP_ASM_N, _FFT_1D_Radix4R_ZP_ASM_S : Called by FFT_1D_Radix4R_ZP_Odd_ASM_N and FFT_1D_Radix4R_ZP_Odd_ASM_S
;	
; 3. ReIFFT_Prepare : r8 input pointer, rdx : N/2
;
; 4. IFFT_1D_Radix4_ASM_N : 
;							r8 input pointer, r9 output pointer, r10 calculation space pointer, rsi = cossin pointer, rax = 0, rbx = 16, rdx = N (must be even power of 2)
; 5. IFFT_1D_Radix4_Odd_ASM_N : 
;							r8 input pointer, r9 output pointer, r10 calculation space pointer, rsi = cossin pointer, rax = 0, rbx = 16, rdx = N (must be odd power of 2)
; Note : For length N data obtainable by FFT from real data, first call ReIFFT_Prepare with rdx set to N/2. After this call IFFT_1D_Radix4_ASM_N (or Odd) with same value in
;		 rdx as this will then do a half-length complex IFFT on the N data points
;
; FFT_1D_Radix4_ASM_N, FFT_1D_Radix4_ASM_S : radix-4 FFT used by the IFFT routines
;-------------------------------------------------------------------------------------------------------------------------------------------------------------------

.DATA

;pointers, data must be 16-byte aligned
FF		DQ ?
FC		DQ ?

;twiddle factors, 16byte aligned pairs
cossin	DQ ?
;real fft factors
refftc	DQ ?

;dimensions
MAXFFT	DD 8192		

nx		DD ?
ny		DD ?
nz		DD ?

nxnynz	DD ?

N		DD ?
M		DD ?
K		DD ?

NMK		DD ?
magpoints DD ?

;mask for even or odd FFT selection : 
;bit 0 : N even power of 2?
;bit 1 : M even power of 2?
;bit 2 : K even power of 2?
evenodd DB 0

maxside	DD ?

.CODE

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Call Functions for Runtime Calculations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PUBLIC SaveXMMRegs
SaveXMMRegs PROC

	movupd [rcx], xmm0
	movupd [rcx+16], xmm1
	movupd [rcx+32], xmm2
	movupd [rcx+48], xmm3
	movupd [rcx+64], xmm4
	movupd [rcx+80], xmm5
	movupd [rcx+96], xmm6
	movupd [rcx+112], xmm7
	movupd [rcx+128], xmm15

	ret

SaveXMMRegs ENDP

PUBLIC RestoreXMMRegs
RestoreXMMRegs PROC

	movupd xmm0, [rcx]
	movupd xmm1, [rcx+16]
	movupd xmm2, [rcx+32]
	movupd xmm3, [rcx+48]
	movupd xmm4, [rcx+64]
	movupd xmm5, [rcx+80]
	movupd xmm6, [rcx+96]
	movupd xmm7, [rcx+112]
	movupd xmm15, [rcx+128]

	ret

RestoreXMMRegs ENDP

;---------------------------------------------------------------------------------------------------------------------------
;Input matrix is N*M with zero padding outside the nx*ny matrix
;Convolution obtained as follows (2D case):
;
; 1. 1D real FFTs of length N with zero padding from nx up on rows 0 to ny-1 and pack data in each row as :
; Re(0), Re(N/2), Re(1), Im(1), ...., Re(N/2-1), Im(N/2-1) - Note : remember symmetries of real FFT to understand why
;
; Notes : 
; Each thread has a FFT calculation space and a FFT space to handle each row (FC : maxside * 8 * 2 bytes per thread),
; thus not an in-place FFT : must leave the magnetization data unaltered.
;
; 2. Transpose each FFTed row into the output space by (Re, Im) pairs. Output space is considered as a M * N/2 array of 
;    complex numbers. The special pair (Re(0), Re(N/2)) considered as a (Re, Im) pair, thus the first row in the transposed
;    array consists of these special pairs.
;
; mx -> Heffx, my -> Heffy, mz -> Heffz
;
; Note : source is read linearly, only suffering from write cache misses - not too bad.
;	Should be faster than an in-place matrix transposition, although this was not tested.
;
; 3. 1D complex FFTs of length M with zero padding from ny up on rows 0 to N/2-1 of the transposed matrix
;
; Note : 1D complex FFT also done on the first row of the transposed matrix.
;        These are in-place FFTs.
;        For 3D convolution a further transpose step and complex FFTs : easy extension of 2D case with one exception - read endnote
;
; Heffx -> Heffx, Heffy -> Heffy, Heffz -> Heffz
;
; 4. Complex multiplication of FFTed magnetization matrices with the pre-FFTed kernels.
;    Note : straightforward complex multiplication only for rows 1 to N/2-1
;           First row needs a more complicated multiplication - see demag.cpp
;           Multiplication implemented in C++, the gain from implementing in ASM is small
;           as most of the time is spent in accesing memory (which is why the parallelization efficiency
;           is so bad on this part of the algorithm).
;
; Heffx -> Heffx, Heffy -> Heffy, Heffz -> Heffz
;
; 5. 1D complex IFFTs of length M on rows 0 to N/2-1. Transpose each FFTed row by (Re, Im) pairs. 
;    Thus a row of M complex numbers becomes two columns of M real numbers each.
;
; Heffx -> FF; Heffy -> FF; Heffz -> FF 
;
; 6. 1D complex IFFT of length N/2 on rows 0 to ny-1. nx values from each IFFTed row is copied into the final output space.
;    Before doing the complex IFFT the values are prepared : classical method of half-length complex IFFT on data which
;    may be obtained by FFT from real data.
;
; FF -> Heffx; FF -> Heffy; FF -> Heffz
;
; Notes: Heffx, Heffy, Heffz must be arrays of size N*M*K but the final output is placed in the first nx*ny*nz qwords so
;        the rest of the program does not need to know about the (N, M, K) space, only when allocating memory. 
;        Calculating the effective field, the demag must be the first to be calculated, all the other contributions just add to it.
;
; Note on 3D extension: After the first set of real FFTs on rows, the special zeroth and (N/2)nd columns (the ones that have zero imaginary
;                       part as follows from real FFT properties) are kept as they are and the zero imaginary part written explicitly.
;                       Thus they are not packed in a single complex FFT input as for the 2D case but an extra complex FFT is used so they are computed
;                       separately. This makes for an easier multiplication with the kernel (pointwise complex multiplication throughout).
;						The transposition step from rows to columns is slightly modified so the (N/2)nd column is placed right after the zeroth column,
;                       each having the zeros in the imaginary parts written explicitly - called them columns but they are rows in the transposed matrix.
;                       When transposing back for the last set of IFFT the transposition step is modified to get rid of the zero imaginary part.
;                       The transposed matrices have (N/2+1) rows (each row has (Re,Im) pairs hence the 1/2).
;
;---------------------------------------------------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_2D_Magnetization_r(int thread, int chunksize, double *input, double *output);}
;FFT on rows
;
;rcx : thread number (0 to MPNP-1 where MPNP is the number of processors)
;rdx : chunk size (rows to do)
;r8 : input (e.g. mx)
;r9 : output (e.g. Heffx)
PUBLIC ReFFT_2D_Magnetization_r
ReFFT_2D_Magnetization_r PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov r12, r9							;set final output

	xor rax, rax
	xor rdi, rdi
	mov eax, maxside
	imul ecx, eax
	shl rcx, 4							;maxside * thread * 16
	shl rax, 3							;maxside * 8
	mov r9, FC
	add r9, rcx							;first thread calculation space
	mov r10, r9
	add r10, rax						;second thread calculation space

	;1D FFT along x
	mov rcx, rdx						;number of FFTs (should be <= ny)
	mov edx, N
	shr rdx, 1							;N/2 for real FFT
	
	mov eax, nx
	shl rax, 3							;nx * 8
	mov r13, rax						;Nzp = nx * 8
	
	mov rbx, 16							;del = 1 (*16)
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 1
	jz ReFFT_2D_Magnetization_r_oddN	;N is odd power of 2
	
	xor rax, rax						;off = 0

ReFFT_2D_Magnetization_r_evenN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1							;pairs of values moved
	mov edi, M							;current vertical length
	call Transpose2Dline
	
	add r8, r13							;increment nx*8 for source	
	add r12, 16							;final output space increment

	dec rcx
	jnz ReFFT_2D_Magnetization_r_evenN_Loop

	jmp ReFFT_2D_Magnetization_r_return

ReFFT_2D_Magnetization_r_oddN:

	xor rax, rax						;off = 0
	
ReFFT_2D_Magnetization_r_oddN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1							;pairs of values moved
	mov edi, M							;current vertical length
	call Transpose2Dline

	add r8, r13							;increment nx*8 for source	
	add r12, 16							;final output space increment

	dec rcx
	jnz ReFFT_2D_Magnetization_r_oddN_Loop

ReFFT_2D_Magnetization_r_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_2D_Magnetization_r ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_3D_Magnetization_r(int thread, int chunksize, double *input, double *output);}
;FFT on rows
;
;rcx : thread number (0 to MPNP-1 where MPNP is the number of processors)
;rdx : chunk size (rows to do)
;r8 : input (e.g. mx)
;r9 : output (FF)
PUBLIC ReFFT_3D_Magnetization_r
ReFFT_3D_Magnetization_r PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov r12, r9							;set final output

	xor rax, rax
	xor rdi, rdi
	mov eax, maxside
	imul ecx, eax
	shl rcx, 4							;maxside * thread * 16
	shl rax, 3							;maxside * 8
	mov r9, FC
	add r9, rcx							;first thread calculation space
	mov r10, r9
	add r10, rax						;second thread calculation space

	;1D FFT along x
	mov rcx, rdx						;number of FFTs (should be <= ny)
	mov edx, N
	shr rdx, 1							;N/2 for real FFT

	mov eax, nx
	shl rax, 3							;nx * 8
	mov r13, rax						;Nzp = nx * 8
	
	mov rbx, 16							;del = 1 (*16)
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 1
	jz ReFFT_3D_Magnetization_r_oddN	;N is odd power of 2
	
	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReFFT_3D_Magnetization_r_evenN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx							;save rows counter

ReFFT_3D_Magnetization_r_evenN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1							;pairs of values moved
	mov edi, M							;current vertical length
	call _Transpose2Dline
	
	add r8, r13							;increment nx*8 for source	
	add r12, 16							;final output space increment

	dec rcx
	jnz ReFFT_3D_Magnetization_r_evenN_Loop

	pop rcx
	pop r12
	pop r8
	
	xor rdi, rdi
	mov edi, nx
	imul edi, ny
	shl rdi, 3							;nx*ny*8
	add r8, rdi							;next plane for source
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M*16
	add r12, rdi						;next plane for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_Magnetization_r_evenN_Loop_planes

	jmp ReFFT_3D_Magnetization_r_return

ReFFT_3D_Magnetization_r_oddN:

	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReFFT_3D_Magnetization_r_oddN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx							;save rows counter

ReFFT_3D_Magnetization_r_oddN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1							;pairs of values moved
	mov edi, M							;current vertical length
	call _Transpose2Dline

	add r8, r13							;increment nx*8 for source	
	add r12, 16							;final output space increment

	dec rcx
	jnz ReFFT_3D_Magnetization_r_oddN_Loop

	pop rcx
	pop r12
	pop r8
	
	xor rdi, rdi
	mov edi, nx
	imul edi, ny
	shl rdi, 3							;nx*ny*8
	add r8, rdi							;next plane for source
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M*16
	add r12, rdi						;next plane for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_Magnetization_r_oddN_Loop_planes

ReFFT_3D_Magnetization_r_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_3D_Magnetization_r ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_2D_Magnetization_c(int thread, int chunksize, double *output);}
;FFT on columns
;
;rcx thread number
;rdx chunksize
;r8 output (e.g. Heffx)
PUBLIC ReFFT_2D_Magnetization_c
ReFFT_2D_Magnetization_c PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov r9, r8							;in-place FFT

	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rcx, rdx						;number of FFTs (<= N)
	shr rcx, 1							;divide by 2 for complex FFTs
	mov rbx, 16							;del = 1 * 16
	mov edx, ny
	mov r13, rdx
	shl r13, 4							;Nzp = ny * 16
	mov edx, M							;FFT length : complex FFT of length M
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 2
	jz ReFFT_2D_Magnetization_c_oddM	;M is odd power of 2

	xor rax, rax						;off = 0

ReFFT_2D_Magnetization_c_evenM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_ASM_N
	pop rcx
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	add r9, r11							

	dec rcx
	jnz ReFFT_2D_Magnetization_c_evenM_Loop

	jmp ReFFT_2D_Magnetization_c_return

ReFFT_2D_Magnetization_c_oddM:

	xor rax, rax						;off = 0

ReFFT_2D_Magnetization_c_oddM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_Odd_ASM_N
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	add r9, r11							

	dec rcx
	jnz ReFFT_2D_Magnetization_c_oddM_Loop

ReFFT_2D_Magnetization_c_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_2D_Magnetization_c ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_3D_Magnetization_c(int thread, int chunksize, double *input, double *output);}
;FFT on columns
;
;rcx thread number
;rdx chunksize
;r8 input (FF)
;r9 output (Heffx)
PUBLIC ReFFT_3D_Magnetization_c
ReFFT_3D_Magnetization_c PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov rsi, rcx						;save thread number

	mov r12, r9							;set final output
	mov r9, r8							;in-place FFT

	xor rax, rax
	xor rdi, rdi
	mov eax, maxside
	imul ecx, eax
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	;1D FFT along x
	mov rcx, rdx						;number of FFTs (should be <= N/2+1)
	mov edx, M							;FFT length : complex FFT of length M
	
	test esi, esi
	jnz ReFFT_3D_Magnetization_c_not0

	inc rcx								;extra FFT for thread zero
	jmp ReFFT_3D_Magnetization_c_0

ReFFT_3D_Magnetization_c_not0:			;not thread zero

	mov rsi, rdx
	shl rsi, 4							;M*16
	add r8, rsi
	mov r9, r8							;if not thread zero then increment start by one row	
	mov esi, K
	shl rsi, 4							;K*16
	add r12, rsi

ReFFT_3D_Magnetization_c_0:

	mov eax, ny
	shl rax, 4							;ny * 16
	mov r13, rax						;Nzp = ny * 16
	
	mov rbx, 16							;del = 1 (*16)
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 2
	jz ReFFT_3D_Magnetization_c_oddN	;M is odd power of 2
	
	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReFFT_3D_Magnetization_c_evenN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx							;save rows counter

ReFFT_3D_Magnetization_c_evenN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4_ZP_ASM_N
	pop rcx
	pop r12

	mov eax, M							;points to transpose
	mov edi, N							
	shr edi, 1
	inc edi
	imul edi, K							;plane : (N/2+1)*K * 2
	shl edi, 1
	call Transpose3Dline
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11							
	mov r9, r8
	mov edi, K
	shl rdi, 4							;K*16
	add r12, rdi						;final output space increment

	dec rcx
	jnz ReFFT_3D_Magnetization_c_evenN_Loop

	pop rcx
	pop r12
	pop r8
	
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M*16
	add r8, rdi							;next plane for source
	mov r9, r8
	add r12, 16							;next plane for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_Magnetization_c_evenN_Loop_planes

	jmp ReFFT_3D_Magnetization_c_return

ReFFT_3D_Magnetization_c_oddN:

	xor rax, rax						;off = 0
	mov edi, nz							;number of planes

ReFFT_3D_Magnetization_c_oddN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx

ReFFT_3D_Magnetization_c_oddN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4_ZP_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, M							;points to transpose
	mov edi, N							
	shr edi, 1
	inc edi
	imul edi, K							;plane : (N/2+1)*K * 2
	shl edi, 1
	call Transpose3Dline

	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11							
	mov r9, r8
	mov edi, K
	shl rdi, 4							;K*16
	add r12, rdi						;final output space increment

	dec rcx
	jnz ReFFT_3D_Magnetization_c_oddN_Loop

	pop rcx
	pop r12
	pop r8
	
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M*16
	add r8, rdi							;next plane for source
	mov r9, r8
	add r12, 16							;next plane for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_Magnetization_c_oddN_Loop_planes

ReFFT_3D_Magnetization_c_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_3D_Magnetization_c ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_3D_Magnetization_h(int thread, int chunksize, double *output);}
;FFT on heights
;
;rcx thread number
;rdx chunksize
;r8 input/output (e.g. Heffx)
PUBLIC ReFFT_3D_Magnetization_h
ReFFT_3D_Magnetization_h PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov rsi, rcx						;save thread number

	mov r9, r8							;in-place FFT

	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rcx, rdx						;number of FFTs (<= N/2+1)
	mov rbx, 16							;del = 1 * 16
	mov edx, nz
	mov r13, rdx
	shl r13, 4							;Nzp = nz * 16
	mov edx, K							;FFT length : complex FFT of length K

	test esi, esi
	jnz ReFFT_3D_Magnetization_h_not0

	inc rcx								;one more FFT to do for thread 0
	jmp ReFFT_3D_Magnetization_h_0

ReFFT_3D_Magnetization_h_not0:

	mov rsi, rdx
	shl rsi, 4
	add r8, rsi
	mov r9, r8

ReFFT_3D_Magnetization_h_0:

	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 4
	jz ReFFT_3D_Magnetization_h_oddM	;K is odd power of 2

	xor rax, rax						;off = 0
	mov edi, M							;number of planes to do
	
ReFFT_3D_Magnetization_h_evenM_Loop_planes:

	push rdi
	push r8								;save source
	push rcx							;save rows counter

ReFFT_3D_Magnetization_h_evenM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_ASM_N
	pop rcx
	
	mov r11, rdx
	shl r11, 4							;increment is K * 16
	add r8, r11
	mov r9, r8

	dec rcx
	jnz ReFFT_3D_Magnetization_h_evenM_Loop

	pop rcx
	pop r8
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, K
	shl rdi, 4							;K*(N/2+1)*16
	add r8, rdi							;next plane
	mov r9, r8

	pop rdi
	dec rdi
	jnz ReFFT_3D_Magnetization_h_evenM_Loop_planes

	jmp ReFFT_3D_Magnetization_h_return

ReFFT_3D_Magnetization_h_oddM:

	xor rax, rax						;off = 0
	mov edi, M							;number of planes to do
	
ReFFT_3D_Magnetization_h_oddM_Loop_planes:

	push rdi
	push r8								;save source
	push rcx

ReFFT_3D_Magnetization_h_oddM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_Odd_ASM_N
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is K * 16
	add r8, r11
	mov r9, r8						

	dec rcx
	jnz ReFFT_3D_Magnetization_h_oddM_Loop

	pop rcx
	pop r8
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, K
	shl rdi, 4							;K*(N/2+1)*16
	add r8, rdi							;next plane
	mov r9, r8

	pop rdi
	dec rdi
	jnz ReFFT_3D_Magnetization_h_oddM_Loop_planes

ReFFT_3D_Magnetization_h_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_3D_Magnetization_h ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReIFFT_3D_Field_h(int thread, int chunksize, double *input, double *output);}
;IFFT on heights
;
;rcx : thread number
;rdx : chunksize
;r8 : input (e.g. Heffx)
;r9 : output (FF)
PUBLIC ReIFFT_3D_Field_h
ReIFFT_3D_Field_h PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov rsi, rcx						;save thread number

	mov r12, r9							;final output
	mov r9, r8							;in-place IFFT

	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rbx, 16							;del = 1 * 16
	mov rcx, rdx						;number of IFFTs (<=N/2+1)	
	mov edx, K							;FFT length : complex IFFT of length K

	test esi, esi
	jnz ReIFFT_3D_Field_h_not0

	inc rcx								;one more IFFT to do for thread zero
	jmp ReIFFT_3D_Field_h_0

ReIFFT_3D_Field_h_not0:

	mov rsi, rdx
	shl rsi, 4							;K*16
	add r8, rsi
	mov r9, r8
	mov esi, M
	shl rsi, 4							;M*16
	add r12, rsi

ReIFFT_3D_Field_h_0:

	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 4
	jz ReIFFT_3D_Field_h_oddM			;K is odd power of 2

	xor rax, rax						;off = 0
	mov edi, M							;number of planes
	
ReIFFT_3D_Field_h_evenM_Loop_planes:

	push rdi
	push r8								;save source
	push r12							;save output
	push rcx							;save rows counter

ReIFFT_3D_Field_h_evenM_Loop:

	push r12
	push rcx
	call IFFT_1D_Radix4_ASM_N
	pop rcx
	pop r12

	mov eax, nz							;points to transpose
	mov edi, N							
	shr edi, 1
	inc edi
	imul edi, M							;plane : (N/2+1)*M * 2
	shl edi, 1
	call Transpose3Dline
	
	mov r11, rdx
	shl r11, 4							;increment is K * 16
	add r8, r11
	mov r9, r8
	mov edi, M
	shl rdi, 4							;M*16
	add r12, rdi						;final output space increment

	dec rcx
	jnz ReIFFT_3D_Field_h_evenM_Loop

	pop rcx
	pop r12
	pop r8

	mov edi, N
	shr edi, 1
	inc edi
	imul edi, K
	shl rdi, 4							;(N/2+1)*K * 16
	add r8, rdi
	mov r9, r8
	add r12, 16							;next plane for output

	pop rdi
	dec rdi
	jnz ReIFFT_3D_Field_h_evenM_Loop_planes

	jmp ReIFFT_3D_Field_h_return

ReIFFT_3D_Field_h_oddM:

	xor rax, rax						;off = 0
	mov edi, M							;number of planes
	
ReIFFT_3D_Field_h_oddM_Loop_planes:

	push rdi
	push r8								;save source
	push r12							;save output
	push rcx

ReIFFT_3D_Field_h_oddM_Loop:

	push r12
	push rcx
	call IFFT_1D_Radix4_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, nz							;points to transpose
	mov edi, N							
	shr edi, 1
	inc edi
	imul edi, M							;plane : (N/2+1)*M * 2
	shl edi, 1
	call Transpose3Dline
	
	mov r11, rdx
	shl r11, 4							;increment is K * 16
	add r8, r11
	mov r9, r8
	mov edi, M
	shl rdi, 4							;M*16
	add r12, rdi						;final output space increment

	dec rcx
	jnz ReIFFT_3D_Field_h_oddM_Loop

	pop rcx
	pop r12
	pop r8

	mov edi, N
	shr edi, 1
	inc edi
	imul edi, K
	shl rdi, 4							;(N/2+1)*K * 16
	add r8, rdi
	mov r9, r8
	add r12, 16							;next plane for output

	pop rdi
	dec rdi
	jnz ReIFFT_3D_Field_h_oddM_Loop_planes

ReIFFT_3D_Field_h_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReIFFT_3D_Field_h ENDP

;---------------------------------------------------------------------------------------------------------------------------
; 
;  Before doing the final IFFT on rows, must call the ReIFFT_Prepare routine for each row. This reverses the transformation introduced when
;  calculating the FFT of real data using a half-length complex FFT.
;
; e.g. let xn be the real input data, Xk = SUM(n) (xn * WNnk) is the FFT of xn
;
; Let zn = x(2n) + i*x(2n+1) and Zk = FFT(zn) the half-length complex FFT
;
; Then Xk = (1 -iWNk)/2 * Zk + (1+WNk)/2 * Z*(N/2-k) for k in 1 to N/2-1 and X0 = ReZ0 + ImZ0, XN/2 = ReZ0 - ImZ0
;
; The reverse of this is as follows :
;
; Given data Xk obtained by DFT from a real sequence xn first transform to Zk then do a half-length complex IFFT.
; Thus ReIFFT_Prepare does the transformation Xk to Zk using these formulas :
;
; ReZk = c1*ReXk + c2*ReXN/2-k - c3*(ImXk + ImXN/2-k)
; ImZk = c1*ImXk - c2*ImXN/2-k + c3*(ReXk - ReXN/2-k)  for k in 1 to N/2-1. Note Zk and ZN/2-k done in pairs as they use same input
; ReZ0 = (ReXn/2 + ReX0)/2
; ImZ0 = (ReX0 - ReXn/2)/2
;
; where c1, c2, c3 are the refftc factors i.e. c1 = (1-sin(2pi n/N))/2, c2 = (1+sin(2pi n/N))/2, c3 = cos(2pi n/N)/2
;
;---------------------------------------------------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReIFFT_2D_Field_c(int thread, int chunksize, double *input, double *output);}
;IFFT on columns
;
;rcx : thread number
;rdx : chunksize
;r8 : input
;r9 : output
PUBLIC ReIFFT_2D_Field_c
ReIFFT_2D_Field_c PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov r12, r9							;final output
	mov r9, r8							;in-place IFFT

	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rsi, cossin						;cossin pointer
	mov rbx, 16							;del = 1 * 16
	mov rcx, rdx						;number of IFFTs (<=N)	
	shr rcx, 1							;divide by 2 for complex IFFTs
	mov edx, M							;FFT length : complex IFFT of length M

	mov al, evenodd
	and al, 2
	jz ReIFFT_2D_Field_c_oddM			;M is odd power of 2

	xor rax, rax						;off = 0

ReIFFT_2D_Field_c_evenM_Loop:

	push r12
	push rcx
	call IFFT_1D_Radix4_ASM_N
	pop rcx
	pop r12

	mov eax, ny							;points to transpose
	mov edi, N							;current vertical length
	shr rdi, 1
	call Transpose2Dline
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	add r9, r11
	add r12, 16

	dec rcx
	jnz ReIFFT_2D_Field_c_evenM_Loop

	jmp ReIFFT_2D_Field_c_return

ReIFFT_2D_Field_c_oddM:

	xor rax, rax						;off = 0

ReIFFT_2D_Field_c_oddM_Loop:

	push r12
	push rcx
	call IFFT_1D_Radix4_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, ny						;points to transpose
	mov edi, N						;current vertical length
	shr rdi, 1
	call Transpose2Dline
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	add r9, r11
	add r12, 16

	dec rcx
	jnz ReIFFT_2D_Field_c_oddM_Loop

ReIFFT_2D_Field_c_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReIFFT_2D_Field_c ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReIFFT_3D_Field_c(int thread, int chunksize, double *input, double *output);}
;IFFT on columns
;
;rcx : thread number
;rdx : chunksize
;r8 : input (FF)
;r9 output (e.g. Heffx)
PUBLIC ReIFFT_3D_Field_c
ReIFFT_3D_Field_c PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov rsi, rcx						;save thread number

	mov r12, r9							;final output
	mov r9, r8							;in-place IFFT

	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rbx, 16							;del = 1 * 16
	mov rcx, rdx						;number of IFFTs (<=N/2+1)	
	mov edx, M							;FFT length : complex IFFT of length M

	test esi, esi
	jnz ReIFFT_3D_Field_c_not0

	inc rcx								;one more to do for thread zero
	mov r13, 2							;first 2 rows are transposed differently
	jmp ReIFFT_3D_Field_c_0

ReIFFT_3D_Field_c_not0:

	mov r13, 0							;normal transpositions
	mov esi, edx
	shl rsi, 4							;M*16
	add r8, rsi
	mov r9, r8

ReIFFT_3D_Field_c_0:

	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 2
	jz ReIFFT_3D_Field_c_oddM			;M is odd power of 2

	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReIFFT_3D_Field_c_evenM_Loop_planes:

	push rdi
	push r8								;save source
	push r12							;save output
	push rcx							;save rows counter
	push r13							;save flag

ReIFFT_3D_Field_c_evenM_Loop:

	push r13
	push r12
	push rcx
	call IFFT_1D_Radix4_ASM_N
	pop rcx
	pop r12

	mov eax, ny							;points to transpose
	mov edi, N							;current vertical length
	shr rdi, 1						
	pop r13
	call Transpose2Dline_
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	mov r9, r8
	add r12, 16							;final output space increment

	dec rcx
	jnz ReIFFT_3D_Field_c_evenM_Loop

	pop r13
	pop rcx
	pop r12
	pop r8

	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M * 16
	add r8, rdi
	mov r9, r8
	add r12, rdi						;next plane for output
	mov edi, M
	shl edi, 4							;M*16
	sub r12, rdi						;N*M*8 increment for output

	pop rdi
	dec rdi
	jnz ReIFFT_3D_Field_c_evenM_Loop_planes

	jmp ReIFFT_3D_Field_c_return

ReIFFT_3D_Field_c_oddM:

	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReIFFT_3D_Field_c_oddM_Loop_planes:

	push rdi
	push r8								;save source
	push r12							;save output
	push rcx
	push r13							;save flag

ReIFFT_3D_Field_c_oddM_Loop:

	push r13
	push r12
	push rcx
	call IFFT_1D_Radix4_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, ny							;points to transpose
	mov edi, N							;current vertical length
	shr rdi, 1						
	pop r13
	call Transpose2Dline_
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	mov r9, r8
	add r12, 16							;final output space increment

	dec rcx
	jnz ReIFFT_3D_Field_c_oddM_Loop

	pop r13				
	pop rcx
	pop r12
	pop r8

	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M * 16
	add r8, rdi
	mov r9, r8
	add r12, rdi						;next plane for output
	mov edi, M
	shl edi, 4							;M*16
	sub r12, rdi						;N*M*8 increment for output

	pop rdi
	dec rdi
	jnz ReIFFT_3D_Field_c_oddM_Loop_planes

ReIFFT_3D_Field_c_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReIFFT_3D_Field_c ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReIFFT_2D_Field_r(int thread, int chunksize, double *input, double *output);}
;IFFT on rows
;
;rcx : thread
;rdx : chunksize
;r8 : input
;r9 : output
PUBLIC ReIFFT_2D_Field_r
ReIFFT_2D_Field_r PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov r12, r9							;final output
	mov r9, r8							;in place IFFT
	
	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rsi, cossin						;cossin pointer
	mov rbx, 16
	mov rcx, rdx						;number of x IFFTs (<=ny)
	mov edx, N							;length of IFFT along x
	shr rdx, 1							;N/2

	mov al, evenodd
	and al, 1
	jz ReIFFT_2D_Field_r_oddN			;N is odd power of 2

	xor rax, rax						;off = 0
	
ReIFFT_2D_Field_r_evenN_Loop:

	push rcx
	push r12	
	call ReIFFT_Prepare
	call IFFT_1D_Radix4_Odd_ASM_N
	pop r12
	call CopyLine						;shift nx points from r9 to Heff space, incrementing r12 in the process
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is N*8
	add r8, r11
	add r9, r11

	dec rcx
	jnz ReIFFT_2D_Field_r_evenN_Loop

	jmp ReIFFT_2D_Field_r_return

ReIFFT_2D_Field_r_oddN:

	xor rax, rax						;off = 0

ReIFFT_2D_Field_r_oddN_Loop:

	push rcx
	push r12
	call ReIFFT_Prepare
	call IFFT_1D_Radix4_ASM_N
	pop r12
	call CopyLine						;shift nx points from r9 to Heff space, incrementing r12 in the process
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is N*8
	add r8, r11
	add r9, r11

	dec rcx
	jnz ReIFFT_2D_Field_r_oddN_Loop
	
ReIFFT_2D_Field_r_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReIFFT_2D_Field_r ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReIFFT_3D_Field_r(int thread, int chunksize, double *output);}
;IFFT on rows
;
;rcx : thread
;rdx : chunksize
;r8 : output
PUBLIC ReIFFT_3D_Field_r
ReIFFT_3D_Field_r PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	mov r9, r8							;in place IFFT
	
	xor rax, rax
	xor rdi, rdi
	imul ecx, maxside
	shl rcx, 4							;maxside * thread * 16
	mov r10, FC
	add r10, rcx						;thread calculation space

	mov rsi, cossin						;cossin pointer
	mov rbx, 16
	mov rcx, rdx						;number of x IFFTs (<=ny)
	mov edx, N							;length of IFFT along x
	shr rdx, 1							;N/2

	mov al, evenodd
	and al, 1
	jz ReIFFT_3D_Field_r_oddN			;N is odd power of 2

	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReIFFT_3D_Field_r_evenN_Loop_planes:

	push rdi
	push r8
	push rcx

ReIFFT_3D_Field_r_evenN_Loop:

	push rcx
	push r12	
	call ReIFFT_Prepare
	call IFFT_1D_Radix4_Odd_ASM_N
	pop r12
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is N*8
	add r8, r11
	mov r9, r8

	dec rcx
	jnz ReIFFT_3D_Field_r_evenN_Loop

	pop rcx
	pop r8
	mov edi, N
	imul edi, M
	shl rdi, 3
	add r8, rdi
	mov r9, r8

	pop rdi
	dec rdi
	jnz ReIFFT_3D_Field_r_evenN_Loop_planes

	jmp ReIFFT_3D_Field_r_return

ReIFFT_3D_Field_r_oddN:

	xor rax, rax						;off = 0
	mov edi, nz							;number of planes
	
ReIFFT_3D_Field_r_oddN_Loop_planes:

	push rdi
	push r8
	push rcx

ReIFFT_3D_Field_r_oddN_Loop:

	push rcx
	push r12
	call ReIFFT_Prepare
	call IFFT_1D_Radix4_ASM_N
	pop r12
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is N*8
	add r8, r11
	mov r9, r8

	dec rcx
	jnz ReIFFT_3D_Field_r_oddN_Loop
	
	pop rcx
	pop r8
	mov edi, N
	imul edi, M
	shl rdi, 3
	add r8, rdi
	mov r9, r8

	pop rdi
	dec rdi
	jnz ReIFFT_3D_Field_r_oddN_Loop_planes

ReIFFT_3D_Field_r_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReIFFT_3D_Field_r ENDP

;intput : r9
;output : r12, incrementing r12 in the process
CopyLine PROC

	mov r11, r9							;save source

	mov ecx, nx
	shr rcx, 1
	mov rdi, 16

CopyLine_loop:

	movapd xmm0, [r11]
	movapd [r12], xmm0

	add r11, rdi
	add r12, rdi
	
	dec rcx
	jnz CopyLine_loop

	ret

CopyLine ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_2D_K(double *kernel);}
;2D FFT on demag tensor
;rcx : input/output
PUBLIC ReFFT_2D_K
ReFFT_2D_K PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	push rcx							;save output address

	mov r12, FF							;set final output

	xor rax, rax
	xor rdi, rdi
	mov r10, FC
	mov r8, rcx
	xor rcx, rcx
	mov r9, r8

	;1D FFT along x
	mov ecx, M							;number of FFTs
	mov edx, N
	shr rdx, 1							;N/2 for real FFT
	
	mov eax, N
	shl rax, 3							;nx * 8
	mov r13, rax						;Nzp = nx * 8
	
	mov rbx, 16							;del = 1 (*16)
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 1
	jz ReFFT_2D_K_r_oddN				;N is odd power of 2
	
	xor rax, rax						;off = 0

ReFFT_2D_K_r_evenN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1
	mov edi, M							;current vertical length
	call Transpose2Dline
	
	add r8, r13							;increment N*8 for source	
	add r9, r13
	add r12, 16							;final output space

	dec rcx
	jnz ReFFT_2D_K_r_evenN_Loop

	jmp ReFFT_2D_K_r_return

ReFFT_2D_K_r_oddN:

	xor rax, rax						;off = 0
	
ReFFT_2D_K_r_oddN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1
	mov edi, M							;current vertical length
	call Transpose2Dline

	add r8, r13							;increment N*8 for source	
	add r9, r13
	add r12, 16

	dec rcx
	jnz ReFFT_2D_K_r_oddN_Loop

ReFFT_2D_K_r_return:

	mov r8, FF
	pop r9								;output

	xor rax, rax
	xor rdi, rdi
	mov r10, FC

	mov ecx, N							;number of FFTs
	shr rcx, 1							;divide by 2 for complex FFTs
	mov rbx, 16							;del = 1 * 16
	mov edx, M
	mov r13, rdx
	shl r13, 4							;Nzp = M * 16
	mov edx, M							;FFT length : complex FFT of length M
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 2
	jz ReFFT_2D_K_c_oddM				;M is odd power of 2

	xor rax, rax						;off = 0

ReFFT_2D_K_c_evenM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_ASM_N
	pop rcx
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	add r9, r11							

	dec rcx
	jnz ReFFT_2D_K_c_evenM_Loop

	jmp ReFFT_2D_K_c_return

ReFFT_2D_K_c_oddM:

	xor rax, rax						;off = 0

ReFFT_2D_K_c_oddM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_Odd_ASM_N
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11
	add r9, r11							

	dec rcx
	jnz ReFFT_2D_K_c_oddM_Loop

ReFFT_2D_K_c_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_2D_K ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extern "C" {void ReFFT_3D_K(double *kernel);}
;3D FFT on demag tensor
;rcx : input
PUBLIC ReFFT_3D_K
ReFFT_3D_K PROC

	push rbx
	push rcx
	push rdx
	push rdi
	push rsi
	push rbp
	push r8
	push r9
	push r12
	push r13
	push r14
	push r15

	push rcx

	mov r8, rcx
	mov r9, r8
	mov r12, FF							;set final output

	xor rax, rax
	xor rdi, rdi
	xor rcx, rcx
	mov r10, FC

	;1D FFT along x
	mov ecx, M							;number of FFTs
	mov edx, N
	shr rdx, 1							;N/2 for real FFT

	mov eax, N
	shl rax, 3							;N * 8
	mov r13, rax						;Nzp = N * 8
	
	mov rbx, 16							;del = 1 (*16)
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 1
	jz ReFFT_3D_K_r_oddN				;N is odd power of 2
	
	xor rax, rax						;off = 0
	mov edi, K							;number of planes
	
ReFFT_3D_K_r_evenN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx							;save rows counter

ReFFT_3D_K_r_evenN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1							;pairs of values moved
	mov edi, M							;current vertical length
	call _Transpose2Dline
	
	add r8, r13							;increment N*8 for source	
	mov r9, r8
	add r12, 16							;final output space increment

	dec rcx
	jnz ReFFT_3D_K_r_evenN_Loop

	pop rcx
	pop r12
	pop r8
	
	xor rdi, rdi
	mov edi, N
	imul edi, M
	shl rdi, 3							;N*M*8
	add r8, rdi							;next plane for source	
	mov r9, r8
	add r12, rdi						;next plane for destination
	mov edi, M
	shl edi, 4							;M*16
	add r12, rdi						;(N/2+1)*M*16 increment for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_K_r_evenN_Loop_planes

	jmp ReFFT_3D_K_r_return

ReFFT_3D_K_r_oddN:

	xor rax, rax						;off = 0
	mov edi, K							;number of planes
	
ReFFT_3D_K_r_oddN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx							;save rows counter

ReFFT_3D_K_r_oddN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4R_ZP_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, N							;points to transpose
	shr rax, 1							;pairs of values moved
	mov edi, M							;current vertical length
	call _Transpose2Dline

	add r8, r13							;increment N*8 for source	
	mov r9, r8
	add r12, 16							;final output space increment

	dec rcx
	jnz ReFFT_3D_K_r_oddN_Loop

	pop rcx
	pop r12
	pop r8
	
	xor rdi, rdi
	mov edi, N
	imul edi, M
	shl rdi, 3							;N*M*8
	add r8, rdi							;next plane for source	
	mov r9, r8
	add r12, rdi						;next plane for destination
	mov edi, M
	shl edi, 4							;M*16
	add r12, rdi						;(N/2+1)*M*16 increment for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_K_r_oddN_Loop_planes

ReFFT_3D_K_r_return:

	pop r12								;set final output
	sub rsp, 8
	mov r8, FF
	mov r9, r8							;in-place FFT

	xor rax, rax
	xor rdi, rdi
	xor rcx, rcx
	mov r10, FC

	;1D FFT along x
	mov ecx, N							;number of FFTs (should be <= N)
	shr rcx, 1							;N/2
	inc rcx								;N/2+1 rows
	mov edx, M							;FFT length : complex FFT of length M

	mov eax, M
	shl rax, 4							;M * 16
	mov r13, rax						;Nzp = M * 16
	
	mov rbx, 16							;del = 1 (*16)
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 2
	jz ReFFT_3D_K_c_oddN				;M is odd power of 2
	
	xor rax, rax						;off = 0
	mov edi, K							;number of planes
	
ReFFT_3D_K_c_evenN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx							;save rows counter

ReFFT_3D_K_c_evenN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4_ZP_ASM_N
	pop rcx
	pop r12

	mov eax, M							;points to transpose
	mov edi, N							
	shr edi, 1
	inc edi
	imul edi, K							;plane : (N/2+1)*K * 2
	shl edi, 1
	call Transpose3Dline
	
	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11							
	mov r9, r8
	mov edi, K
	shl rdi, 4							;K*16
	add r12, rdi						;final output space increment

	dec rcx
	jnz ReFFT_3D_K_c_evenN_Loop

	pop rcx
	pop r12
	pop r8
	
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M*16
	add r8, rdi							;next plane for source
	mov r9, r8
	add r12, 16							;next plane for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_K_c_evenN_Loop_planes

	jmp ReFFT_3D_K_c_return

ReFFT_3D_K_c_oddN:

	xor rax, rax						;off = 0
	mov edi, K							;number of planes

ReFFT_3D_K_c_oddN_Loop_planes:

	push rdi
	push r8								;save input pointer
	push r12							;save output pointer
	push rcx

ReFFT_3D_K_c_oddN_Loop:
	
	push r12
	push rcx
	call FFT_1D_Radix4_ZP_Odd_ASM_N
	pop rcx
	pop r12

	mov eax, M							;points to transpose
	mov edi, N							
	shr edi, 1
	inc edi
	imul edi, K							;plane : (N/2+1)*K * 2
	shl edi, 1
	call Transpose3Dline

	mov r11, rdx
	shl r11, 4							;increment is M * 16
	add r8, r11							
	mov r9, r8
	mov edi, K
	shl rdi, 4							;K*16
	add r12, rdi						;final output space increment

	dec rcx
	jnz ReFFT_3D_K_c_oddN_Loop

	pop rcx
	pop r12
	pop r8
	
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, M
	shl rdi, 4							;(N/2+1)*M*16
	add r8, rdi							;next plane for source
	mov r9, r8
	add r12, 16							;next plane for destination

	pop rdi								;planes counter
	dec rdi
	jnz ReFFT_3D_K_c_oddN_Loop_planes

ReFFT_3D_K_c_return:

	pop r8
	mov r9, r8							;in-place FFT

	xor rax, rax
	xor rdi, rdi
	xor rcx, rcx
	mov r10, FC

	mov ecx, N							;number of FFTs (<= N)
	shr rcx, 1							;divide by 2 for complex FFTs
	inc rcx								;N/2+1
	mov rbx, 16							;del = 1 * 16
	mov edx, K							;FFT length : complex FFT of length K
	mov r13, rdx
	shl r13, 4							;Nzp = K * 16
	mov rsi, cossin						;cossin pointer

	mov al, evenodd
	and al, 4
	jz ReFFT_3D_K_h_oddM				;K is odd power of 2

	xor rax, rax						;off = 0
	mov edi, M							;number of planes to do
	
ReFFT_3D_K_h_evenM_Loop_planes:

	push rdi
	push r8								;save source
	push rcx							;save rows counter

ReFFT_3D_K_h_evenM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_ASM_N
	pop rcx
	
	mov r11, rdx
	shl r11, 4							;increment is K * 16
	add r8, r11
	mov r9, r8

	dec rcx
	jnz ReFFT_3D_K_h_evenM_Loop

	pop rcx
	pop r8
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, K
	shl rdi, 4							;K*(N/2+1)*16
	add r8, rdi							;next plane
	mov r9, r8

	pop rdi
	dec rdi
	jnz ReFFT_3D_K_h_evenM_Loop_planes

	jmp ReFFT_3D_K_h_return

ReFFT_3D_K_h_oddM:

	xor rax, rax						;off = 0
	mov edi, M							;number of planes to do
	
ReFFT_3D_K_h_oddM_Loop_planes:

	push rdi
	push r8								;save source
	push rcx

ReFFT_3D_K_h_oddM_Loop:

	push rcx
	call FFT_1D_Radix4_ZP_Odd_ASM_N
	pop rcx

	mov r11, rdx
	shl r11, 4							;increment is K * 16
	add r8, r11
	mov r9, r8						

	dec rcx
	jnz ReFFT_3D_K_h_oddM_Loop

	pop rcx
	pop r8
	mov edi, N
	shr edi, 1
	inc edi
	imul edi, K
	shl rdi, 4							;K*(N/2+1)*16
	add r8, rdi							;next plane
	mov r9, r8

	pop rdi
	dec rdi
	jnz ReFFT_3D_K_h_oddM_Loop_planes

ReFFT_3D_K_h_return:

	pop r15
	pop r14
	pop r13
	pop r12
	pop r9
	pop r8
	pop rbp
	pop rsi
	pop rdi
	pop rdx
	pop rcx
	pop rbx

	ret

ReFFT_3D_K ENDP

;;;;;;;;;;;;;;;;;;;;;Transpose Routines
;row<->column transpose routine, mapping a set of memory-contiguous double-float pairs to a set of pairs spaced rdi*16 bytes apart
;r9 : source
;r12 : output - must not be the same as the input
;rax : points to transpose from this line
;rdi : current vertical length - which will become the horizontal in the transposed array, thus this is the new spacing between points (e.g. M)
;Note : this routine leaves rax zero, needed for the routines calling this
Transpose2Dline PROC

	mov r11, r9					;source pointer
	mov rbp, r12				;destination pointer
	shl rdi, 4					;length*16

Transpose2Dline_loop:

	movapd xmm0, [r11]
	movapd [rbp], xmm0

	add rbp, rdi				;+=length(*16) for destination
	add r11, 16					;next point for source
	dec rax
	jnz Transpose2Dline_loop

	ret

Transpose2Dline ENDP

;modified forward row<->column transpose routine, mapping a set of memory-contiguous double-float pairs to a set of pairs spaced rdi*16 bytes apart
;except for the first two double-floats which are transposed each to a (double-float,0) pair.
;r9 : source
;r12 : output - must not be the same as the input
;rax : points to transpose from this line
;rdi : current vertical length - which will become the horizontal in the transposed array, thus this is the new spacing between points (e.g. M)
;Note : this routine leaves rax zero, needed for the routines calling this
_Transpose2Dline PROC

	mov r11, r9					;source pointer
	mov rbp, r12				;destination pointer
	shl rdi, 4					;length*16

	xorpd xmm1, xmm1
	movapd xmm0, [r11]
	movapd xmm2, xmm0

	shufpd xmm0, xmm1, 0
	shufpd xmm2, xmm1, 1
	
	movapd [rbp], xmm0
	add rbp, rdi
	movapd [rbp], xmm2
	add rbp, rdi
	add r11, 16
	dec rax 

_Transpose2Dline_loop:

	movapd xmm0, [r11]
	movapd [rbp], xmm0

	add rbp, rdi				;+=length(*16) for destination
	add r11, 16					;next point for source
	dec rax
	jnz _Transpose2Dline_loop

	ret

_Transpose2Dline ENDP

;modified back row<->column transpose routine, mapping a set of memory-contiguous double-float pairs to a set of pairs spaced rdi*16 bytes apart
;except for the first two rows which have their imaginary parts stripped.
;r9 : source
;r12 : output - must not be the same as the input
;rax : points to transpose from this line
;rdi : vertical length (e.g. M)
;changed registers : rax (zero at the end), rdi, rbp, r11, r14
;r13 : flag
Transpose2Dline_ PROC

	mov r11, r9					;source pointer
	mov rbp, r12				;destination pointer
	shl rdi, 4					;length*16

	test r13, r13
	jnz Transpose2Dline__loop2

Transpose2Dline__loop:

	movapd xmm0, [r11]
	movapd [rbp], xmm0

	add rbp, rdi				;+=length(*16) for destination
	add r11, 16					;next point for source
	dec rax
	jnz Transpose2Dline__loop

	ret

Transpose2Dline__loop2:

	mov r14, [r11]
	mov [rbp], r14

	add rbp, rdi
	add r11, 16
	dec rax
	jnz Transpose2Dline__loop2

	dec r13
	sub r12, 8					

	ret

Transpose2Dline_ ENDP

;Does basically the same thing as Transpose2Dline except the spacing is specified as the number of double-floats in a plane e.g. (N/2+1)*M*2
;r9 : source
;r12 : output - must not be the same as the input
;rax : points to transpose from this line
;rdi : plane size (e.g. N*M)
;changed registers : rax (zero at the end), rdi, rbp, r11, r14
Transpose3Dline PROC

	mov r11, r9					;source pointer
	mov rbp, r12				;destination pointer
	shl rdi, 3					;plane*8

Transpose3Dline_loop:

	movapd xmm0, [r11]
	movapd [rbp], xmm0

	add rbp, rdi				;+=plane(*8) for destination
	add r11, 16					;next point for source
	dec rax
	jnz Transpose3Dline_loop

	ret

Transpose3Dline ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;External Call Initialization Functions;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_FFTSpaces_Old(double *FF, double *FC);}

PUBLIC InitASMDemag_FFTSpaces_Old
InitASMDemag_FFTSpaces_Old PROC

	mov FF, rcx
	mov FC, rdx

	ret

InitASMDemag_FFTSpaces_Old ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_twiddles(ReIm *cossin, double *refftc, int _MAXFFT);}

PUBLIC InitASMDemag_twiddles
InitASMDemag_twiddles PROC

	mov cossin, rcx
	mov refftc, rdx
	mov rax, r8
	mov MAXFFT, eax

	ret

InitASMDemag_twiddles ENDP

;;;;;;;;;;;;;;;;;;;;;;extern "C" {void InitASMDemag_sizes_Old(int nx, int ny, int nz, int N, int M, int K, int magpoints, int maxside);}

PUBLIC InitASMDemag_sizes_Old
InitASMDemag_sizes_Old PROC

	push rbp
	mov rbp, rsp

	push rbx
	push rcx
	push rdx

	mov nx, ecx
	mov ny, edx
	mov rax, r8
	mov nz, eax
	
	imul ecx, edx
	imul ecx, eax			;nx * ny * nz
	mov nxnynz, ecx

	mov rax, r9
	mov N, eax
	push rax

	mov ebx, eax			;N

	mov rax, [rbp + 48]		;M
	mov M, eax
	imul ebx, eax			;N*M
	push rax

	mov rax, [rbp + 56]		;K
	mov K, eax
	imul ebx, eax			;N*M*K
	mov NMK, ebx
	push rax

	mov rax, [rbp + 64]		;magpoints
	mov magpoints, eax

	mov rbx, 4

InitASMDemag_sizes_Next:

	dec rbx
	jz InitASMDemag_sizes_return

	mov al, evenodd
	shl al, 1
	mov evenodd, al

	xor rdx, rdx
	xor rax, rax
	pop rcx					;get K, then M, then N

InitASMDemag_sizes_SHR:

	inc rax
	shr ecx, 1
	jz InitASMDemag_sizes_Next					;for size 1, keep it marked as odd

	cmp ecx, 1
	jne InitASMDemag_sizes_SHR					;keep dividing

	and rax, 1
	jnz InitASMDemag_sizes_Next					;if remainder then it's odd - keep it marked as odd

InitASMDemag_sizes_Even:

	or evenodd, 1								;mark it as even
	jmp InitASMDemag_sizes_Next

InitASMDemag_sizes_return:

	mov rax, [rbp + 72]		;maxside
	mov maxside, eax

	pop rdx
	pop rcx
	pop rbx

	mov rsp, rbp
	pop rbp

	ret

InitASMDemag_sizes_Old ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 ODD, REAL INPUT, ZERO PADDED;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;---------------------------------------------------------------------------------------------------------------------------
;radix 4 FFT, real input, complex output, length odd power of 2. Zero padding from Nzp to N-1
;does the usual radix-4 (length N/2) complex FFT on the input data then recombine
;Uses modified radix-4 routines for zero padding (artificially pruned input), _FFT_1D_Radix4R_ZP_ASM_S, _FFT_1D_Radix4R_ZP_ASM_N
;---------------------------------------------------------------------------------------------------------------------------
;r8 : input (data has to be 16-byte aligned and have even number of dfp values)
;r9 : output (as above)
;r10 : calculation space
;rsi : cossin table
;rax : offset (set to zero)
;rbx : del (set to 16)
;rdx : length N/2 (N must be odd power of 2)
;r13 : Nzp * 8 value (points Nzp to N-1 are taken as zero - Nzp must be even)
;registers changed : rcx, rdi, rbp, r11, r12, r14, r15
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5
FFT_1D_Radix4R_ZP_Odd_ASM_N PROC

	;the following call will use the values x(n*2) as real part and x(n*2+1) as imaginary part where r8 points to x
	;the fft of this will be in the first (N/2) Re|Im pairs of array pointed by r10 - calculation space
	call _FFT_1D_Radix4R_ZP_ASM_S

	;final output (fft of real input r8) in array pointed by r9 after the following transformation

	mov rdi, rax						;save rax
	mov rcx, rdx						;rcx has N/2
	xor rdx, rdx
	mov eax, MAXFFT						;MAXFFT
	idiv rcx					
	imul rax, 12						;t1 = (MAXFFT/N) * 3 * 8
	mov rdx, rcx						;restore rdx

	mov rbp, refftc						;m

	mov r11, r10						;incrementing input index
	mov r12, rcx
	imul r12, rbx
	add r12, r10						;decrementing input index

	mov r14, r9							;incrementing output index
		
	;k = 0 and k = N/2 case treated separately:

	movapd xmm0, [r11]					; ImFCx[0] | ReFCx[0]
	movapd xmm1, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	shufpd xmm0, xmm0, 1				; Re-Im | Re + Im
	movapd [r14], xmm0					; store first point as : Fx[N/2] | Fx[0]
	
	;the rest

	dec rcx

FFT_1D_Radix4R_ZP_Odd_ASM_N_Loop:

	add rbp, rax						;m+=t1
	add r11, rbx
	sub r12, rbx
	add r14, rbx

	movapd xmm0, [r11]					; ImFCx[k] | ReFCx[k]
	movapd xmm1, [r12]					; ImFCx[(N/2)-k] | ReFCx[(N/2)-k]
	movapd xmm2, xmm1
	addsubpd xmm2, xmm0					; SIm | DRe
	movddup xmm3, QWORD PTR [rbp]
	mulpd xmm3, xmm0					; refftc*ImFCx[k] | refftc*ReFCx[k]
	movddup xmm4, QWORD PTR [rbp+8]
	mulpd xmm4, xmm1					; refftc*ImFCx[(N/2)-k] | refftc*ReFCx[(N/2)-k]
	movddup xmm5, QWORD PTR [rbp+16]
	mulpd xmm5, xmm2					; refftc*SIm | refftc*DRe
	shufpd xmm3, xmm3, 1
	shufpd xmm4, xmm4, 1
	addsubpd xmm3, xmm4
	addpd xmm3, xmm5
	shufpd xmm3, xmm3, 1				; ImFx[k] | ReFx[k]
	movapd [r14], xmm3					; store Fx[k] as Im | Re pair

	dec rcx
	jnz FFT_1D_Radix4R_ZP_Odd_ASM_N_Loop

	mov rax, rdi						;restore rax
	ret

FFT_1D_Radix4R_ZP_Odd_ASM_N ENDP


;r9, r10 swapped
FFT_1D_Radix4R_ZP_Odd_ASM_S PROC

	cmp rdx, 2
	je FFT_1D_Radix4R_ZP_Odd_ASM_S_Base

	shl rbx, 1										;del*2, N/2
	shr rdx, 1
	call _FFT_1D_Radix4R_ZP_ASM_N					;Recurse FFT on even part (offset not changed)
	shr rbx, 1
	add rax, rbx									;offset for odd part
	shl rbx, 1
	call _FFT_1D_Radix4R_ZP_ASM_N					;Recurse FFT on odd part (offset increased by del)
	shr rbx, 1
	shl rdx, 1
	sub rax, rbx									;restore offset

	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11
	shr rcx, 1										;N/2
	mov r14, rcx
	imul r14, rbx									;store (N/2)*del
	
	add r11, r9										;i1 - FC[i1]
	mov r12, rax									
	add r12, r10									;o1 - F[o1]	
	add r14, rax									
	add r14, r10									;o2 - F[o2]
	mov r15, rsi									;m
	mov rbp, rbx
	shl rbp, 1										;del*2

	;k = 0 handled separately - no multiplications required
	movapd xmm0, [r11]								;ImFC[i1] | ReFC[i1]
	movapd xmm1, [r11+rbx]							;ImFC[i2] | ReFC[i2]

	movapd xmm2, xmm0
	addpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o1]
	subpd xmm2, xmm1
	movapd [r14], xmm2								;store F[o1]

	dec rcx											;rcx has loop count
	
	;start loop

FFT_1D_Radix4R_ZP_Odd_ASM_S_Loop:

	add r15, rdi									;m+=t
	add r11, rbp									;i1+=2*del
	add r12, rbx									;o1+=del
	add r14, rbx									;o2+=del	

	movddup xmm0, QWORD PTR [r15]					;cos	 | cos
	movddup xmm1, QWORD PTR [r15+8]					;sin	 | sin   
	movapd xmm2, [r11+rbx]							;ImFC[i2] | ReFC[i2] 
	movapd xmm3, [r11]								;ImFC[i1] | ReFC[i1]

	mulpd xmm1, xmm2								;ImFC[i2] * sin | ReFC[i2] * sin
	shufpd xmm2, xmm2, 1							;ReFC[i2] | ImFC[i2]
	mulpd xmm0, xmm2								;ReFC[i2] * cos | ImFC[i2] * cos

	addsubpd xmm0, xmm1								;tmp0 | tmp1
	shufpd xmm0, xmm0, 1							;tmp1 | tmp0

	movapd xmm1, xmm3
	addpd xmm1, xmm0								;ImF[o1] | ReF[o1]
	subpd xmm3, xmm0								;ImF[o2] | ReF[o2]

	movapd [r12], xmm1								;store F[o1]
	movapd [r14], xmm3								;store F[o2]

	dec rcx
	jnz FFT_1D_Radix4R_ZP_Odd_ASM_S_Loop						

	ret

FFT_1D_Radix4R_ZP_Odd_ASM_S_Base:

;2-point fft

	mov r11, r8
	add r11, rax									;i1 = off
	mov r12, r10
	add r12, rax									;o1 = off, o2 = i2
	mov r14, r12
	add r14, rbx									;o2 = off + del
	mov r15, rax

	cmp r15, r13
	jge FFT_1D_Radix4R_ZP_Odd_ASM_S_Base_Z1
	movapd xmm0, [r11]								; Imx[i1] | Rex[i1]
	add r15, rbx
	cmp r15, r13
	jge FFT_1D_Radix4R_ZP_Odd_ASM_S_Base_Z2
	movapd xmm1, [r11+rbx]							; Imx[i2] | Rex[i2]

	jmp FFT_1D_Radix4R_ZP_Odd_ASM_S_Base_Continue

FFT_1D_Radix4R_ZP_Odd_ASM_S_Base_Z1:
	xorpd xmm0, xmm0
FFT_1D_Radix4R_ZP_Odd_ASM_S_Base_Z2:
	xorpd xmm1, xmm1

FFT_1D_Radix4R_ZP_Odd_ASM_S_Base_Continue:

	movapd xmm2, xmm0
	addpd xmm0, xmm1								
	movapd [r12], xmm0								; store in F[o1]
	subpd xmm2, xmm1
	movapd [r14], xmm2								; store in F[o2]

FFT_1D_Radix4R_ZP_Odd_ASM_S_Return:

	ret
FFT_1D_Radix4R_ZP_Odd_ASM_S ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 FFT ROUTINES WITH ZERO PADDING;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r8 : input (data has to be 16-byte aligned and have even number of dfp values)
;r9 : output (as above)
;r10 : calculation space
;rsi : cossin table
;rax : offset (set to zero)
;rbx : del (set to 16)
;rdx : length N/2 (N must be odd power of 2)
;r13 : Nzp * 8 value (points Nzp to N-1 are taken as zero - Nzp must be even)
;registers changed : rcx, rdi, rbp, r11, r12, r14, r15
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5
_FFT_1D_Radix4R_ZP_ASM_N PROC

	cmp rdx, 4										;reached 4-point fft ?
	je FFT_1D_Radix4R_ZP_ASM_N_Base

	shl rbx, 2										;del*4, N/4
	shr rdx, 2
	call _FFT_1D_Radix4R_ZP_ASM_S					;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call _FFT_1D_Radix4R_ZP_ASM_S					;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call _FFT_1D_Radix4R_ZP_ASM_S					;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call _FFT_1D_Radix4R_ZP_ASM_S					;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx									;offset restored
	
	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11									
	
	add r11, r10									;set i
	mov r12, rax									
	add r12, r9										;set o
	
	shr rcx, 2										;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx									;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jz FFT_1D_Radix4R_ZP_ASM_N_Return

	xor r15, r15									

	;start loop

FFT_1D_Radix4R_ZP_ASM_N_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	add r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi									;save rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz FFT_1D_Radix4R_ZP_ASM_N_Loop						

	ret					

FFT_1D_Radix4R_ZP_ASM_N_Base:

;4-point FFT

	mov r14, rax									
	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r9										;set o

	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	cmp r14, r13									;off>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_N_Base_Z1				;if yes then zero input
	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r14, rbx
	cmp r14, r13									;off+del>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_N_Base_Z2				;if yes then zero input
	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r14, rbx
	cmp r14, r13									;off+2*del>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_N_Base_Z3				;if yes then zero input
	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r14, rbx
	cmp r14, r13									;off+3*del>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_N_Base_Z4				;if yes then zero input
	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	jmp FFT_1D_Radix4R_ZP_ASM_N_Base_Continue

FFT_1D_Radix4R_ZP_ASM_N_Base_Z1:
	xorpd xmm0, xmm0
FFT_1D_Radix4R_ZP_ASM_N_Base_Z2:
	xorpd xmm1, xmm1
FFT_1D_Radix4R_ZP_ASM_N_Base_Z3:
	xorpd xmm3, xmm3
FFT_1D_Radix4R_ZP_ASM_N_Base_Z4:
	xorpd xmm5, xmm5

FFT_1D_Radix4R_ZP_ASM_N_Base_Continue:

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

FFT_1D_Radix4R_ZP_ASM_N_Return:

	ret

_FFT_1D_Radix4R_ZP_ASM_N ENDP

;r9, r10 swapped
_FFT_1D_Radix4R_ZP_ASM_S PROC

	cmp rdx, 4										;reached 4-point fft ?
	je FFT_1D_Radix4R_ZP_ASM_S_Base

	shl rbx, 2										;del*4, N/4
	shr rdx, 2
	call _FFT_1D_Radix4R_ZP_ASM_N					;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call _FFT_1D_Radix4R_ZP_ASM_N					;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call _FFT_1D_Radix4R_ZP_ASM_N					;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call _FFT_1D_Radix4R_ZP_ASM_N					;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2									
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx									;offset restored

	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11									
	
	add r11, r9										;set i
	mov r12, rax									
	add r12, r10									;set o
	
	shr rcx, 2										;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx									;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jz FFT_1D_Radix4R_ZP_ASM_S_Return

	xor r15, r15									

	;start loop

FFT_1D_Radix4R_ZP_ASM_S_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	add r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;m+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;m+m+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz FFT_1D_Radix4R_ZP_ASM_S_Loop						

	ret

FFT_1D_Radix4R_ZP_ASM_S_Base:

;4-point FFT

	mov r14, rax									
	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r10									;set o

	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	cmp r14, r13									;off>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_S_Base_Z1				;if yes then zero input
	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r14, rbx
	cmp r14, r13									;off+del>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_S_Base_Z2				;if yes then zero input
	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r14, rbx
	cmp r14, r13									;off+2*del>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_S_Base_Z3				;if yes then zero input
	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r14, rbx
	cmp r14, r13									;off+3*del>=Nzp ?
	jge FFT_1D_Radix4R_ZP_ASM_S_Base_Z4				;if yes then zero input
	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	jmp FFT_1D_Radix4R_ZP_ASM_S_Base_Continue

FFT_1D_Radix4R_ZP_ASM_S_Base_Z1:
	xorpd xmm0, xmm0
FFT_1D_Radix4R_ZP_ASM_S_Base_Z2:
	xorpd xmm1, xmm1
FFT_1D_Radix4R_ZP_ASM_S_Base_Z3:
	xorpd xmm3, xmm3
FFT_1D_Radix4R_ZP_ASM_S_Base_Z4:
	xorpd xmm5, xmm5

FFT_1D_Radix4R_ZP_ASM_S_Base_Continue:

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]
	
FFT_1D_Radix4R_ZP_ASM_S_Return:

	ret

_FFT_1D_Radix4R_ZP_ASM_S ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 EVEN, REAL INPUT, ZERO PADDED;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;---------------------------------------------------------------------------------------------------------------------------
;radix 4 FFT, real input, complex output, length even power of 2. Zero padding from Nzp to N-1
;first does a half-length radix-2 complex fft (FFT_1D_Radix4R_ZP_Odd_ASM_S - which then breaks down into radix-4)
;then recombine at the end to obtain the output
;---------------------------------------------------------------------------------------------------------------------------
;r8 : input (data has to be 16-byte aligned and have even number of dfp values)
;r9 : output (as above)
;r10 : calculation space
;rsi : cossin table
;rax : offset (set to zero)
;rbx : del (set to 16)
;rdx : length N/2 (N must be even power of 2)
;r13 : Nzp * 8 value (points Nzp to N-1 are taken as zero - Nzp must be even)
;registers changed : rcx, rdi, rbp, r11, r12, r14, r15
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5
FFT_1D_Radix4R_ZP_ASM_N PROC

	;the following call will use the values x(n*2) as real part and x(n*2+1) as imaginary part where r8 points to x
	;the fft of this will be in the first (N/2) Re|Im pairs of array pointed by r10 - calculation space
	cmp rdx, 1
	je FFT_1D_Radix4R_ZP_ASM_N_2point
	call FFT_1D_Radix4R_ZP_Odd_ASM_S 

	;final output (fft of real input r8) in array pointed by r9 after the following transformation

	mov rdi, rax						;save rax
	mov rcx, rdx						;rcx has N/2
	xor rdx, rdx
	mov eax, MAXFFT						;MAXFFT
	idiv rcx					
	imul rax, 12						;t1 = (MAXFFT/N) * 3 * 8
	mov rdx, rcx						;restore rdx

	mov rbp, refftc						;m
			
	mov r11, r10						;incrementing input index
	mov r12, rcx
	imul r12, rbx
	add r12, r10						;decrementing input index

	mov r14, r9							;incrementing output index
	
	;k = 0 and k = N/2 case treated separately:

	movapd xmm0, [r11]					; ImFCx[0] | ReFCx[0]
	movapd xmm1, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	shufpd xmm0, xmm0, 1				; Re-Im | Re + Im
	movapd [r14], xmm0					; store first point as : Fx[N/2] | Fx[0]
	
	;the rest

	dec rcx

FFT_1D_Radix4R_ZP_ASM_N_Loop:

	add rbp, rax						;m+=t1
	add r11, rbx
	sub r12, rbx
	add r14, rbx

	movapd xmm0, [r11]					; ImFCx[k] | ReFCx[k]
	movapd xmm1, [r12]					; ImFCx[(N/2)-k] | ReFCx[(N/2)-k]
	movapd xmm2, xmm1
	addsubpd xmm2, xmm0					; SIm | DRe
	movddup xmm3, QWORD PTR [rbp]
	mulpd xmm3, xmm0					; refftc*ImFCx[k] | refftc*ReFCx[k]
	movddup xmm4, QWORD PTR [rbp+8]
	mulpd xmm4, xmm1					; refftc*ImFCx[(N/2)-k] | refftc*ReFCx[(N/2)-k]
	movddup xmm5, QWORD PTR [rbp+16]
	mulpd xmm5, xmm2					; refftc*SIm | refftc*DRe
	shufpd xmm3, xmm3, 1
	shufpd xmm4, xmm4, 1
	addsubpd xmm3, xmm4
	addpd xmm3, xmm5
	shufpd xmm3, xmm3, 1				; ImFx[k] | ReFx[k]
	movapd [r14], xmm3					; store Fx[k]

	dec rcx
	jnz FFT_1D_Radix4R_ZP_ASM_N_Loop

	mov rax, rdi						;restore rax
	ret

FFT_1D_Radix4R_ZP_ASM_N_2point:

	movapd xmm0, [r8]					; x[1] | x[0]
	movapd xmm1, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd xmm3, xmm0
	xorpd xmm2, xmm2
	shufpd xmm0, xmm2, 1
	shufpd xmm3, xmm2, 0
	movapd [r9], xmm0
	movapd [r9+16], xmm3

	ret

FFT_1D_Radix4R_ZP_ASM_N ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;ReIFFT Prepare;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;-----------------------------------------------------------------------------------------------------------------------------
;Prepares a sequence of data for half-length complex IFFT. The data is obtainable from a real sequence by FFT
;See the explanation before the IFFT_2D_Field routine
;-----------------------------------------------------------------------------------------------------------------------------
;r8 input pointer
;outputs to r8
;rdx : (N/2)
;registers changed : rax, rcx, rbp, r11, r12
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6
;rax is set to zero before returning as IFFT routine should follow this routine
ReIFFT_Prepare PROC

;do k = 0 first

	movapd xmm0, [r8]			;ReXn/2 | ReX0
	mov rbp, refftc				;m
		
	movapd xmm2, [rbp]			; 0.5 | 0.5

	mulpd xmm0, xmm2			;divide by 2
	movapd xmm1, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1			;(ReXn/2 + ReX0)/2 | (ReX0 - ReXn/2)/2
	shufpd xmm0, xmm0, 1
	movapd [r8], xmm0			;store Im Z0 | Re Z0

;k = 1 to N/4-1

	mov r11, rdx				; N/2
	mov rcx, r11
	shl rcx, 1					; rcx = N
	xor rdx, rdx
	xor rax, rax
	mov eax, MAXFFT				;MAXFFT
	idiv rcx					
	imul rax, 24				;t1 = (MAXFFT/N) * 3 * 8
	
	mov rdx, r11				;restore rdx
	shr rcx, 2					;counter : N/4

	jz ReIFFT_Prepare_return	;N = 2 : nothing more to do

	mov r11, r8					;incrementing parsing pointer
	mov r12, rdx				
	shl r12, 4					;(N/2) * 16
	add r12, r8					;decrementing parsing pointer

	add r11, 16					;
	sub r12, 16					;
	add rbp, rax				;m+=t1

	dec rcx
	jz ReIFFT_Prepare_kn4		;N = 4 : just case k = N/4 (1) left to do

ReIFFT_Prepare_loop:

	movapd xmm0, [r11]			; ImXk | ReXk
	movapd xmm1, [r12]			; ImXN/2-k | ReXN/2-k
	movapd xmm5, xmm0			; ImXk | ReXk
	movapd xmm6, xmm1			; ImXN/2-k | ReXN/2-k

	movddup xmm2, QWORD PTR [rbp]		; c1 | c1
	movddup xmm3, QWORD PTR [rbp+8]		; c2 | c2
	movddup xmm4, QWORD PTR [rbp+16]	; c3 | c3

	mulpd xmm0, xmm2			;c1 * ImXk | c1 * ReXk
	mulpd xmm1, xmm3			;c2 * ImXN/2-k | c2 * ReXN/2-k
	
	mulpd xmm3, xmm5			;c2 * ImXk | c2 * ReXk
	mulpd xmm2, xmm6			;c1 * ImXN/2-k | c1 * ReXN/2-k

	addsubpd xmm5, xmm6			;ImXk+ImXN/2-k|ReXk-ReXN/2-k
	mulpd xmm4, xmm5			;c3*(ImXk+ImXN/2-k)|c3*(ReXk-ReXN/2-k)
	
	shufpd xmm4, xmm4, 1
	addsubpd xmm0, xmm4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1			
	shufpd xmm0, xmm0, 1		;Im Zk | Re Zk

	movapd [r11], xmm0			;store Zk

	addpd xmm2, xmm4
	shufpd xmm2, xmm2, 1
	shufpd xmm3, xmm3, 1
	addsubpd xmm2, xmm3
	shufpd xmm2, xmm2, 1		;Im ZN/2-k | Re ZN/2-k

	movapd [r12], xmm2			;store ZN/2-k

	add r11, 16					;
	sub r12, 16					;
	add rbp, rax				;m+=t1

	dec rcx
	jnz ReIFFT_Prepare_loop

;k = N/4 special case

ReIFFT_Prepare_kn4:

	movapd xmm0, [r11]			; ImXN/4 | ReXN/4
	shufpd xmm0, xmm0, 1
	xorpd xmm2, xmm2
	addsubpd xmm2, xmm0
	shufpd xmm2, xmm2, 1
	movapd [r11], xmm2

ReIFFT_Prepare_return:

	xor rax, rax

	ret

ReIFFT_Prepare ENDP


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;IFFT RADIX 4 EVEN;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;---------------------------------------------------------------------------------------------------------------------------
;radix 4 complex IFFT, length even power of 2.
;Does a normal FFT using FFT_1D_Radix4_ASM_S / N routines then a divide by N and reshuffle at last stage
;---------------------------------------------------------------------------------------------------------------------------
;r8 : input (data has to be 16-byte aligned and have even number of dfp values)
;r9 : output (as above)
;r10 : calculation space
;rsi : cossin table
;rax : offset (set to zero)
;rbx : del (set to 16)
;rdx : length N (N must be even power of 2)
;registers changed : rcx, rdi, rbp, r11, r12, r14, r15
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6
IFFT_1D_Radix4_ASM_N PROC

	cmp rdx, 4										;reached 4-point fft ?
	je IFFT_1D_Radix4_ASM_N_Base

	shl rbx, 2										;del*4, N/4
	shr rdx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx									;offset restored

	;initialize loop

	mov [rsp-4], edx
	cvtsi2sd xmm6, DWORD PTR [rsp-4]
	movddup xmm6, xmm6							; N | N for ifft division

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx
	mov rax, r11									
	
	add r11, r10									;set i
	mov r12, rax									
	add r12, r9										;set o
	
	shr rcx, 2										;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx									;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	divpd xmm3, xmm6
	movapd [r12], xmm3								;store F[o1]

	mov r14, rdx
	shl r14, 4
	add r14, r12
	mov r12, r14
	sub r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	divpd xmm3, xmm6
	movapd [r12], xmm3								;store F[o2]

	sub r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	divpd xmm4, xmm6
	movapd [r12], xmm4								;store F[o3]

	sub r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	divpd xmm0, xmm6
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jz IFFT_1D_Radix4_ASM_N_Return

	xor r15, r15									

	;start loop

IFFT_1D_Radix4_ASM_N_LastStage_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	sub r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi									;save rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	divpd xmm3, xmm6
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	sub r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	divpd xmm3, xmm6
	movapd [r12], xmm3								;store F[o2]

	sub r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	divpd xmm4, xmm6
	movapd [r12], xmm4								;store F[o3]

	sub r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	divpd xmm0, xmm6
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz IFFT_1D_Radix4_ASM_N_LastStage_Loop						

	ret

IFFT_1D_Radix4_ASM_N_Base:

;4-point IFFT

	mov [rsp-4], edx
	cvtsi2sd xmm6, DWORD PTR [rsp-4]
	movddup xmm6, xmm6								; N | N for ifft division

	mov rcx, rdx
	shl rcx, 4										;N*16

	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r9										;set o
	
	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	divpd xmm3, xmm6
	movapd [r12], xmm3								;store F[o1]

	add r12, rcx
	sub r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	divpd xmm3, xmm6
	movapd [r12], xmm3								;store F[o2]

	sub r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	divpd xmm4, xmm6
	movapd [r12], xmm4								;store F[o3]

	sub r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	divpd xmm0, xmm6
	movapd [r12], xmm0								;store F[o4]

IFFT_1D_Radix4_ASM_N_Return:

	ret

IFFT_1D_Radix4_ASM_N ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 IFFT ODD;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;---------------------------------------------------------------------------------------------------------------------------
;radix 4 complex IFFT, length odd power of 2.
;Does a normal FFT using first a radix-2 step then radix-4 using FFT_1D_Radix4_ASM_S / N routines.
;Divide by N and reshuffle at last stage
;---------------------------------------------------------------------------------------------------------------------------
;r8 : input (data has to be 16-byte aligned and have even number of dfp values)
;r9 : output (as above)
;r10 : calculation space
;rsi : cossin table
;rax : offset (set to zero)
;rbx : del (set to 16)
;rdx : length N (N must be odd power of 2)
;registers changed : rcx, rdi, rbp, r11, r12, r14, r15
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6
IFFT_1D_Radix4_Odd_ASM_N PROC

	cmp rdx, 2
	je IFFT_1D_Radix4_Odd_ASM_N_Base

	shl rbx, 1										;del*2, N/2
	shr rdx, 1
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on even part (offset not changed)
	shr rbx, 1
	add rax, rbx									;offset for odd part
	shl rbx, 1
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on odd part (offset increased by del)
	shr rbx, 1
	shl rdx, 1
	sub rax, rbx									;restore offset

	;initialize loop

	mov [rsp-4], edx
	cvtsi2sd xmm6, DWORD PTR [rsp-4]
	movddup xmm6, xmm6								; N | N for ifft division

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11
	shr rcx, 1										;N/2
	mov rbp, rcx
	imul rbp, rbx									;store (N/2)*del
	
	mov r14, rax
	add r14, r9
	sub r14, rbp									;use (N/2)*del in rbp to from o2 output index

	mov rbp, rdx
	shl rbp, 4										;N * 16

	add r11, r10									;i1 - FC[i1]
	mov r12, rax
	add r12, r10									;i2 - FC[i2]
	add r12, rbx										
	mov r13, rax									
	add r13, r9										;o1 - F[o1]	
	add r14, rbp									;o2 - F[o2]
	mov r15, rsi									;m

	;k = 0 handled separately - no multiplications required

	movapd xmm0, [r11]								;ImFC[i1] | ReFC[i1]
	movapd xmm1, [r12]								;ImFC[i2] | ReFC[i2]

	movapd xmm2, xmm0
	addpd xmm0, xmm1
	divpd xmm0, xmm6
	movapd [r13], xmm0								;store F[o1]
	subpd xmm2, xmm1
	divpd xmm2, xmm6
	movapd [r14], xmm2								;store F[o1]

	dec rcx											;rcx has loop count
	add r13, rbp									;o1 index for ifft
	mov rbp, rbx
	shl rbp, 1										;del*2

	;start loop

IFFT_1D_Radix4_Odd_ASM_N_Loop:

	add r15, rdi									;m+=t
	add r11, rbp									;i1+=2*del
	add r12, rbp									;i2+=2*del
	sub r13, rbx									;o1-=del
	sub r14, rbx									;o2-=del	

	movddup xmm0, QWORD PTR [r15]					;cos	 | cos
	movddup xmm1, QWORD PTR [r15+8]					;sin	 | sin   
	movapd xmm2, [r12]								;ImFC[i2] | ReFC[i2] 
	movapd xmm3, [r11]								;ImFC[i1] | ReFC[i1]

	mulpd xmm1, xmm2								;ImFC[i2] * sin | ReFC[i2] * sin
	shufpd xmm2, xmm2, 1							;ReFC[i2] | ImFC[i2]
	mulpd xmm0, xmm2								;ReFC[i2] * cos | ImFC[i2] * cos

	addsubpd xmm0, xmm1								;tmp0 | tmp1
	shufpd xmm0, xmm0, 1							;tmp1 | tmp0

	movapd xmm1, xmm3
	addpd xmm1, xmm0								;ImF[o1] | ReF[o1]
	subpd xmm3, xmm0								;ImF[o2] | ReF[o2]

	divpd xmm1, xmm6
	movapd [r13], xmm1								;store F[o1]
	divpd xmm3, xmm6
	movapd [r14], xmm3								;store F[o2]

	dec rcx
	jnz IFFT_1D_Radix4_Odd_ASM_N_Loop						

	ret

IFFT_1D_Radix4_Odd_ASM_N_Base:

;2-point fft

	mov [rsp-4], edx
	cvtsi2sd xmm6, DWORD PTR [rsp-4]
	movddup xmm6, xmm6

	mov r11, r8
	add r11, rax									;i1 = off
	mov r12, r11
	add r12, rbx									;i2 = off + del
	mov r13, r9
	add r13, rax									;o1 = off, o2 = i2
	mov r14, r13
	add r14, rbx									;o2 = off + del

	movapd xmm0, [r11]								; Imx[i1] | Rex[i1]
	movapd xmm1, [r12]								; Imx[i2] | Rex[i2]
	movapd xmm2, xmm0

	addpd xmm0, xmm1								
	divpd xmm0, xmm6
	movapd [r13], xmm0								;store in F[o1]
	subpd xmm2, xmm1
	divpd xmm2, xmm6
	movapd [r14], xmm2								;store in F[o2]

IFFT_1D_Radix4_Odd_ASM_N_Return:

	ret
IFFT_1D_Radix4_Odd_ASM_N ENDP


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 FFT;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;---------------------------------------------------------------------------------------------------------------------------
;radix 4 FFT, complex input, complex output, length even power of 2
;FFT_1D_Radix4_ASM_N and FFT_1D_Radix4_ASM_S work in pairs
;---------------------------------------------------------------------------------------------------------------------------
;r8 : input (data has to be 16-byte aligned and have even number of dfp values)
;r9 : output (as above)
;r10 : calculation space
;rsi : cossin table
;rax : offset (set to zero)
;rbx : del (set to 16)
;rdx : length N (must be even power of 2)
;registers changed : rcx, rdi, rbp, r11, r12, r13, r14, r15
;sse registers changed : xmm0, xmm1, xmm2, xmm3, xmm4, xmm5
FFT_1D_Radix4_ASM_N PROC

	cmp rdx, 4										;reached 4-point fft ?
	je FFT_1D_Radix4_ASM_N_Base

	shl rbx, 2										;del*4, N/4
	shr rdx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_S						;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx									;offset restored
	
	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11									
	
	add r11, r10									;set i
	mov r12, rax									
	add r12, r9										;set o
	
	shr rcx, 2										;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx									;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jz FFT_1D_Radix4_ASM_N_Return

	xor r15, r15									

	;start loop

FFT_1D_Radix4_ASM_N_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	add r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi									;save rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz FFT_1D_Radix4_ASM_N_Loop						

	ret					

FFT_1D_Radix4_ASM_N_Base:

;4-point FFT

	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r9										;set o
	
	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

FFT_1D_Radix4_ASM_N_Return:

	ret

FFT_1D_Radix4_ASM_N ENDP

;swapped r9, r10
FFT_1D_Radix4_ASM_S PROC

	cmp rdx, 4										;reached 4-point fft ?
	je FFT_1D_Radix4_ASM_S_Base

	shl rbx, 2										;del*4, N/4
	shr rdx, 2
	call FFT_1D_Radix4_ASM_N						;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_N						;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_N						;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ASM_N						;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2									
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx									;offset restored

	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11									
	
	add r11, r9										;set i
	mov r12, rax									
	add r12, r10									;set o
	
	shr rcx, 2										;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx									;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jz FFT_1D_Radix4_ASM_S_Return

	xor r15, r15									

	;start loop

FFT_1D_Radix4_ASM_S_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	add r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;m+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;m+m+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz FFT_1D_Radix4_ASM_S_Loop						

	ret

FFT_1D_Radix4_ASM_S_Base:

;4-point FFT

	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r10									;set o
	
	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

FFT_1D_Radix4_ASM_S_Return:

	ret

FFT_1D_Radix4_ASM_S ENDP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 EVEN ZERO PADDED;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;normal radix4 fft with zero-padded input, even power of N (r9, r10 not swapped)
;Before calling this routine set:
;Input pointer in r8 (Re|Im 16byte-aligned doubles)
;Output pointer in r9 (Re|Im 16byte-aligned doubles)
;Pointer to calculation space in r10
;Pointer to cossin values in rsi
;offset in rax (= 0)
;del in rbx (= 1 * 16)
;N in rdx (must be even power of 2)
;Nzp *16 in r13 <= N : input is only Nzp points, Nzp to N-1 are set to zero inside this routine
;uses all registers but the following not changed : 
;r8 (input), r9 (output), r10 (calculation space), rsi (cossin),  rax (offset), rbx (del), rdx (N), r13 (Nzp * 16)
FFT_1D_Radix4_ZP_ASM_N PROC

	cmp rdx, 4											;reached 4-point fft ?
	je FFT_1D_Radix4_ZP_ASM_N_Base

	shl rbx, 2											;del*4, N/4
	shr rdx, 2
	call FFT_1D_Radix4_ZP_ASM_S							;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx										;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ZP_ASM_S							;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx										;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ZP_ASM_S							;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx										;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ZP_ASM_S							;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx										;offset restored
	
	;initialize loop

	mov r11, rax										;save offset				
	mov rcx, rdx										;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax										;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx										;restore rdx and rax
	mov rax, r11									
	
	add r11, r10										;set i
	mov r12, rax									
	add r12, r9											;set o
	
	shr rcx, 2											;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx										;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1								;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1								;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1								;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1								;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3									;cr_i2 | si_i2
	subpd xmm0, xmm3									;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5									;cr_i4 | si_i4
	subpd xmm1, xmm5									;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3									;store F[o1]

	mov r14, r12										;save o1
	add r12, rbp										;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3									;store F[o2]

	add r12, rbp										;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4									;store F[o3]

	add r12, rbp										;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0									;store F[o4]

	dec rcx
	jz FFT_1D_Radix4_ZP_ASM_N_Return

	xor r15, r15									

	;start loop

FFT_1D_Radix4_ZP_ASM_N_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	add r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi									;save rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz FFT_1D_Radix4_ZP_ASM_N_Loop						

	ret					

FFT_1D_Radix4_ZP_ASM_N_Base:

;4-point FFT

	mov r14, rax									
	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r9										;set o
	
	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	cmp r14, r13									;off>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_N_Base_Z1				;if yes then zero input
	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r14, rbx
	cmp r14, r13									;off+del>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_N_Base_Z2				;if yes then zero input
	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r14, rbx
	cmp r14, r13									;off+2*del>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_N_Base_Z3				;if yes then zero input
	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r14, rbx
	cmp r14, r13									;off+3*del>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_N_Base_Z4				;if yes then zero input
	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	jmp FFT_1D_Radix4_ZP_ASM_N_Base_Continue

FFT_1D_Radix4_ZP_ASM_N_Base_Z1:
	xorpd xmm0, xmm0
FFT_1D_Radix4_ZP_ASM_N_Base_Z2:
	xorpd xmm1, xmm1
FFT_1D_Radix4_ZP_ASM_N_Base_Z3:
	xorpd xmm3, xmm3
FFT_1D_Radix4_ZP_ASM_N_Base_Z4:
	xorpd xmm5, xmm5

FFT_1D_Radix4_ZP_ASM_N_Base_Continue:

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

FFT_1D_Radix4_ZP_ASM_N_Return:

	ret

FFT_1D_Radix4_ZP_ASM_N ENDP

;swapped radix4 fft with zero padding (r9, r10 swapped)
FFT_1D_Radix4_ZP_ASM_S PROC

	cmp rdx, 4										;reached 4-point fft ?
	je FFT_1D_Radix4_ZP_ASM_S_Base

	shl rbx, 2										;del*4, N/4
	shr rdx, 2
	call FFT_1D_Radix4_ZP_ASM_N						;Recurse FFT on part 1 (offset not changed)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ZP_ASM_N						;Recurse FFT on part 2 (offset increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ZP_ASM_N						;Recurse FFT on part 3 (offset further increased by del)
	shr rbx, 2
	add rax, rbx									;offset for next part
	shl rbx, 2
	call FFT_1D_Radix4_ZP_ASM_N						;Recurse FFT on part 4 (offset further increased by del)
	shr rbx, 2
	shl rdx, 2									
	sub rax, rbx
	sub rax, rbx
	sub rax, rbx									;offset restored

	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t1 = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11									
	
	add r11, r9										;set i
	mov r12, rax									
	add r12, r10									;set o
	
	shr rcx, 2										;N/4, loop count
	mov rbp, rcx
	imul rbp, rbx									;t2 = (N/4)*del

;case k=0 treated separately - faster

	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]

	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;ReFC[i3] | ImFC[i3]

	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;ReFC[i4] | ImFC[i4]

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jz FFT_1D_Radix4_ZP_ASM_S_Return

	xor r15, r15									

	;start loop

FFT_1D_Radix4_ZP_ASM_S_Loop:

	add r11, rbx									;i = i1 for next loop
	mov r12, r14									;restore o1
	add r12, rbx									;o = o1 for next loop
	add r15, rdi									;m+=t1 

	movapd xmm0, [r11]								
	shufpd xmm0, xmm0, 1							;ReFC[i1] | ImFC[i1]

	mov r14, rsi
	add r14, r15
	add r11, rbx									;i = i2
	movapd xmm1, [r11]								;ImFC[i2] | ReFC[i2]
	movddup xmm2, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm2, xmm1								;sin*ImFC | sin*ReFC
	shufpd xmm1, xmm1, 1							;ReFC[i2] | ImFC[i2]
	movddup xmm3, QWORD PTR [r14]						;cos | cos
	mulpd xmm1, xmm3								;cos*ReFC | cos*ImFC
	addsubpd xmm1, xmm2								; c1 | s1
	
	add r14, r15									;m+m
	add r11, rbx									;i = i3
	movapd xmm2, [r11]								;ImFC[i3] | ReFC[i3]
	movddup xmm3, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm3, xmm2								;sin*ImFC | sin*ReFC
	shufpd xmm2, xmm2, 1							;ReFC[i3] | ImFC[i3]
	movddup xmm4, QWORD PTR [r14]						;cos | cos
	mulpd xmm2, xmm4								;cos*ReFC | cos*ImFC
	addsubpd xmm2, xmm3								; c2 | s2
	
	add r14, r15									;m+m+m
	add r11, rbx									;i = i4
	movapd xmm3, [r11]								;ImFC[i4] | ReFC[i4]
	movddup xmm4, QWORD PTR [r14+8]					;sin | sin
	mulpd xmm4, xmm3								;sin*ImFC | sin*ReFC
	shufpd xmm3, xmm3, 1							;ReFC[i4] | ImFC[i4]
	movddup xmm5, QWORD PTR [r14]						;cos | cos
	mulpd xmm3, xmm5								;cos*ReFC | cos*ImFC
	addsubpd xmm3, xmm4								; c3 | s3

	movapd xmm4, xmm2
	addpd xmm4, xmm0								;cr_i2 | si_i2
	subpd xmm0, xmm2								;cr_i3 | si_i3
	movapd xmm2, xmm1
	addpd xmm2, xmm3								;cr_i4 | si_i4
	subpd xmm1, xmm3								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	mov r14, r12									;save o1
	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]

	dec rcx
	jnz FFT_1D_Radix4_ZP_ASM_S_Loop						

	ret

FFT_1D_Radix4_ZP_ASM_S_Base:

;4-point FFT

	mov r14, rax									
	mov r11, r8
	add r11, rax									;set i for x inputs
	mov r12, rax									
	add r12, r10									;set o
	
	mov rbp, rdx									
	shr rbp, 2										;N/4
	imul rbp, rbx									;t2 = (N/4)*del

	cmp r14, r13									;off>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_S_Base_Z1				;if yes then zero input
	movapd xmm0, [r11]
	shufpd xmm0, xmm0, 1							;Rex[i1] | Imx[i1]

	add r14, rbx
	cmp r14, r13									;off+del>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_S_Base_Z2				;if yes then zero input
	add r11, rbx
	movapd xmm1, [r11]
	shufpd xmm1, xmm1, 1							;Rex[i2] | Imx[i2]

	add r14, rbx
	cmp r14, r13									;off+2*del>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_S_Base_Z3				;if yes then zero input
	add r11, rbx
	movapd xmm3, [r11]
	shufpd xmm3, xmm3, 1							;Rex[i3] | Imx[i3]

	add r14, rbx
	cmp r14, r13									;off+3*del>=Nzp ?
	jge FFT_1D_Radix4_ZP_ASM_S_Base_Z4				;if yes then zero input
	add r11, rbx
	movapd xmm5, [r11]
	shufpd xmm5, xmm5, 1							;Rex[i4] | Imx[i4]

	jmp FFT_1D_Radix4_ZP_ASM_S_Base_Continue

FFT_1D_Radix4_ZP_ASM_S_Base_Z1:
	xorpd xmm0, xmm0
FFT_1D_Radix4_ZP_ASM_S_Base_Z2:
	xorpd xmm1, xmm1
FFT_1D_Radix4_ZP_ASM_S_Base_Z3:
	xorpd xmm3, xmm3
FFT_1D_Radix4_ZP_ASM_S_Base_Z4:
	xorpd xmm5, xmm5

FFT_1D_Radix4_ZP_ASM_S_Base_Continue:

	movapd xmm4, xmm0
	addpd xmm4, xmm3								;cr_i2 | si_i2
	subpd xmm0, xmm3								;cr_i3 | si_i3

	movapd xmm2, xmm1
	addpd xmm2, xmm5								;cr_i4 | si_i4
	subpd xmm1, xmm5								;ci_i4 | sr_i4

	movapd xmm3, xmm4
	addpd xmm3, xmm2
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o1]

	add r12, rbp									;o = o2
	movapd xmm3, xmm0
	shufpd xmm1, xmm1, 1
	addsubpd xmm3, xmm1
	shufpd xmm3, xmm3, 1
	movapd [r12], xmm3								;store F[o2]

	add r12, rbp									;o = o3
	subpd xmm4, xmm2
	shufpd xmm4, xmm4, 1
	movapd [r12], xmm4								;store F[o3]

	add r12, rbp									;o = o4
	shufpd xmm0, xmm0, 1
	shufpd xmm1, xmm1, 1
	addsubpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o4]
	
FFT_1D_Radix4_ZP_ASM_S_Return:

	ret

FFT_1D_Radix4_ZP_ASM_S ENDP


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;RADIX 4 ODD ZERO PADDED;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;Initial radix2 step for odd powers of 2 when using radix4 fft with zero padding(r9, r10 not swapped)
;Before calling this routine set:
;Input pointer in r8 (Re|Im 16byte-aligned doubles)
;Output pointer in r9 (Re|Im 16byte-aligned doubles)
;Pointer to calculation space in r10
;Pointer to cossin values in rsi
;offset in rax (= 0)
;del in rbx (= 1 * 16)
;N in rdx (must be odd power of 2)
;Nzp *16 in r13, Nzp <= N
;uses all registers but the following not changed : 
;r8 (input), r9 (output), r10 (calculation space), rsi (cossin),  rax (offset), rbx (del), rdx (N), r13 (Nzp * 16)
FFT_1D_Radix4_ZP_Odd_ASM_N PROC

	cmp rdx, 2
	je FFT_1D_Radix4_ZP_Odd_ASM_N_Base

	shl rbx, 1										;del*2, N/2
	shr rdx, 1
	call FFT_1D_Radix4_ZP_ASM_S						;Recurse FFT on even part (offset not changed)
	shr rbx, 1
	add rax, rbx									;offset for odd part
	shl rbx, 1
	call FFT_1D_Radix4_ZP_ASM_S						;Recurse FFT on odd part (offset increased by del)
	shr rbx, 1
	shl rdx, 1
	sub rax, rbx									;restore offset

	;initialize loop

	mov r11, rax									;save offset				
	mov rcx, rdx									;rcx(ecx) has N
	xor edx, edx
	mov eax, MAXFFT
	idiv ecx
	mov rdi, rax									;t = MAXFFT/N (*16)
	shl rdi, 4
	mov rdx, rcx									;restore rdx and rax
	mov rax, r11
	shr rcx, 1										;N/2
	mov r14, rcx
	imul r14, rbx									;store (N/2)*del
	
	add r11, r10									;i1 - FC[i1]
	mov r12, rax									
	add r12, r9										;o1 - F[o1]	
	add r14, rax									
	add r14, r9										;o2 - F[o2]
	mov r15, rsi									;m
	mov rbp, rbx
	shl rbp, 1										;del*2

	;k = 0 handled separately - no multiplications required
	movapd xmm0, [r11]								;ImFC[i1] | ReFC[i1]
	movapd xmm1, [r11+rbx]							;ImFC[i2] | ReFC[i2]

	movapd xmm2, xmm0
	addpd xmm0, xmm1
	movapd [r12], xmm0								;store F[o1]
	subpd xmm2, xmm1
	movapd [r14], xmm2								;store F[o1]

	dec rcx											;rcx has loop count
	
	;start loop

FFT_1D_Radix4_ZP_Odd_ASM_N_Loop:

	add r15, rdi									;m+=t
	add r11, rbp									;i1+=2*del
	add r12, rbx									;o1+=del
	add r14, rbx									;o2+=del	

	movddup xmm0, QWORD PTR [r15]					;cos	 | cos
	movddup xmm1, QWORD PTR [r15+8]					;sin	 | sin   
	movapd xmm2, [r11+rbx]							;ImFC[i2] | ReFC[i2] 
	movapd xmm3, [r11]								;ImFC[i1] | ReFC[i1]

	mulpd xmm1, xmm2								;ImFC[i2] * sin | ReFC[i2] * sin
	shufpd xmm2, xmm2, 1							;ReFC[i2] | ImFC[i2]
	mulpd xmm0, xmm2								;ReFC[i2] * cos | ImFC[i2] * cos

	addsubpd xmm0, xmm1								;tmp0 | tmp1
	shufpd xmm0, xmm0, 1							;tmp1 | tmp0

	movapd xmm1, xmm3
	addpd xmm1, xmm0								;ImF[o1] | ReF[o1]
	subpd xmm3, xmm0								;ImF[o2] | ReF[o2]

	movapd [r12], xmm1								;store F[o1]
	movapd [r14], xmm3								;store F[o2]

	dec rcx
	jnz FFT_1D_Radix4_ZP_Odd_ASM_N_Loop						

	ret

FFT_1D_Radix4_ZP_Odd_ASM_N_Base:

;2-point fft

	mov r11, r8
	add r11, rax									;i1 = off
	mov r12, r9
	add r12, rax									;o1 = off, o2 = i2
	mov r14, r12
	add r14, rbx									;o2 = off + del
	mov r15, rax

	cmp r15, r13
	jge FFT_1D_Radix4_ZP_Odd_ASM_N_Base_Z1
	movapd xmm0, [r11]								; Imx[i1] | Rex[i1]
	add r15, rbx
	cmp r15, r13
	jge FFT_1D_Radix4_ZP_Odd_ASM_N_Base_Z2
	movapd xmm1, [r11+rbx]							; Imx[i2] | Rex[i2]

	jmp FFT_1D_Radix4_ZP_Odd_ASM_N_Base_Continue

FFT_1D_Radix4_ZP_Odd_ASM_N_Base_Z1:
	xorpd xmm0, xmm0
FFT_1D_Radix4_ZP_Odd_ASM_N_Base_Z2:
	xorpd xmm1, xmm1

FFT_1D_Radix4_ZP_Odd_ASM_N_Base_Continue:

	movapd xmm2, xmm0
	addpd xmm0, xmm1								
	movapd [r12], xmm0								;store in F[o1]
	subpd xmm2, xmm1
	movapd [r14], xmm2								;store in F[o2]

FFT_1D_Radix4_ZP_Odd_ASM_N_Return:

	ret
FFT_1D_Radix4_ZP_Odd_ASM_N ENDP

_TEXT   ENDS
END