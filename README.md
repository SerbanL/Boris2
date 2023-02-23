# Boris2
Boris Computational Spintronics.

C++17 used. The codebase is currently contained in 918 files (.h, .cpp, .cu, .cuh, .py), ~164k non-trivial loc, and can be compiled on Windows or Linux-based OS with MSVC compiler or g++ compiler respectively.

# Download
Latest compiled version with installer, including source code with makefile for Linux-based OS, found here : https://boris-spintronics.uk/download

# Manual
Latest manual rolled in with installer, also found here in the Manual directory together with examples.

# External Dependencies
CUDA 9.2 or newer : https://developer.nvidia.com/cuda-92-download-archive
Python3 development version : https://www.python.org/downloads/
FFTW3 : http://www.fftw.org/download.html

# OS
The full code can be compiled on Windows 7 or Windows 10 using the MSVC compiler.
The code has also been ported to Linux (I've tested on Ubuntu 20.04) and compiled with g++, but with restrictions:

1) The graphical interface was originally written using DirectX11 so when compiling on Linux the GRAPHICS 0 flag needs to be set (see below). In the near future I plan to re-write the graphical interface in SFML.

# Building From Source
<b>Windows:</b>

1. Clone the project.
2. Open the Visual Studio solution file (I use Visual Studio 2017).
3. Make sure all external dependencies are updated - see above.
4. Configure the compilation as needed - see CompileFlags.h, BorisLib_Config.h, and cuBLib_Flags.h, should be self explanatory.
5. Compile!

<b>Linux (tested on Ubuntu 20.04):</b>

Extract the archive. On Linux-based OS the program needs to be compiled from source using the provided makefile in the extracted BorisLin directory.
Make sure you have all the required updates and dependencies:

<b>Step 0: Updates.</b>

1. Get latest g++ compiler: $ sudo apt install build-essential
2. Get OpenMP: $ sudo apt-get install libomp-dev
3. Get LibTBB: $ sudo apt install libtbb-dev
4. Get latest CUDA Toolkit (see manual for further details)
5. Get and install FFTW3: Instructions at http://www.fftw.org/fftw2_doc/fftw_6.html
6. Get Python3 development version, required for running Python scripts in embedded mode. To get Python3 development version:
$ sudo apt-get install python-dev

Open terminal and go to extracted BorisLin directory.


<b>Step 1: Configuration.</b>

<i>$ make configure (arch=xx) (sprec=0/1) (python=x.x) (cuda=x.x)</i>

Before compiling you need to set the correct CUDA architecture for your NVidia GPU.

For a list of architectures and more details see: https://en.wikipedia.org/wiki/CUDA.

Possible values for arch are:

• arch=50 is required for Maxwell architecture; translates to -arch=sm_50 in nvcc compilation.

• arch=60 is required for Pascal architecture; translates to -arch=sm_60 in nvcc compilation.

• arch=70 is required for Volta (and Turing) architecture; translates to -arch=sm_70 in nvcc compilation.

• arch=80 is required for Ampere architecture; translates to -arch=sm_80 in nvcc compilation.

• arch=90 is required for Ada (and Hopper) architecture; translates to -arch=sm_90 in nvcc compilation.

sprec sets either single precision (1) or double precision (0) for CUDA computations.

python is the Python version installed, e.g. 3.8

cuda is the CUDA Toolkit version installed, e.g. 12.0.

<b>Example: $ make configure arch=80 sprec=1 python=3.8 cuda=12.0</b>



<b>Step 2: Compilation.</b>

<i>$ make compile -j N</i>

(replace N with the number of logical cores on your CPU for multi-processor compilation)

<b>Example: $ make compile -j 16</b>



<b>Step 3: Installation.</b>

<i>$ make install</i>



<b>Step4: Run.</b>

<i>$ ./BorisLin</i>


# Publications

There are a number of articles which cover various parts of the software.

<b>General</b> (if using Boris for published works please use this as a reference)

•	S. Lepadatu, “Boris computational spintronics — High performance multi-mesh magnetic and spin transport modeling software”, Journal of Applied Physics 128, 243902 (2020)

<b>Differential equation solvers</b>

•	S. Lepadatu “Speeding Up Explicit Numerical Evaluation Methods for Micromagnetic Simulations Using Demagnetizing Field Polynomial Extrapolation” IEEE Transactions on Magnetics 58, 1 (2022)

<b>Multilayered convolution</b>

•	S. Lepadatu, “Efficient computation of demagnetizing fields for magnetic multilayers using multilayered convolution” Journal of Applied Physics 126, 103903 (2019)

<b>Parallel Monte Carlo algorithm</b>

•	S. Lepadatu, G. Mckenzie, T. Mercer, C.R. MacKinnon, P.R. Bissell, “Computation of magnetization, exchange stiffness, anisotropy, and susceptibilities in large-scale systems using GPU-accelerated atomistic parallel Monte Carlo algorithms” Journal of Magnetism and Magnetic Materials 540, 168460 (2021)

<b>Micromagnetic Monte Carlo algorithm (with demagnetizing field parallelization)</b>

•	S. Lepadatu “Micromagnetic Monte Carlo method with variable magnetization length based on the Landau–Lifshitz–Bloch equation for computation of large-scale thermodynamic equilibrium states” Journal of Applied Physics 130, 163902 (2021)

<b>Roughness effective field</b>

•	S. Lepadatu, “Effective field model of roughness in magnetic nano-structures” Journal of Applied Physics 118, 243908 (2015)

<b>Heat flow solver, LLB and 2-sublattice LLB</b>

•	S. Lepadatu, “Interaction of Magnetization and Heat Dynamics for Pulsed Domain Wall Movement with Joule Heating” Journal of Applied Physics 120, 163908 (2016)

•	S. Lepadatu “Emergence of transient domain wall skyrmions after ultrafast demagnetization” Physical Review B 102, 094402 (2020)

<b>Spin transport solver</b>

•	S. Lepadatu, “Unified treatment of spin torques using a coupled magnetisation dynamics and three-dimensional spin current solver” Scientific Reports 7, 12937 (2017)

•	S. Lepadatu, “Effect of inter-layer spin diffusion on skyrmion motion in magnetic multilayers” Scientific Reports 9, 9592 (2019)

•	C.R. MacKinnon, S. Lepadatu, T. Mercer, and P.R. Bissell “Role of an additional interfacial spin-transfer torque for current-driven skyrmion dynamics in chiral magnetic layers” Physical Review B 102, 214408 (2020)

•	C.R. MacKinnon, K. Zeissler, S. Finizio, J. Raabe, C.H. Marrows, T. Mercer, P.R. Bissell, and S. Lepadatu, “Collective skyrmion motion under the influence of an additional interfacial spin transfer torque” Scientific Reports 12, 10786 (2022)

•	S. Lepadatu and A. Dobrynin, “Self-consistent computation of spin torques and magneto-resistance in tunnel junctions and magnetic read-heads with metallic pinhole defects” Journal of Physics: Condensed Matter 35, 115801 (2023)

<b>Elastodynamics solver with thermoelastic effect and magnetostriction</b>

•	S. Lepadatu, “All-optical magento-thermo-elastic skyrmion motion” arXiv:2301.07034v2 (2023)
