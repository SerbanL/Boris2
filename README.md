# Boris2
Boris Computational Spintronics.

C++17 used. The codebase is contained in 821 files (.h, .cpp, .cu, .cuh, .py) and can be compiled on Windows or Linux-based OS with MSVC compiler or g++ compiler respectively.

# Download
Latest compiled version with installer, including source code with makefile for Linux-based OS, found here : https://boris-spintronics.uk/download

# Manual
Latest manual rolled in with installer, also found here in the Manual directory together with examples.

# External Dependencies
CUDA 9.2 or newer : https://developer.nvidia.com/cuda-92-download-archive

FFTW3 : http://www.fftw.org/download.html

SFML : https://www.sfml-dev.org/download.php

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

Make sure you have all the required updates and dependencies:

Updates:

1.	Get latest g++ compiler: $ sudo apt install build-essential
2.	Get OpenMP: $ sudo apt-get install libomp-dev
3.  Get LibTBB: $ sudo apt-get install libtbb-dev
3.	Get CUDA: $ sudo apt install nvidia-cuda-toolkit
4.	Get SFML: $ sudo apt-get install libsfml-dev
5.	Get FFTW3: Instructions at http://www.fftw.org/fftw2_doc/fftw_6.html

Open terminal and go to extracted BorisLin directory.

<b>Step 1: Configuration.</b>

$ make configure (arch=xx) (sprec=0/1)

Before compiling you need to set the correct CUDA architecture for your NVidia GPU.

For a list of architectures and more details see: https://en.wikipedia.org/wiki/CUDA.

Possible values for arch are:

<li>arch=50 is required for Maxwell architecture; translates to                              -arch=sm_50 in nvcc compilation.</li>
<li>arch=60 is required for Pascal architecture; translates to                                 -arch=sm_60 in nvcc compilation.</li>
<li>arch=70 is required for Volta (and Turing) architecture; translates to                 -arch=sm_70 in nvcc compilation.</li>


Example: $ make configure arch=70

If arch is not specified a default value of 50 is used.
You can also compile CUDA code with single or double precision floating point. The default value, if not specified, is sprec=1 (single precision – recommended for most users). If you have a GPU capable of handling double precision floating point efficiently you can configure with sprec=0.

<b>Step 2: Compilation.</b>

$ make compile -j N

(replace N with the number of logical cores on your CPU for multi-processor compilation, e.g. $ make compile -j 16)

<b>Step 3: Installation.</b>

$ make install

<b>Run</b>

$ ./BorisLin

# Publication

A technical paper on Boris has been published: 
S. Lepadatu, "Boris computational spintronics — High performance multi-mesh magnetic and spin transport modeling software", J. Appl. Phys. 128, 243902 (2020).

If using Boris for published works please use this as a reference.
