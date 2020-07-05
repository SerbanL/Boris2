# Boris2
Boris Computational Spintronics.

C++17 used. The codebase is contained in 781 files (.h, .cpp, .cu, .cuh, .py) and can be compiled on Windows or Linux-based OS with MSVC compiler or g++ compiler respectively.

# Download
Latest compiled version with installer found here : https://boris-spintronics.uk/download

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
<u>Windows:</u>
1. Clone the project.
2. Open the Visual Studio solution file (I use Visual Studio 2017).
3. Make sure all external dependencies are updated - see above.
4. Configure the compilation as needed - see CompileFlags.h, BorisLib_Config.h, and cuBLib_Flags.h, should be self explanatory.
5. Compile!

<u>Linux (tested on Ubuntu 20.04):</u>
1. Clone the project.
2. Make sure you have all the updates required and external dependencies:

    Get latest g++ compiler: <i>$ sudo apt install build-essential</i>
    
    Get OpenMP: <i>$ sudo apt-get install libomp-dev</i>
    
    Get CUDA: $ sudo apt install nvidia-cuda-toolkit
    
    Get SFML: <i>$ sudo apt-get install libsfml-dev</i>
    
    Get FFTW3: Instructions at http://www.fftw.org/fftw2_doc/fftw_6.html

3. Navigate to cloned directory and use the makefile to compile from source:
  
  <b>Notes</b>: before compiling on Linux-based OS you need to have the correct compilation flags in the source code (if not already set as below):

    i) Find CompileFlags.h file in Boris directory.

    i.i) Set #define OPERATING_SYSTEM	OS_LIN

    i.ii) Set #define GRAPHICS	0
  
      Graphics not currently supported on Linux, so only a basic text console is available, which is really just the graphical console text output but with all the text formatting specifiers stripped out. A portable graphical interface is on the to-do list. Simulations can still be run using Python scripts.

    ii) Find BorisLib_Config.h file in BorisLib directory.
    
    ii.i) Set #define OPERATING_SYSTEM	OS_LIN
    
    iii) If compiling CUDA code (COMPILECUDA 1 in CompileFlags.h) find cuBLib_Flags.h file in BorisCUDALib directory.
    
    iii.i) Set CUDA_ARCH to the correct value. This also has to match the nvcc compilation flag in makefile.
    
    See https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html, in particular table in section 5.2. GPU Feature List.
    
    You will need to identify the GPU architecture of your card and set the correct compilation flag for nvcc, which must match the setting in cuBLib_Flags.h. e.g. for -arch=sm_50 you need to set CUDA_ARCH to 500, etc.
    
    You can also compile in single or double floating point precision by setting the #define SINGLEPRECISION value (1 for single precision, 0 for double precision).
  
  <b>Compile:</b>
  
    $ make compile -j N
  
    (replace N with the number of logical cores on your CPU for multi-processor compilation, e.g. make compile -j 16)

    $ make install
    
    $ ./BorisLin

# Publication
A technical peer-reviewed publication on Boris to follow soon.
