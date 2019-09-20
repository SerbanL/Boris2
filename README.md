# Boris2

Boris Computational Spintronics.

C++14 used. 
As a rough idea of codebase size, at the last count there were 129491 lines of code in 426 files. A full compilation without precompiled headers takes around 30 mins.

# Download

Latest compiled version with installer found here : https://boris-spintronics.uk/download

# Manual

Latest manual rolled in with installer, also found here in the Manual directory together with examples. 

# External Dependencies

CUDA 9.2 or newer : https://developer.nvidia.com/cuda-92-download-archive

FFTW3 : http://www.fftw.org/download.html

If compiling with GRAPHICS 0 flag, then you also need the latest SFML version : https://www.sfml-dev.org/download.php
With GRAPHICS 1 flag currently using DirectX 11, but plan to switching to SFML entirely in the near future for portability.

# OS

Currently only Windows 7 and Windows 10 supported. Porting to Linux is planned (there are a number of Windows-specific functions which need to be replaced but not too difficult; graphics code will need to be rewritten in SFML, but for now GRAPHICS 0 compilation flag can be set, enabling only a basic text console but otherwise with full functionality; other than this the rest of the code should port without problems - in theory!).

# Contributions

Contributions are welcome. 

The most straightforward type of contribution is in the form of new computational modules. A simple procedure can be followed to add self-contained modules. Similar procedures are in place for adding new types of computational meshes and material parameters. Documentation on these will be uploaded soon.

# Building From Source

Instructions will be provided soon.

# Publication

A technical peer-reviewed publication on Boris to follow soon.

# Bugs Status

NOT SOLVED:

1. Saving simulation file sometimes sets dT to zero (to a floating point error). I've only seen it happen with CUDA enabled. Very rare, no apparent cause found yet.

LIKELY SOLVED:

1. Using a python script may result in program hanging if issuing a flood of commands.

2. If heat solver diverges (e.g. due to too high a time step), and at least 1 material parameter has a temperature dependence, when in CUDA mode out of gpu memory errors can result requiring a program restart. 
I seem to have fixed it using extra checks on Temperature when getting updated parameter values, but I don't understand why this happens without the checks so the solution seems like a hack. Not happy with this!

SOLVED:

1. Drag and drop simulation file sometimes crashes program. Found bad conversion function - I'm certain that was the problem, so consider this solved but keep an eye on this for a while.

