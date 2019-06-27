# Boris2

Boris Computational Spintronics.

C++14 used.

# Download

Latest compiled version with installer found here : https://boris-spintronics.uk/download

# Manual

Latest manual rolled in with installer. The older v2.2 manual is found here: https://www.researchgate.net/publication/331715880_Boris_Computational_Spintronics_User_manual_version_220

# External Dependencies

CUDA 9.2 or newer : https://developer.nvidia.com/cuda-92-download-archive

FFTW3 : http://www.fftw.org/download.html

If compiling with GRAPHICS 0 flag, then you also need the latest SFML version : https://www.sfml-dev.org/download.php
With GRAPHICS 1 flag currently using DirectX 11, but plan to switching to SFML entirely in the near future for portability.

# OS

Currently only Windows 7 and Windows 10 supported. Porting to Linux is planned (there are a number of Windows-specific functions which need to be replaced; graphics code will need to be rewritten in SFML, but for now GRAPHICS 0 compilation flag can be set, enabling only a basic text console but otherwise with full functionality)

# Contributions

Contributions are welcome. 

The most straightforward type of contribution is in the form of new computational modules. A simple procedure can be followed to add self-contained modules. Similar procedures are in place for adding new types of computational meshes and material parameters. Documentation on these will be uploaded soon.