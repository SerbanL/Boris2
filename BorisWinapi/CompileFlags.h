#pragma once

//----------------------------------------------------------------- OPERATING SYSTEM

//These flags not currently used - will port to Linux eventually
#define OS_WIN	1
#define OS_LIN	2

#define OPERATING_SYSTEM	OS_WIN

//----------------------------------------------------------------- INTERFACE

//turn graphical interface on or off. If turned off then only use a basic text console.
#define GRAPHICS	1

//----------------------------------------------------------------- CUDA

//include cuda code in compilation. If set to 0 then only C++ CPU code is compiled so the executable will not need any cuda dlls.
#define COMPILECUDA	1

//compile with cuda single precision (float types) or double precision (double types). Set SINGLEPRECISION to 1 for single, otherwise (0) for double.
#define SINGLEPRECISION 1

//----------------------------------------------------------------- MODULES

//1 : compile all modules
//2 : minimal configuration
//3 : all disabled
//4 : custom

#define COMPILE_ALL		1
#define COMPILE_MIN		2
#define COMPILE_NONE	3
#define COMPILE_CUST	4

#define MODULE_COMPILATION	COMPILE_ALL

//full
#if MODULE_COMPILATION == COMPILE_ALL

//modules enable / disable flags.
//By disabling modules you can make compilation faster, leaving only the ones you need for development / testing.
//Flags for new modules go here
#define MODULE_ANIUNI
#define MODULE_ANICUBI
#define MODULE_DEMAG
#define MODULE_DEMAG_N
#define MODULE_DMEXCHANGE
#define MODULE_EXCHANGE
#define MODULE_HEAT
#define MODULE_IDMEXCHANGE
#define MODULE_ROUGHNESS
#define MODULE_SOTFIELD
#define MODULE_SURFEXCHANGE
#define MODULE_TRANSPORT
#define MODULE_ZEEMAN

#define MODULE_OERSTED
#define MODULE_STRAYFIELD
#define MODULE_SDEMAG

//minimal
#elif MODULE_COMPILATION == COMPILE_MIN

#define MODULE_DEMAG
#define MODULE_EXCHANGE
#define MODULE_ZEEMAN

//all disabled
#elif MODULE_COMPILATION == COMPILE_NONE

//custom
#elif MODULE_COMPILATION == COMPILE_CUST

#define MODULE_DEMAG
#define MODULE_EXCHANGE
#define MODULE_ZEEMAN
#define MODULE_ANIUNI
#define MODULE_IDMEXCHANGE

#define MODULE_SDEMAG

#endif

//-----------------------------------------------------------------



