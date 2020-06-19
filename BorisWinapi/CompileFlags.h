#pragma once

//Compilation Flags

//Most of these are used to disable various parts of the program when needed.

//This is useful for debugging and testing, but most importantly for reducing compilation time when developing new modules.
//A full compilation now takes up to half an hour and rising. A bare bones compilation just for a new module can be done in under a minute.

//First develop new module in CPU code only with COMPILECUDA 0 - CUDA compilation is by far the most costly as it cannot be done using multiple cores.
//Making CUDA versions of modules should be mostly a cut-and-paste job after the cpu version is done.

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
//Much faster compilation as well, so turn this off when developing new modules.
#define COMPILECUDA	1

//To change precision of floating point variables modify value cuBLib_FLags.h (BorisCUDALib)

//----------------------------------------------------------------- MULTISCALE

//set to 0 to disable all atomistic computations
#define ATOMISTIC	1

//----------------------------------------------------------------- MESHES

#define MESH_COMPILATION_ALL	1
#define MESH_COMPILATION_TEST	2
#define MESH_COMPILATION_CUST	3

//Set this
#define MESH_COMPILATION	MESH_COMPILATION_ALL

//full
#if MESH_COMPILATION == MESH_COMPILATION_ALL

#define MESH_COMPILATION_FERROMAGNETIC
#define MESH_COMPILATION_ANTIFERROMAGNETIC
#define MESH_COMPILATION_DIAMAGNETIC
#define MESH_COMPILATION_DIPOLE
#define MESH_COMPILATION_METAL
#define MESH_COMPILATION_INSULATOR

#if ATOMISTIC == 1
#define MESH_COMPILATION_ATOM_CUBIC
#endif

#elif MESH_COMPILATION == MESH_COMPILATION_TEST

#define MESH_COMPILATION_FERROMAGNETIC

#if ATOMISTIC == 1
#define MESH_COMPILATION_ATOM_CUBIC
#endif

#elif MESH_COMPILATION == MESH_COMPILATION_CUST

#define MESH_COMPILATION_FERROMAGNETIC
#define MESH_COMPILATION_ANTIFERROMAGNETIC
#define MESH_COMPILATION_INSULATOR
#define MESH_COMPILATION_METAL

#endif

//----------------------------------------------------------------- ODE EVALS

//1 : compile all ode eval methods
//2 : compile only the rk4 method (testing only, faster compilation)
//3 : custom

#define ODE_EVAL_ALL	1
#define ODE_EVAL_TEST	2
#define ODE_EVAL_CUST	3

//Set this
#define ODE_EVAL_COMPILATION	ODE_EVAL_ALL

//full
#if ODE_EVAL_COMPILATION == ODE_EVAL_ALL

#define ODE_EVAL_EULER
#define ODE_EVAL_TEULER
#define ODE_EVAL_AHEUN
#define ODE_EVAL_ABM
#define ODE_EVAL_RK23
#define ODE_EVAL_RK4
#define ODE_EVAL_RKF
#define ODE_EVAL_RKCK
#define ODE_EVAL_RKDP
#define ODE_EVAL_SD

#elif ODE_EVAL_COMPILATION == ODE_EVAL_TEST

#define ODE_EVAL_EULER

#elif ODE_EVAL_COMPILATION == ODE_EVAL_CUST

#define ODE_EVAL_EULER
#define ODE_EVAL_TEULER
#define ODE_EVAL_AHEUN
#define ODE_EVAL_ABM
//#define ODE_EVAL_RK23
#define ODE_EVAL_RK4
#define ODE_EVAL_RKF
#define ODE_EVAL_RKCK
//#define ODE_EVAL_RKDP
//#define ODE_EVAL_SD

#endif

//----------------------------------------------------------------- MODULES

//1 : compile all modules
//2 : minimal configuration
//3 : all disabled
//4 : custom

#define COMPILE_MODULES_ALL		1
#define COMPILE_MODULES_MIN		2
#define COMPILE_MODULES_NONE	3
#define COMPILE_MODULES_CUST	4

//Set this
#define MODULE_COMPILATION	COMPILE_MODULES_ALL

//full
#if MODULE_COMPILATION == COMPILE_MODULES_ALL

//modules enable / disable flags.
//By disabling modules you can make compilation faster, leaving only the ones you need for development / testing.
//Flags for new modules go here
#define MODULE_ANIUNI
#define MODULE_ANICUBI
#define MODULE_DEMAG
#define MODULE_DEMAG_N
#define MODULE_DMEXCHANGE
#define MODULE_EXCHANGE
#define MODULE_IDMEXCHANGE
#define MODULE_MELASTIC
#define MODULE_HEAT
#define MODULE_ROUGHNESS
#define MODULE_SOTFIELD
#define MODULE_SURFEXCHANGE
#define MODULE_TRANSPORT
#define MODULE_ZEEMAN
#define MODULE_MOPTICAL

#if ATOMISTIC == 1
#define MODULE_ATOM_DIPOLEDIPOLE
#endif

#define MODULE_OERSTED
#define MODULE_STRAYFIELD
#define MODULE_SDEMAG

//minimal
#elif MODULE_COMPILATION == COMPILE_MODULES_MIN

#define MODULE_DEMAG
#define MODULE_EXCHANGE
#define MODULE_ZEEMAN

//all disabled
#elif MODULE_COMPILATION == COMPILE_MODULES_NONE

//custom
#elif MODULE_COMPILATION == COMPILE_MODULES_CUST

#define MODULE_DEMAG
#define MODULE_DEMAG_N
#define MODULE_EXCHANGE
#define MODULE_DMEXCHANGE
#define MODULE_IDMEXCHANGE
#define MODULE_ANIUNI
#define MODULE_ANICUBI
#define MODULE_ZEEMAN
#define MODULE_HEAT
#define MODULE_MOPTICAL
#define MODULE_SDEMAG

#if ATOMISTIC == 1
#define MODULE_ATOM_DIPOLEDIPOLE
#endif

#endif

//-----------------------------------------------------------------



