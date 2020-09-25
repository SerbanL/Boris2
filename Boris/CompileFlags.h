#pragma once

//Compilation Flags

//Most of these are used to disable various parts of the program when needed.

//This is useful for debugging and testing, and also for reducing compilation time when developing new modules.

//First develop new module in CPU code only with COMPILECUDA 0 - CUDA compilation is by far the most costly.
//Making CUDA versions of modules should be mostly a cut-and-paste job after the cpu version is done.

//----------------------------------------------------------------- OPERATING SYSTEM

#define OS_WIN	1
#define OS_LIN	2

#define OPERATING_SYSTEM	OS_LIN

//----------------------------------------------------------------- INTERFACE

//turn graphical interface on or off. If turned off then only use a basic text console.
#define GRAPHICS	0

//----------------------------------------------------------------- CUDA

//include cuda code in compilation. If set to 0 then only C++ CPU code is compiled so the executable will not need any cuda dlls.
//Much faster compilation as well, so turn this off when developing new modules.
#define COMPILECUDA	1

//To change precision of floating point variables modify value cuBLib_Flags.h (BorisCUDALib)

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

#if ATOMISTIC == 1
#define MESH_COMPILATION_ATOM_CUBIC
#endif

#elif MESH_COMPILATION == MESH_COMPILATION_CUST

#define MESH_COMPILATION_FERROMAGNETIC

#endif

//----------------------------------------------------------------- ODE EVALS

//1 : compile all ode eval methods
//2 : compile only the rk4 method (testing only, faster compilation)
//3 : custom

#define ODE_EVAL_COMPILATION_ALL	1
#define ODE_EVAL_COMPILATION_TEST	2
#define ODE_EVAL_COMPILATION_CUST	3

//Set this
#define ODE_EVAL_COMPILATION	ODE_EVAL_COMPILATION_ALL

//full
#if ODE_EVAL_COMPILATION == ODE_EVAL_COMPILATION_ALL

#define ODE_EVAL_COMPILATION_EULER
#define ODE_EVAL_COMPILATION_TEULER
#define ODE_EVAL_COMPILATION_AHEUN
#define ODE_EVAL_COMPILATION_ABM
#define ODE_EVAL_COMPILATION_RK23
#define ODE_EVAL_COMPILATION_RK4
#define ODE_EVAL_COMPILATION_RKF
#define ODE_EVAL_COMPILATION_RKCK
#define ODE_EVAL_COMPILATION_RKDP
#define ODE_EVAL_COMPILATION_SD

#elif ODE_EVAL_COMPILATION == ODE_EVAL_COMPILATION_TEST

#define ODE_EVAL_COMPILATION_EULER

#elif ODE_EVAL_COMPILATION == ODE_EVAL_COMPILATION_CUST

#define ODE_EVAL_COMPILATION_TEULER
#define ODE_EVAL_COMPILATION_RK4
#define ODE_EVAL_COMPILATION_RKF


#endif

//----------------------------------------------------------------- MODULES

//1 : compile all modules
//2 : minimal configuration
//3 : all disabled
//4 : custom

#define MODULE_COMPILATION_ALL		1
#define MODULE_COMPILATION_MIN		2
#define MODULE_COMPILATION_NONE	3
#define MODULE_COMPILATION_CUST	4

//Set this
#define MODULE_COMPILATION	MODULE_COMPILATION_ALL

//full
#if MODULE_COMPILATION == MODULE_COMPILATION_ALL

//modules enable / disable flags.
//By disabling modules you can make compilation faster, leaving only the ones you need for development / testing.
//Flags for new modules go here
#define MODULE_COMPILATION_ANIUNI
#define MODULE_COMPILATION_ANICUBI
#define MODULE_COMPILATION_DEMAG
#define MODULE_COMPILATION_DEMAG_N
#define MODULE_COMPILATION_DMEXCHANGE
#define MODULE_COMPILATION_EXCHANGE
#define MODULE_COMPILATION_IDMEXCHANGE
#define MODULE_COMPILATION_MELASTIC
#define MODULE_COMPILATION_HEAT
#define MODULE_COMPILATION_ROUGHNESS
#define MODULE_COMPILATION_SOTFIELD
#define MODULE_COMPILATION_STFIELD
#define MODULE_COMPILATION_SURFEXCHANGE
#define MODULE_COMPILATION_TRANSPORT
#define MODULE_COMPILATION_ZEEMAN
#define MODULE_COMPILATION_MOPTICAL

#if ATOMISTIC == 1
#define MODULE_COMPILATION_ATOM_DIPOLEDIPOLE
#endif

#define MODULE_COMPILATION_OERSTED
#define MODULE_COMPILATION_STRAYFIELD
#define MODULE_COMPILATION_SDEMAG

//minimal
#elif MODULE_COMPILATION == MODULE_COMPILATION_MIN

#define MODULE_COMPILATION_DEMAG
#define MODULE_COMPILATION_EXCHANGE
#define MODULE_COMPILATION_ZEEMAN

//all disabled
#elif MODULE_COMPILATION == MODULE_COMPILATION_NONE

//custom
#elif MODULE_COMPILATION == MODULE_COMPILATION_CUST

#define MODULE_COMPILATION_DEMAG
#define MODULE_COMPILATION_EXCHANGE
#define MODULE_COMPILATION_ANIUNI
#define MODULE_COMPILATION_ZEEMAN

#if ATOMISTIC == 1
#define MODULE_COMPILATION_ATOM_DIPOLEDIPOLE
#endif

#endif

//-----------------------------------------------------------------



