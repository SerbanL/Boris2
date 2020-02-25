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

//To change precision of floating point variables modify value cuBLib_FLags.h (BorisCUDALib)

//----------------------------------------------------------------- ODE EVALS

//1 : compile all ode eval methods
//2 : compile only the rk4 method (testing only, faster compilation)
//3 : custom

#define ODE_EVAL_ALL	1
#define ODE_EVAL_TEST	2
#define ODE_EVAL_CUST	3

//Set this
#define ODE_EVAL_COMPILATION	ODE_EVAL_TEST

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

#define ODE_EVAL_RK4
#define ODE_EVAL_SD
#define ODE_EVAL_RKF

#elif ODE_EVAL_COMPILATION == ODE_EVAL_CUST

#define ODE_EVAL_RK4
#define ODE_EVAL_RK23
#define ODE_EVAL_RKF
#define ODE_EVAL_ABM
#define ODE_EVAL_SD

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
#define MODULE_EXCHANGE
#define MODULE_IDMEXCHANGE
#define MODULE_ZEEMAN
#define MODULE_ANIUNI
#define MODULE_ANICUBI
#define MODULE_TRANSPORT
#define MODULE_MELASTIC
#define MODULE_ROUGHNESS
#define MODULE_SOTFIELD
#define MODULE_SDEMAG

#endif

//-----------------------------------------------------------------



