#pragma once

//Conventions :

//Files ending in "_base.h" are self-contained headers, i.e. they do not include any other BorisLib headers (there may be exceptions but the exceptions are trivial). They are mainly meant for BorisLib internal use.

//Files starting in "BLib_" define a library functionality subset : include these in your project to access functionality without having to include the entire library.

//Files ending in _Windows.h or _Linux.h will only compile on given OS, and are guarded by defines. All other files have portable code. The exception is WinSocks.h which only compiles on Windows.

//To include entire BorisLib include this header

#include "BLib_Conversions.h"
#include "BLib_Types.h"
#include "BLib_Funcs.h"
#include "BLib_Network.h"
#include "BLib_prng.h"
#include "BLib_OmpReduction.h"
#include "BLib_ProgramState.h"
#include "BLib_TEquation.h"
#include "BLib_vector_lut.h"
#include "BLib_Threads.h"
#include "BLib_VEC.h"
#include "BLib_VEC_VC.h"
#include "BLib_CurveFitting.h"

