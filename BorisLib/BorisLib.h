#pragma once

//Order below to help checks for circular inclusion

//---------------------------------------------------------------------------------------
//base headers (do not include other BorisLib headers)

#include "Introspection_base.h"	//None
#include "Tuples.h"				//None
#include "WinSocks.h"			//None
#include "Funcs_Aux_base.h"		//None
#include "Funcs_Math_base.h"	//None
#include "prng.h"				//None

//---------------------------------------------------------------------------------------
//1st order : include at least one base header

//Collects various Types_.... files (apart from Types_Any.h) -> see below
#include "Types.h"				//Introspection_base.h indirectly
/* these are collected in Types.h
#include "Types_Conversion.h"	//Introspection_base.h, Funcs_Aux_base.h
#include "Types_VAL.h"			//Types_Conversion.h, Funcs_Conv.h (this is higher order but no circular inclusion arises), Funcs_Strings.h (this is higher order but no circular inclusion arises)
#include "Types_ReIm.h"			//Types_Conversion.h
#include "Types_Rect.h"			//Types_VAL.h
#include "Types_Sequences.h"	//Types_VAL.h
*/

//---------------------------------------------------------------------------------------
//2nd order : include at least one 1st order header, etc.

#include "Introspection.h"		//Types.h, Introspection_base.h
#include "Funcs_Conv.h"			//Types_Conversion.h, Introspection_base.h
#include "OmpReduction.h"		//Types.h
#include "Funcs_Algorithms.h"	//Types_Rect.h

//---------------------------------------------------------------------------------------
//3rd order

#include "Funcs_Aux.h"			//Types.h
#include "Funcs_Math.h"			//Types.h, Funcs_Algorithms.h
#include "FFT.h"				//Types.h
#include "FFT_mng.h"			//extends FFT (don't count it as higher order)
#include "FFT_shuf.h"			//extends FFT (don't count it as higher order)
#include "FFT_R2.h"				//extends FFT (don't count it as higher order)
#include "FFT_R4.h"				//extends FFT (don't count it as higher order)
#include "FFT_SR.h"				//extends FFT (don't count it as higher order)
#include "Threads.h"			//Types.h
#include "vector_lut.h"			//Types.h
#include "Funcs_Vectors.h"		//Introspection_base.h, Funcs_Conv.h
#include "Funcs_Net.h"			//Funcs_Conv.h, WinSocks.h

//---------------------------------------------------------------------------------------
//4th order

#include "Types_Any.h"			//Types.h, Funcs_Conv.h
#include "Funcs_Strings.h"		//Funcs_Vectors.h

//---------------------------------------------------------------------------------------
//5th order

#include "ProgramState.h"		//Types.h, Types_Any.h, Introspection.h, Funcs_Vectors.h
#include "VEC.h"				//Funcs_Vectors.h
#include "VEC_mng.h"			//extends VEC (don't count it as higher order)
#include "VEC_aux.h"			//extends VEC (don't count it as higher order)
#include "VEC_generate.h"		//extends VEC (don't count it as higher order)
#include "VEC_Voronoi.h"		//extends VEC (don't count it as higher order)
#include "VEC_oper.h"			//extends VEC (don't count it as higher order)
#include "VEC_trans.h"			//extends VEC (don't count it as higher order)
#include "VEC_matops.h"			//extends VEC (don't count it as higher order)
#include "VEC_MeshTransfer.h"	//extends VEC (don't count it as higher order)

//---------------------------------------------------------------------------------------
//6th order

#include "CurveFitting.h"	//VEC.h

#include "VEC_VC.h"				//VEC.h, ProgramState.h
#include "VEC_VC_mng.h"			//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_flags.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_cmbnd.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_shape.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_genshape.h"	//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_Voronoi.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_oper.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_Grad.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_Div.h"			//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_Curl.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_Del.h"			//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_diff2.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_Solve.h"		//extends VEC_VC (don't count it as higher order)
#include "VEC_VC_CGSolve.h"		//extends VEC_VC (don't count it as higher order)

#include "Funcs_Files.h"		//Funcs_Conv.h, Funcs_Vectors.h, VEC.h

//---------------------------------------------------------------------------------------
//7th order

#include "Funcs_Windows.h"		//Funcs_Files.h



