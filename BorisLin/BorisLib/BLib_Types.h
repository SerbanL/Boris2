#pragma once

//Include new types used throughout BorisLib

#include "Types_VAL.h"
#include "Types_ReIm.h"
#include "Types_Rect.h"
#include "Types_Sequences.h"
#include "Types_Any.h"

//CIRCULAR INCLUSION CHECK : PASSED 

/*

#include "Types_VAL.h"
-> Funcs_Conv
	-> Introspection_base
	-> Types_Conversion
	-> Introspection_base
	-> Funcs_Math_base
-> Funcs_Strings
	-> Funcs_Vectors
-> Funcs_Aux_base

#include "Types_ReIm.h"
-> Funcs_Conv
	-> Introspection_base
	-> Types_Conversion
	-> Introspection_base
	-> Funcs_Math_base

#include "Types_Rect.h"
-> Types_VAL (see above)
-> Funcs_Math
	-> Funcs_Algorithms
	-> Funcs_Math_base.h
	-> Types_VAL (see above)

#include "Types_Sequences.h"
-> Types_VAL (see above)
-> Funcs_Math
	-> Funcs_Algorithms
	-> Funcs_Math_base.h
	-> Types_VAL (see above)
-> Funcs_Vectors
	-> Funcs_Conv
		-> Types_Conversion
			-> Introspection_base
			-> Funcs_Math_base
-> Funcs_Files
	-> Funcs_Conv
		-> Introspection_base
		-> Types_Conversion
			-> Introspection_base
			-> Funcs_Math_base
	-> Funcs_Strings
	-> Introspection
		-> Introspection_base
		-> Types_VAL
			-> Funcs_Conv
				-> Introspection_base
				-> Types_Conversion
				-> Introspection_base
				-> Funcs_Math_base
			-> Funcs_Strings
				-> Funcs_Vectors
			-> Funcs_Aux_base

#include "Types_Any.h"
	-> Types_Conversion
	-> Types_VAL
	-> Types_ReIm
	-> Types_Rect
	-> Types_Sequences
	-> Types_Info
*/