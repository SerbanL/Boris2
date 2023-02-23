#pragma once

//Include various auxiliary functions

#include "Funcs_Strings.h"
#include "Funcs_Vectors.h"
#include "Funcs_Algorithms.h"
#include "Funcs_Math.h"
#include "Funcs_Files_Windows.h"
#include "Funcs_Files_Linux.h"
#include "Funcs_Files.h"
#include "Funcs_Aux_base.h"
#include "Funcs_Aux_Windows.h"
#include "Funcs_Aux_Linux.h"
#include "Funcs_Net_Windows.h"
#include "Funcs_Net_Linux.h"

//CIRCULAR INCLUSION CHECK : PASSED 

/*

#include "Funcs_Strings.h"
-> Funcs_Vectors

#include "Funcs_Vectors.h"
-> Funcs_Conv
	-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base

#include "Funcs_Algorithms.h"
-> Types_VAL
	-> Funcs_Conv
		-> Introspection_base
		-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base
	-> Funcs_Strings
		-> Funcs_Vectors
	-> Funcs_Aux_base

#include "Funcs_Math.h"
-> Funcs_Algorithms
-> Funcs_Math_base.h
-> Types_VAL
	-> Funcs_Conv
		-> Introspection_base
		-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base
	-> Funcs_Strings
		-> Funcs_Vectors
	-> Funcs_Aux_base

#include "Funcs_Files_Windows.h"
-> Funcs_Conv_Windows

#include "Funcs_Files.h"
-> Funcs_Conv
	-> Introspection_base
	-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base
-> Funcs_Strings
	-> Funcs_Vectors
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


#include "Funcs_Aux_base.h"

#include "Funcs_Aux_Windows.h"
-> Funcs_Aux_base
-> Types_VAL
	-> Funcs_Conv
		-> Introspection_base
		-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base
	-> Funcs_Strings
		-> Funcs_Vectors
	-> Funcs_Aux_base

#include "Funcs_Net_Windows.h"
-> Funcs_Conv
	-> Introspection_base
	-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base
-> Funcs_Conv_Windows
-> WinSocks
*/
