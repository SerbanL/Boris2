#pragma once

namespace EqComp {

	enum FUNC_ {

		//BASIS FUNCTIONS

		//user variable (with and without param multiplication)
		FUNC_BVAR, FUNC_BVAR_PMUL,

		//user constant
		FUNC_CONST,

		//UNARY FUNCTIONS

		//sin(...), cos(...), tan(...)
		FUNC_SIN, FUNC_SIN_PMUL, FUNC_SINC, FUNC_SINC_PMUL, FUNC_COS, FUNC_COS_PMUL, FUNC_TAN, FUNC_TAN_PMUL,
		//sinh(...), cosh(...), tanh(...)
		FUNC_SINH, FUNC_SINH_PMUL, FUNC_COSH, FUNC_COSH_PMUL, FUNC_TANH, FUNC_TANH_PMUL,
		//sqrt(...), exp(...)
		FUNC_SQRT, FUNC_SQRT_PMUL, FUNC_EXP, FUNC_EXP_PMUL,

		//asin(...), acos(...), atan(...)
		FUNC_ASIN, FUNC_ASIN_PMUL, FUNC_ACOS, FUNC_ACOS_PMUL, FUNC_ATAN, FUNC_ATAN_PMUL,
		//asinh(...), acosh(...), atanh(...)
		FUNC_ASINH, FUNC_ASINH_PMUL, FUNC_ACOSH, FUNC_ACOSH_PMUL, FUNC_ATANH, FUNC_ATANH_PMUL,
		//log(...), ln(...) 
		FUNC_LOG, FUNC_LOG_PMUL, FUNC_LN, FUNC_LN_PMUL,

		//modulus
		FUNC_ABS, FUNC_ABS_PMUL,
		//step function: step(t) = 0 for t < 0, = 1 to t>= 0
		FUNC_STEP, FUNC_STEP_PMUL,
		//square and triangular waves (period 2*PI, range -1 to 1)
		FUNC_SWAV, FUNC_SWAV_PMUL, FUNC_TWAV, FUNC_TWAV_PMUL,

		//power functions where either the exponent or base is a constant (initially create FUNC_POW, then optimize later to one of these if possible)
		FUNC_POWER_EXPCONST, FUNC_POWER_EXPCONST_PMUL, FUNC_POWER_BASECONST, FUNC_POWER_BASECONST_PMUL,

		//Definite sum function : not actually evaluated, but used to expand the input std::string
		FUNC_SUM,

		//SPECIAL FUNCTIONS

		//Curie-Weiss law
		FUNC_CURIEWEISS,
		FUNC_LONGRELSUS,

		//Curie-Weiss law for 2-sublattice model
		FUNC_CURIEWEISS1, FUNC_CURIEWEISS2,
		FUNC_LONGRELSUS1, FUNC_LONGRELSUS2,
		FUNC_ALPHA1, FUNC_ALPHA2,

		//BINARY OPERATORS (order here is deliberate : lower position in enum means higher order of operation precedence)

		FUNC_POW, FUNC_POW_PMUL, FUNC_DIV, FUNC_DIV_PMUL, FUNC_MUL, FUNC_MUL_PMUL, FUNC_SUB, FUNC_SUB_PMUL, FUNC_ADD, FUNC_ADD_PMUL,

		FUNC_NUMTYPES
	};

	//struct specifying an equation component, e.g. a function, an operator, a variable, a constant.
	struct FSPEC {

		FUNC_ type;

		//user variable index in a parameter pack
		int varlevel = 0;

		//parameter for this copmonent
		double param = 1.0;

		//for FUNC_POWER_EXPCONST or FUNC_POWER_BASECONST functions this is the base or exponent value.
		//cannot set it in param since that could be used with FUNC_POWER_EXPCONST_PMUL or FUNC_POWER_BASECONST_PMUL
		double base_or_exponent = 0;

		//multiply with param?
		bool do_pmul = false;

		//indexes used by binary operators : in a vector of FSPEC objects these are indexes in that vector where we find other FPSEC objects on which the operator is applied.
		int bin_idx1 = 0;
		int bin_idx2 = 0;

		//comparison operators used for sorting
		bool operator<(const FSPEC& rhs) const { return type < rhs.type; }
		bool operator<=(const FSPEC& rhs) const { return type <= rhs.type; }
		bool operator>(const FSPEC& rhs) const { return type > rhs.type; }
		bool operator>=(const FSPEC& rhs) const { return type >= rhs.type; }
		bool operator==(const FSPEC& rhs) const { return type == rhs.type; }

		bool is_binary_operator(void) const
		{
			if (type == FUNC_ADD || type == FUNC_ADD_PMUL ||
				type == FUNC_SUB || type == FUNC_SUB_PMUL ||
				type == FUNC_MUL || type == FUNC_MUL_PMUL ||
				type == FUNC_POW || type == FUNC_POW_PMUL ||
				type == FUNC_DIV || type == FUNC_DIV_PMUL) return true;

			else return false;
		}

		bool is_basis_function(void) const
		{
			if (type == FUNC_BVAR || type == FUNC_BVAR_PMUL ||
				type == FUNC_CONST) return true;

			else return false;
		}

		bool is_constant(void) const
		{
			return (type == FUNC_CONST);
		}

		bool is_unary_function(void) const
		{
			//unary if neither binary nor basis : easier this way since there are a lot of unary functions
			return (!is_binary_operator() && !is_basis_function());
		}

		//for a version without parameter make into version with parameter multiplication (by default parameter multiplication is not set)
		void make_pmul(double param_)
		{
			if (!do_pmul) {

				param = param_;
				do_pmul = true;

				//adding 1 changes a FUNC_ into a FUNC_ with PMUL (must be applicable)
				type = (FUNC_)(type + 1);
			}
			//if pmul already set then further multiply with currently set param
			else param *= param_;
		}

		FSPEC(FUNC_ type_) :
			type(type_)
		{}

		FSPEC(FUNC_ type_, int varlevel_) :
			type(type_), varlevel(varlevel_)
		{
		}

		FSPEC(FUNC_ type_, double param_) :
			type(type_), param(param_)
		{
		}
	};
}
