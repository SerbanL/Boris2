#pragma once

#include "TEquation_FSPEC.h"

#include "Funcs_Math.h"
#include "Obj_Math_Special.h"

#include <memory>

namespace EqComp {

	template <typename ... BVarType>
	class Function {

	private:

		std::function<double(const Function&, BVarType...)> func = &Function::F_const;

		//further functions for unary functions and binary operators
		Function* pFunc1 = nullptr;
		Function* pFunc2 = nullptr;

		//for functions of an equation variable (user variable) we need to know which one it is out of a variadic list; count from 0 up.
		//used by F_bvar to return the correct value from a template parameter pack
		int varlevel = 0;

		//some functions need a parameter
		double param = 1.0;

		//for FUNC_POWER_EXPCONST or FUNC_POWER_BASECONST functions this is the base or exponent value.
		//cannot set it in param since that could be used with FUNC_POWER_EXPCONST_PMUL or FUNC_POWER_BASECONST_PMUL
		double base_or_exponent = 0;

		//Special Functions

		std::shared_ptr<Funcs_Special> pCurieWeiss = nullptr;
		std::shared_ptr<Funcs_Special> pLongRelSus = nullptr;
		std::shared_ptr<Funcs_Special> pCurieWeiss1 = nullptr;
		std::shared_ptr<Funcs_Special> pCurieWeiss2 = nullptr;
		std::shared_ptr<Funcs_Special> pLongRelSus1 = nullptr;
		std::shared_ptr<Funcs_Special> pLongRelSus2 = nullptr;
		std::shared_ptr<Funcs_Special> pAlpha1 = nullptr;
		std::shared_ptr<Funcs_Special> pAlpha2 = nullptr;

	private:

		//BASIS FUNCTIONS

		double __F_bvar(int varlevelsearch, double bvar) const
		{
			if (varlevelsearch == varlevel) {

				return bvar;
			}
			else return 0.0;
		}

		template <typename ... __BVarType>
		double __F_bvar(int varlevelsearch, double bvar, __BVarType... bvars) const
		{
			if (varlevelsearch == varlevel) {

				return bvar;
			}
			else {

				return __F_bvar(++varlevelsearch, bvars...);
			}
		}

		//variable index 0
		template <typename ... __BVarType>
		double get_0(const double& bvar0, __BVarType... bvars) const { return bvar0; }
		double get_0(const double& bvar0) const { return bvar0; }
		double get_0(...) const { return 0.0; }

		//variable index 1
		template <typename ... __BVarType>
		double get_1(const double& bvar0, const double& bvar1, __BVarType... bvars) const { return bvar1; }
		double get_1(const double& bvar0, const double& bvar1) const { return bvar1; }
		double get_1(...) const { return 0.0; }

		//variable index 2
		template <typename ... __BVarType>
		double get_2(const double& bvar0, const double& bvar1, const double& bvar2, __BVarType... bvars) const { return bvar2; }
		double get_2(const double& bvar0, const double& bvar1, const double& bvar2) const { return bvar2; }
		double get_2(...) const { return 0.0; }

		//variable index 3
		template <typename ... __BVarType>
		double get_3(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, __BVarType... bvars) const { return bvar3; }
		double get_3(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3) const { return bvar3; }
		double get_3(...) const { return 0.0; }

		//variable index 4
		template <typename ... __BVarType>
		double get_4(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, __BVarType... bvars) const { return bvar4; }
		double get_4(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4) const { return bvar4; }
		double get_4(...) const { return 0.0; }

		//variable index 5
		template <typename ... __BVarType>
		double get_5(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, __BVarType... bvars) const { return bvar5; }
		double get_5(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5) const { return bvar5; }
		double get_5(...) const { return 0.0; }

		//variable index 6
		template <typename ... __BVarType>
		double get_6(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, __BVarType... bvars) const { return bvar6; }
		double get_6(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6) const { return bvar6; }
		double get_6(...) const { return 0.0; }

		//variable index 7
		template <typename ... __BVarType>
		double get_7(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, const double& bvar7, __BVarType... bvars) const { return bvar7; }
		double get_7(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, const double& bvar7) const { return bvar7; }
		double get_7(...) const { return 0.0; }

		//variable index 8
		template <typename ... __BVarType>
		double get_8(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, const double& bvar7, const double& bvar8, __BVarType... bvars) const { return bvar8; }
		double get_8(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, const double& bvar7, const double& bvar8) const { return bvar8; }
		double get_8(...) const { return 0.0; }

		//variable index 9
		template <typename ... __BVarType>
		double get_9(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, const double& bvar7, const double& bvar8, const double& bvar9, __BVarType... bvars) const { return bvar9; }
		double get_9(const double& bvar0, const double& bvar1, const double& bvar2, const double& bvar3, const double& bvar4, const double& bvar5, const double& bvar6, const double& bvar7, const double& bvar8, const double& bvar9) const { return bvar9; }
		double get_9(...) const { return 0.0; }

		//user variable : search and return the correct one depending on varlevel value
		double F_bvar(BVarType... bvars) const
		{
			//hard-coded up to 10 variables; anything above this (really!?) you'll need function recursion (slower) to get them.
			switch (varlevel) {

			case 0:
				return get_0(bvars...);

			case 1:
				return get_1(bvars...);

			case 2:
				return get_2(bvars...);

			case 3:
				return get_3(bvars...);

			case 4:
				return get_4(bvars...);

			case 5:
				return get_5(bvars...);

			case 6:
				return get_6(bvars...);

			case 7:
				return get_7(bvars...);

			case 8:
				return get_8(bvars...);

			case 9:
				return get_9(bvars...);

			default:
				return __F_bvar(10, bvars...);
			}
		}

		double F_bvar_pmul(BVarType... bvars) const
		{
			//hard-coded up to 10 variables; anything above this (really!?) you'll need function recursion (slower) to get them.
			switch (varlevel) {

			case 0:
				return param * get_0(bvars...);

			case 1:
				return param * get_1(bvars...);

			case 2:
				return param * get_2(bvars...);

			case 3:
				return param * get_3(bvars...);

			case 4:
				return param * get_4(bvars...);

			case 5:
				return param * get_5(bvars...);

			case 6:
				return param * get_6(bvars...);

			case 7:
				return param * get_7(bvars...);

			case 8:
				return param * get_8(bvars...);

			case 9:
				return param * get_9(bvars...);

			default:
				return param * __F_bvar(10, bvars...);
			}
		}

		//constant
		double F_const(BVarType... bvars) const { return param; }

		//UNARY FUNCTIONS

		double F_sin(BVarType... bvars) const
		{
			return sin(pFunc1->func(*pFunc1, bvars...));
		}

		double F_sin_pmul(BVarType... bvars) const
		{
			return param * sin(pFunc1->func(*pFunc1, bvars...));
		}

		double F_sinc(BVarType... bvars) const
		{
			double arg = pFunc1->func(*pFunc1, bvars...);

			if (arg) return sin(arg) / arg;
			else return 1.0;
		}

		double F_sinc_pmul(BVarType... bvars) const
		{
			double arg = pFunc1->func(*pFunc1, bvars...);

			if (arg) return param * sin(arg) / arg;
			else return param;
		}

		double F_cos(BVarType... bvars) const
		{
			return cos(pFunc1->func(*pFunc1, bvars...));
		}

		double F_cos_pmul(BVarType... bvars) const
		{
			return param * cos(pFunc1->func(*pFunc1, bvars...));
		}

		double F_tan(BVarType... bvars) const
		{
			return tan(pFunc1->func(*pFunc1, bvars...));
		}

		double F_tan_pmul(BVarType... bvars) const
		{
			return param * tan(pFunc1->func(*pFunc1, bvars...));
		}

		double F_sinh(BVarType... bvars) const
		{
			return sinh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_sinh_pmul(BVarType... bvars) const
		{
			return param * sinh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_cosh(BVarType... bvars) const
		{
			return cosh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_cosh_pmul(BVarType... bvars) const
		{
			return param * cosh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_tanh(BVarType... bvars) const
		{
			return tanh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_tanh_pmul(BVarType... bvars) const
		{
			return param * tanh(pFunc1->func(*pFunc1, bvars...));
		}

		//square root
		double F_sqrt(BVarType... bvars) const
		{
			return sqrt(pFunc1->func(*pFunc1, bvars...));
		}

		double F_sqrt_pmul(BVarType... bvars) const
		{
			return param * sqrt(pFunc1->func(*pFunc1, bvars...));
		}

		double F_exp(BVarType... bvars) const
		{
			return exp(pFunc1->func(*pFunc1, bvars...));
		}

		double F_exp_pmul(BVarType... bvars) const
		{
			return param * exp(pFunc1->func(*pFunc1, bvars...));
		}

		double F_asin(BVarType... bvars) const
		{
			return asin(pFunc1->func(*pFunc1, bvars...));
		}

		double F_asin_pmul(BVarType... bvars) const
		{
			return param * asin(pFunc1->func(*pFunc1, bvars...));
		}

		double F_acos(BVarType... bvars) const
		{
			return acos(pFunc1->func(*pFunc1, bvars...));
		}

		double F_acos_pmul(BVarType... bvars) const
		{
			return param * acos(pFunc1->func(*pFunc1, bvars...));
		}

		double F_atan(BVarType... bvars) const
		{
			return atan(pFunc1->func(*pFunc1, bvars...));
		}

		double F_atan_pmul(BVarType... bvars) const
		{
			return param * atan(pFunc1->func(*pFunc1, bvars...));
		}

		double F_asinh(BVarType... bvars) const
		{
			return asinh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_asinh_pmul(BVarType... bvars) const
		{
			return param * asinh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_acosh(BVarType... bvars) const
		{
			return acosh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_acosh_pmul(BVarType... bvars) const
		{
			return param * acosh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_atanh(BVarType... bvars) const
		{
			return atanh(pFunc1->func(*pFunc1, bvars...));
		}

		double F_atanh_pmul(BVarType... bvars) const
		{
			return param * atanh(pFunc1->func(*pFunc1, bvars...));
		}

		//natual logarithm
		double F_ln(BVarType... bvars) const
		{
			return log(pFunc1->func(*pFunc1, bvars...));
		}

		double F_ln_pmul(BVarType... bvars) const
		{
			return param * log(pFunc1->func(*pFunc1, bvars...));
		}

		//logarithm base 10
		double F_log(BVarType... bvars) const
		{
			return log10(pFunc1->func(*pFunc1, bvars...));
		}

		double F_log_pmul(BVarType... bvars) const
		{
			return param * log10(pFunc1->func(*pFunc1, bvars...));
		}

		//modulus
		double F_abs(BVarType... bvars) const
		{
			return fabs(pFunc1->func(*pFunc1, bvars...));
		}

		double F_abs_pmul(BVarType... bvars) const
		{
			return param * fabs(pFunc1->func(*pFunc1, bvars...));
		}

		//step function: step(t) = 0 for t < 0, = 1 to t>= 0
		double F_step(BVarType... bvars) const
		{
			if (pFunc1->func(*pFunc1, bvars...) < 0) return 0.0;
			else return 1.0;
		}

		double F_step_pmul(BVarType... bvars) const
		{
			if (pFunc1->func(*pFunc1, bvars...) < 0) return 0.0;
			else return param;
		}

		//square wave period 2*PI range -1 to 1 : swav(0+) = +1, swav(0-) = -1
		double F_swav(BVarType... bvars) const
		{
			double value = pFunc1->func(*pFunc1, bvars...);
			return param * (-2 * (((int)floor(fabs(value) / PI) + (get_sign(value) < 0)) % 2) + 1);
		}

		double F_swav_pmul(BVarType... bvars) const
		{
			double value = pFunc1->func(*pFunc1, bvars...);
			return param * (-2 * (((int)floor(fabs(value) / PI) + (get_sign(value) < 0)) % 2) + 1);
		}

		//triangular wave period 2*PI range -1 to 1 : twav(0 to PI) = +1 to -1, twav(-PI to 0) = -1 to +1
		double F_twav(BVarType... bvars) const
		{
			double value = pFunc1->func(*pFunc1, bvars...);

			if ((int)floor(fabs(value) / PI) % 2) return (2 * fmod(fabs(value), PI) / PI - 1);
			else return (1 - 2 * fmod(fabs(value), PI) / PI);
		}

		double F_twav_pmul(BVarType... bvars) const
		{
			double value = pFunc1->func(*pFunc1, bvars...);

			if ((int)floor(fabs(value) / PI) % 2) return param * (2 * fmod(fabs(value), PI) / PI - 1);
			else return param * (1 - 2 * fmod(fabs(value), PI) / PI);
		}

		//power function : constant exponent
		double F_pow_expconst(BVarType... bvars) const
		{
			return pow(pFunc1->func(*pFunc1, bvars...), base_or_exponent);
		}

		double F_pow_expconst_pmul(BVarType... bvars) const
		{
			return param * pow(pFunc1->func(*pFunc1, bvars...), base_or_exponent);
		}

		//power function : constant base
		double F_pow_baseconst(BVarType... bvars) const
		{
			return pow(base_or_exponent, pFunc1->func(*pFunc1, bvars...));
		}

		double F_pow_baseconst_pmul(BVarType... bvars) const
		{
			return param * pow(base_or_exponent, pFunc1->func(*pFunc1, bvars...));
		}

		//SPECIAL FUNCTIONS

		//me : Curie-Weiss law
		double F_CurieWeiss(BVarType... bvars) const
		{
			if (pCurieWeiss) {

				return param * pCurieWeiss->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 1.0;
		}

		//me1 : Curie-Weiss law 2-sublattice model component 1
		double F_CurieWeiss1(BVarType... bvars) const
		{
			if (pCurieWeiss1) {

				return param * pCurieWeiss1->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 1.0;
		}

		//me2 : Curie-Weiss law 2-sublattice model component 2
		double F_CurieWeiss2(BVarType... bvars) const
		{
			if (pCurieWeiss2) {

				return param * pCurieWeiss2->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 1.0;
		}

		//chi : longitudinal relative susceptibility
		double F_LongRelSus(BVarType... bvars) const
		{
			if (pLongRelSus) {

				return param * pLongRelSus->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 0.0;
		}

		//chi1 : longitudinal relative susceptibility, 2-sublattice model component 1
		double F_LongRelSus1(BVarType... bvars) const
		{
			if (pLongRelSus1) {

				return param * pLongRelSus1->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 0.0;
		}

		//chi2 : longitudinal relative susceptibility, 2-sublattice model component 2
		double F_LongRelSus2(BVarType... bvars) const
		{
			if (pLongRelSus2) {

				return param * pLongRelSus2->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 0.0;
		}

		//alpha1 : damping scaling in 2-sublattice model
		double F_Alpha1(BVarType... bvars) const
		{
			if (pAlpha1) {

				return param * pAlpha1->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 1.0;
		}

		//alpha2 : damping scaling in 2-sublattice model
		double F_Alpha2(BVarType... bvars) const
		{
			if (pAlpha2) {

				return param * pAlpha2->evaluate(pFunc1->func(*pFunc1, bvars...));
			}
			else return 1.0;
		}

		//BINARY FUNCTIONS
		double F_add(BVarType... bvars) const
		{
			return pFunc1->func(*pFunc1, bvars...) + pFunc2->func(*pFunc2, bvars...);
		}

		double F_add_pmul(BVarType... bvars) const
		{
			return param * (pFunc1->func(*pFunc1, bvars...) + pFunc2->func(*pFunc2, bvars...));
		}

		double F_sub(BVarType... bvars) const
		{
			return pFunc1->func(*pFunc1, bvars...) - pFunc2->func(*pFunc2, bvars...);
		}

		double F_sub_pmul(BVarType... bvars) const
		{
			return param * (pFunc1->func(*pFunc1, bvars...) - pFunc2->func(*pFunc2, bvars...));
		}

		double F_mul(BVarType... bvars) const
		{
			return pFunc1->func(*pFunc1, bvars...) * pFunc2->func(*pFunc2, bvars...);
		}

		double F_mul_pmul(BVarType... bvars) const
		{
			return param * pFunc1->func(*pFunc1, bvars...) * pFunc2->func(*pFunc2, bvars...);
		}

		double F_div(BVarType... bvars) const
		{
			return pFunc1->func(*pFunc1, bvars...) / pFunc2->func(*pFunc2, bvars...);
		}

		double F_div_pmul(BVarType... bvars) const
		{
			return param * (pFunc1->func(*pFunc1, bvars...) / pFunc2->func(*pFunc2, bvars...));
		}

		double F_pow(BVarType... bvars) const
		{
			return pow(pFunc1->func(*pFunc1, bvars...), pFunc2->func(*pFunc2, bvars...));
		}

		double F_pow_pmul(BVarType... bvars) const
		{
			return param * pow(pFunc1->func(*pFunc1, bvars...), pFunc2->func(*pFunc2, bvars...));
		}

	public:

		//BASIS FUNCTIONS
		Function(FSPEC fspec)
		{
			param = fspec.param;
			varlevel = fspec.varlevel;

			switch (fspec.type) {

			case FUNC_BVAR:
				func = &Function::F_bvar;
				break;

			case FUNC_BVAR_PMUL:
				func = &Function::F_bvar_pmul;
				break;

			case FUNC_CONST:
				func = &Function::F_const;
				break;
			}
		}

		//UNARY FUNCTIONS or SPECIAL FUNCTIONS
		Function(FSPEC fspec, Function& Func)
		{
			pFunc1 = &Func;

			param = fspec.param;
			base_or_exponent = fspec.base_or_exponent;

			switch (fspec.type) {

			case FUNC_SIN:
				func = &Function::F_sin;
				break;

			case FUNC_SIN_PMUL:
				func = &Function::F_sin_pmul;
				break;

			case FUNC_SINC:
				func = &Function::F_sinc;
				break;

			case FUNC_SINC_PMUL:
				func = &Function::F_sinc_pmul;
				break;

			case FUNC_COS:
				func = &Function::F_cos;
				break;

			case FUNC_COS_PMUL:
				func = &Function::F_cos_pmul;
				break;

			case FUNC_TAN:
				func = &Function::F_tan;
				break;

			case FUNC_TAN_PMUL:
				func = &Function::F_tan_pmul;
				break;

			case FUNC_SINH:
				func = &Function::F_sinh;
				break;

			case FUNC_SINH_PMUL:
				func = &Function::F_sinh_pmul;
				break;

			case FUNC_COSH:
				func = &Function::F_cosh;
				break;

			case FUNC_COSH_PMUL:
				func = &Function::F_cosh_pmul;
				break;

			case FUNC_TANH:
				func = &Function::F_tanh;
				break;

			case FUNC_TANH_PMUL:
				func = &Function::F_tanh_pmul;
				break;

			case FUNC_SQRT:
				func = &Function::F_sqrt;
				break;

			case FUNC_SQRT_PMUL:
				func = &Function::F_sqrt_pmul;
				break;

			case FUNC_EXP:
				func = &Function::F_exp;
				break;

			case FUNC_EXP_PMUL:
				func = &Function::F_exp_pmul;
				break;

			case FUNC_ASIN:
				func = &Function::F_asin;
				break;

			case FUNC_ASIN_PMUL:
				func = &Function::F_asin_pmul;
				break;

			case FUNC_ACOS:
				func = &Function::F_acos;
				break;

			case FUNC_ACOS_PMUL:
				func = &Function::F_acos_pmul;
				break;

			case FUNC_ATAN:
				func = &Function::F_atan;
				break;

			case FUNC_ATAN_PMUL:
				func = &Function::F_atan_pmul;
				break;

			case FUNC_ASINH:
				func = &Function::F_asinh;
				break;

			case FUNC_ASINH_PMUL:
				func = &Function::F_asinh_pmul;
				break;

			case FUNC_ACOSH:
				func = &Function::F_acosh;
				break;

			case FUNC_ACOSH_PMUL:
				func = &Function::F_acosh_pmul;
				break;

			case FUNC_ATANH:
				func = &Function::F_atanh;
				break;

			case FUNC_ATANH_PMUL:
				func = &Function::F_atanh_pmul;
				break;

			case FUNC_LN:
				func = &Function::F_ln;
				break;

			case FUNC_LN_PMUL:
				func = &Function::F_ln_pmul;
				break;

			case FUNC_LOG:
				func = &Function::F_log;
				break;

			case FUNC_LOG_PMUL:
				func = &Function::F_log_pmul;
				break;

			case FUNC_ABS:
				func = &Function::F_abs;
				break;

			case FUNC_ABS_PMUL:
				func = &Function::F_abs_pmul;
				break;

			case FUNC_STEP:
				func = &Function::F_step;
				break;

			case FUNC_STEP_PMUL:
				func = &Function::F_step_pmul;
				break;

			case FUNC_SWAV:
				func = &Function::F_swav;
				break;

			case FUNC_SWAV_PMUL:
				func = &Function::F_swav_pmul;
				break;

			case FUNC_TWAV:
				func = &Function::F_twav;
				break;

			case FUNC_TWAV_PMUL:
				func = &Function::F_twav_pmul;
				break;

			case FUNC_POWER_EXPCONST:
				func = &Function::F_pow_expconst;
				break;

			case FUNC_POWER_EXPCONST_PMUL:
				func = &Function::F_pow_expconst_pmul;
				break;

			case FUNC_POWER_BASECONST:
				func = &Function::F_pow_baseconst;
				break;

			case FUNC_POWER_BASECONST_PMUL:
				func = &Function::F_pow_baseconst_pmul;
				break;

			case FUNC_CURIEWEISS:
				func = &Function::F_CurieWeiss;
				break;

			case FUNC_CURIEWEISS1:
				func = &Function::F_CurieWeiss1;
				break;

			case FUNC_CURIEWEISS2:
				func = &Function::F_CurieWeiss2;
				break;

			case FUNC_LONGRELSUS:
				func = &Function::F_LongRelSus;
				break;

			case FUNC_LONGRELSUS1:
				func = &Function::F_LongRelSus1;
				break;

			case FUNC_LONGRELSUS2:
				func = &Function::F_LongRelSus2;
				break;

			case FUNC_ALPHA1:
				func = &Function::F_Alpha1;
				break;

			case FUNC_ALPHA2:
				func = &Function::F_Alpha2;
				break;
			}
		}

		//BINARY FUNCTIONS
		Function(FSPEC fspec, Function& Func1, Function& Func2)
		{
			pFunc1 = &Func1;
			pFunc2 = &Func2;

			param = fspec.param;

			switch (fspec.type) {

			case FUNC_ADD:
				func = &Function::F_add;
				break;

			case FUNC_ADD_PMUL:
				func = &Function::F_add_pmul;
				break;

			case FUNC_SUB:
				func = &Function::F_sub;
				break;

			case FUNC_SUB_PMUL:
				func = &Function::F_sub_pmul;
				break;

			case FUNC_MUL:
				func = &Function::F_mul;
				break;

			case FUNC_MUL_PMUL:
				func = &Function::F_mul_pmul;
				break;

			case FUNC_DIV:
				func = &Function::F_div;
				break;

			case FUNC_DIV_PMUL:
				func = &Function::F_div_pmul;
				break;

			case FUNC_POW:
				func = &Function::F_pow;
				break;

			case FUNC_POW_PMUL:
				func = &Function::F_pow_pmul;
				break;
			}
		}

		void set_param(double param_value) { param = param_value; }
		void set_varlevel(int level) { varlevel = level; }

		double evaluate(BVarType... bvars) const
		{
			return this->func(*this, bvars...);
		}

		//Set special functions

		void Set_SpecialFunction(EqComp::FUNC_ type, std::shared_ptr<Funcs_Special> pSpecialFunc)
		{
			switch (type) {

			case EqComp::FUNC_CURIEWEISS:
				pCurieWeiss = pSpecialFunc;
				break;

			case EqComp::FUNC_CURIEWEISS1:
				pCurieWeiss1 = pSpecialFunc;
				break;

			case EqComp::FUNC_CURIEWEISS2:
				pCurieWeiss2 = pSpecialFunc;
				break;

			case EqComp::FUNC_LONGRELSUS:
				pLongRelSus = pSpecialFunc;
				break;

			case EqComp::FUNC_LONGRELSUS1:
				pLongRelSus1 = pSpecialFunc;
				break;

			case EqComp::FUNC_LONGRELSUS2:
				pLongRelSus2 = pSpecialFunc;
				break;

			case EqComp::FUNC_ALPHA1:
				pAlpha1 = pSpecialFunc;
				break;

			case EqComp::FUNC_ALPHA2:
				pAlpha2 = pSpecialFunc;
				break;
			}
		}
	};
};