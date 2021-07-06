#pragma once

#include "TEquation_Component.h"

#include "ProgramState.h"

//Build an equation from a user-supplied std::string, which can be efficiently evaluated at run-time.
//
//The equation can be typed in a natural format, detailed below.

//The usual arithmetic operators are used: ^, /, *, -, +

//Fundamental functions may be used with arguments between round brackets, e.g. sin(...)
//The allowed fundamental functions are:
//sin, cos, tan, sinh, cosh, tanh, sqrt, exp, asin, acos, atan, asinh, acosh, atanh, ln, log

//Any number of user variables may be defined, e.g. x, y, z, t

//User supplied constants may be defined, which must be alphanumeric, start with a letter, and cannot clash with any other reserved names.

//Numbers may be used in the equations, as integer, floating point, or scientific notation (e.g. 1.2e-0.3)

//Various reserved special constants are defined which may be used in the equation:
//pi : PI
//free space permeability : mu0
//Bohr magneton	: muB
//electron charge : ec
//reduced Planck constant : hbar
//Boltzmann constant : kB
//mu0 * |gamma_e|, gamma_e is the electron gyromagnetic ratio : gamma

//Example use 1:

/*

//define a TEquation object with 4 input variables named x, y, z, t
TEquation<double, double, double, double> equation({ "x", "y", "z", "t" });

//make the equation H * sin(2*PI*t) * exp(-x) * exp(-y) * exp(-z), where H is a user-defined constant
equation.make_from_string("H * sin(2*PI*t) * exp(-x) * exp(-y) * exp(-z)", { {"H", 1e6} })

//adjust the user constant later
equation.set_constant("H", 2e6);

//evaluate the equation for particular values of x, y, z, t (in this order)
double value = equation.evaluate(5, 10, 15, 3.5)

*/

//A vector equation may be specified with 3 components by using comma as a separator between 3 equations.

template <typename ... BVarType>
class TEquation :
	public ProgramState<TEquation<BVarType...>, std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::vector<double>>, std::tuple<>>
{

	bool error_on_create = false;

private:

	Equation_Component<BVarType...> eq_component_1;
	Equation_Component<BVarType...> eq_component_2;
	Equation_Component<BVarType...> eq_component_3;

	///////////////////////////////////////////////////////////

	//the user equation stored as a std::string
	std::string text_equation;

	///////////////////////////////////////////////////////////

	std::vector<std::pair<std::string, double>> reserved_constants =
	{
		//pi
		{"PI", 3.1415926535897932384626433833},
		//free space permeability
		{"mu0", 1.256637061435916e-6},
		//Bohr magneton	
		{"muB", 9.27400968e-24},
		//electron charge
		{"ec", 1.60217662e-19},
		//reduced Planck constant
		{"hbar", 1.054571817e-34},
		//Boltzmann constant
		{"kB", 1.3806488e-23},
		//mu0 * |gamma_e|, gamma_e is the electron gyromagnetic ratio
		{"gamma", 2.212761569e5}
	};

	///////////////////////////////////////////////////////////

	//text specifiers for unary functions
	std::vector<std::pair<std::string, EqComp::FSPEC>> unary_funcs =
	{
		{"sin(", EqComp::FSPEC(EqComp::FUNC_SIN)},
		{"sinc(", EqComp::FSPEC(EqComp::FUNC_SINC)},
		{"cos(", EqComp::FSPEC(EqComp::FUNC_COS)},
		{"tan(", EqComp::FSPEC(EqComp::FUNC_TAN)},
		{"sinh(", EqComp::FSPEC(EqComp::FUNC_SINH)},
		{"cosh(", EqComp::FSPEC(EqComp::FUNC_COSH)},
		{"tanh(", EqComp::FSPEC(EqComp::FUNC_TANH)},
		{"sqrt(", EqComp::FSPEC(EqComp::FUNC_SQRT)},
		{"exp(", EqComp::FSPEC(EqComp::FUNC_EXP)},
		{"asin(", EqComp::FSPEC(EqComp::FUNC_ASIN)},
		{"acos(", EqComp::FSPEC(EqComp::FUNC_ACOS)},
		{"atan(", EqComp::FSPEC(EqComp::FUNC_ATAN)},
		{"asinh(", EqComp::FSPEC(EqComp::FUNC_ASINH)},
		{"acosh(", EqComp::FSPEC(EqComp::FUNC_ACOSH)},
		{"atanh(", EqComp::FSPEC(EqComp::FUNC_ATANH)},
		{"ln(", EqComp::FSPEC(EqComp::FUNC_LN)},
		{"log(", EqComp::FSPEC(EqComp::FUNC_LOG)},
		{"abs(", EqComp::FSPEC(EqComp::FUNC_ABS)},
		{"sgn(", EqComp::FSPEC(EqComp::FUNC_SGN)},
		{"step(", EqComp::FSPEC(EqComp::FUNC_STEP)},
		{"swav(", EqComp::FSPEC(EqComp::FUNC_SWAV)},
		{"twav(", EqComp::FSPEC(EqComp::FUNC_TWAV)},
		{"me(", EqComp::FSPEC(EqComp::FUNC_CURIEWEISS)},
		{"chi(", EqComp::FSPEC(EqComp::FUNC_LONGRELSUS)},
		{"me1(", EqComp::FSPEC(EqComp::FUNC_CURIEWEISS1)},
		{"me2(", EqComp::FSPEC(EqComp::FUNC_CURIEWEISS2)},
		{"chi1(", EqComp::FSPEC(EqComp::FUNC_LONGRELSUS1)},
		{"chi2(", EqComp::FSPEC(EqComp::FUNC_LONGRELSUS2)},
		{"alpha1(", EqComp::FSPEC(EqComp::FUNC_ALPHA1)},
		{"alpha2(", EqComp::FSPEC(EqComp::FUNC_ALPHA2)}
	};

	///////////////////////////////////////////////////////////

	//used to expand the equation before processing it
	std::vector<std::pair<std::string, EqComp::FUNC_>> equation_expanders =
	{
		{"sum(", EqComp::FUNC_SUM}
	};

	///////////////////////////////////////////////////////////

	//text specifiers for user-defined function variables
	std::vector<std::string> variables;

	///////////////////////////////////////////////////////////

	//text specifiers for user-defined constants
	//must consist of alphanumeric characters only, must start with a letter
	//cannot clash with any other text specifiers
	std::vector<std::string> userconstants;
	std::vector<double> userconstants_values;

	///////////////////////////////////////////////////////////

	//text specifiers for binary operators
	//the order here specifies the precedence order (rank)
	std::vector<std::pair<std::string, EqComp::FSPEC>> binop =
	{
		{"^", EqComp::FSPEC(EqComp::FUNC_POW)},
		{"/", EqComp::FSPEC(EqComp::FUNC_DIV)},
		{"*", EqComp::FSPEC(EqComp::FUNC_MUL)},
		{"-", EqComp::FSPEC(EqComp::FUNC_SUB)},
		{"+", EqComp::FSPEC(EqComp::FUNC_ADD)}
	};

	///////////////////////////////////////////////////////////

	//Special Functions (set in TEquation_Function when remaking equation)

	std::shared_ptr<Funcs_Special> pCurieWeiss = nullptr;
	std::shared_ptr<Funcs_Special> pLongRelSus = nullptr;
	std::shared_ptr<Funcs_Special> pCurieWeiss1 = nullptr;
	std::shared_ptr<Funcs_Special> pCurieWeiss2 = nullptr;
	std::shared_ptr<Funcs_Special> pLongRelSus1 = nullptr;
	std::shared_ptr<Funcs_Special> pLongRelSus2 = nullptr;
	std::shared_ptr<Funcs_Special> pAlpha1 = nullptr;
	std::shared_ptr<Funcs_Special> pAlpha2 = nullptr;

private:

	/////////////////////////////////////////////////////////
	//
	// AUXILIARY

	//expand equation if any terms in equation_expanders appear; return true if expansion done, false otherwise (including if there was an equation error - this will be caught later anyway)
	bool expand_equation_string(std::string& eq_str)
	{
		for (auto entry : equation_expanders) {

			switch (entry.second) {

			case EqComp::FUNC_SUM:
			{
				std::string sum_str = entry.first;

				size_t sum_pos = eq_str.find(sum_str);
				if (sum_pos == std::string::npos) return false;

				size_t pos_sc1 = eq_str.find(";");
				if (pos_sc1 == std::string::npos) return false;

				size_t pos_sc2 = eq_str.find(";", pos_sc1 + 1);
				if (pos_sc2 == std::string::npos) return false;

				size_t pos_sc3 = eq_str.find(";", pos_sc2 + 1);
				if (pos_sc3 == std::string::npos) return false;

				std::string var = std::string("<") + eq_str.substr(sum_pos + sum_str.length(), pos_sc1 - sum_str.length() - sum_pos) + std::string(">");
				int lo = ToNum(eq_str.substr(pos_sc1 + 1, pos_sc2 - pos_sc1 - 1));
				int hi = ToNum(eq_str.substr(pos_sc2 + 1, pos_sc3 - pos_sc2 - 1));
				if (lo >= hi) return false;

				int num_open = 1;
				int idx_start = pos_sc3 + 1;

				for (int idx_str = idx_start; idx_str < eq_str.length(); idx_str++) {

					if (eq_str[idx_str] == '(') { num_open++; continue; }
					if (eq_str[idx_str] == ')') {

						num_open--;

						if (!num_open) {

							//found substring between brackets : work on this by going a level deeper to build branch
							std::string eq_substr = std::string("(") + eq_str.substr(idx_start, idx_str - idx_start) + std::string(")");

							//now expand eq_substr as a sum, by replacing occurences of <var> by numerical values
							std::string eq_expanded = std::string("(");

							for (int idx = lo; idx <= hi; idx++) {

								std::string sum_term = eq_substr;
								replaceall(sum_term, var, ToString(idx));

								if (idx != lo) eq_expanded += std::string("+") + sum_term;
								else eq_expanded += sum_term;
							}

							eq_expanded += std::string(")");

							eq_str.replace(sum_pos, idx_str - sum_pos + 1, eq_expanded);
							//replacement done : there may be more so return true
							return true;
						}
					}
				}

				//If we are here then number of brackets didn't match, so equation is wrong
				return false;
			}
			break;

			//other cases here if any needed in the future
			default:
				return false;
				break;
			}
		}

		return false;
	}

	void Set_StoredSpecialFunctions(void)
	{
		auto set = [&](EqComp::FUNC_ type, std::shared_ptr<Funcs_Special> pSpecialFunction) -> void {

			if (pSpecialFunction) {

				eq_component_1.Set_SpecialFunction(type, pSpecialFunction);
				eq_component_2.Set_SpecialFunction(type, pSpecialFunction);
				eq_component_3.Set_SpecialFunction(type, pSpecialFunction);
			}
		};

		set(EqComp::FUNC_CURIEWEISS, pCurieWeiss);
		set(EqComp::FUNC_LONGRELSUS, pLongRelSus);
		set(EqComp::FUNC_CURIEWEISS1, pCurieWeiss1);
		set(EqComp::FUNC_CURIEWEISS2, pCurieWeiss2);
		set(EqComp::FUNC_LONGRELSUS1, pLongRelSus1);
		set(EqComp::FUNC_LONGRELSUS2, pLongRelSus2);
		set(EqComp::FUNC_ALPHA1, pAlpha1);
		set(EqComp::FUNC_ALPHA2, pAlpha2);
	}

	//clear allocated Function objects
	void clear_Funcs(void)
	{
		eq_component_1.clear_Funcs();
		eq_component_2.clear_Funcs();
		eq_component_3.clear_Funcs();
	}

	bool is_unary_function(std::string text, int *pidx = nullptr)
	{
		for (int idx = 0; idx < unary_funcs.size(); idx++) {

			if (text == unary_funcs[idx].first) {

				if (pidx) *pidx = idx;
				return true;
			}
		}

		return false;
	}

	bool is_binary_operator(std::string text, int *pidx = nullptr)
	{
		for (int idx = 0; idx < binop.size(); idx++) {

			if (text == binop[idx].first) {

				if (pidx) *pidx = idx;
				return true;
			}
		}

		return false;
	}

	int count_bracket_pairs(std::string eq_str)
	{
		int num_open = 0;

		for (int idx_str = 0; idx_str < eq_str.length(); idx_str++) {

			if (eq_str[idx_str] == '(') num_open++;
			if (eq_str[idx_str] == ')') num_open--;
		}

		return num_open;
	}

	/////////////////////////////////////////////////////////
	//
	// MAKE USER CONSTANTS

	void create_constants(const std::vector<std::pair<std::string, double>>& constants)
	{
		//check all constants names
		std::vector<std::string> constants_names(constants.size());
		std::vector<double> constants_values(constants.size());

		bool success = true;

		for (int idx = 0; idx < constants_names.size(); idx++) {

			constants_names[idx] = constants[idx].first;
			constants_values[idx] = constants[idx].second;

			//the constants names have to be alphanumberic and start with a letter
			//cannot clash with any other specifiers (function names and user variables)

			//alphanumeric and must start with letter
			success &= (std::find_if(constants_names[idx].begin(), constants_names[idx].end(), [](const char& c) { return !isalnum(c); }) == constants_names[idx].end()) && constants_names[idx].length() && isalpha(constants_names[idx][0]);

			//cannot clash with named variables
			success &= (std::find_if(variables.begin(), variables.end(), [&](const std::string& entry) { return (entry == constants_names[idx]); }) == variables.end());

			//cannot clash with reserved constants names
			success &= (std::find_if(reserved_constants.begin(), reserved_constants.end(), [&](const std::pair<std::string, double>& entry) { return (entry.first == constants_names[idx]); }) == reserved_constants.end());

			//cannot clash with reserved function names
			success &= (std::find_if(unary_funcs.begin(), unary_funcs.end(), [&](const std::pair<std::string, EqComp::FSPEC>& entry) { return (entry.first.substr(0, entry.first.length() - 1) == constants_names[idx]); }) == unary_funcs.end());

			//add entry if not already present, else modify value
			if (success) {

				int idx_uc = search_vector(userconstants, constants_names[idx]);

				if (idx_uc >= 0) {

					userconstants_values[idx_uc] = constants_values[idx];
				}
				else {

					userconstants.push_back(constants_names[idx]);
					userconstants_values.push_back(constants_values[idx]);
				}
			}
		}
	}

	/////////////////////////////////////////////////////////
	//
	// CREATE fspec FROM TEXT EQUATION BRANCH BY BRANCH

	//make equation branch, either by continuing an existing branch (EqComp::FSPEC_branch) or by starting a new one. In the end add the branch to EqComp::FSPEC.
	bool make_branch(std::string eq_str, std::vector< std::vector<EqComp::FSPEC> >& fspec, std::vector<EqComp::FSPEC>& fspec_branch)
	{
		//first identify binary operators which join branches at this level, i.e. which have only matched pairs of brackets to the left
		//a binary operator defines a node which join two branches, thus if we identify binary operators we have to start new branches
		std::vector<std::string> eq_branches;
		std::vector<EqComp::FSPEC> eq_binops;

		std::function<void(std::string, int)> split_equation = [&](std::string eq_substr, int start_idx) -> void {

			//do not look at first character : if it contains an operator it will be '-' meaning a negative constant.
			for (int idx = start_idx + 1; idx < eq_substr.length(); idx++) {

				int op_idx = 0;
				if (is_binary_operator(eq_substr.substr(idx, 1), &op_idx)) {

					//found a binary operator. The exception is when '-' is used with a scientific notation number; in this case the pattern is ae-b, where a and b are digits.
					if (binop[op_idx].second.type == EqComp::FUNC_SUB) {

						if (idx >= 2 && idx < eq_substr.size() - 1 &&
							eq_substr[idx - 1] == 'e' &&
							std::string("0123456789").find(eq_substr.substr(idx - 2, 1)) != std::string::npos &&
							std::string("0123456789").find(eq_substr.substr(idx + 1, 1)) != std::string::npos) continue;
					}

					std::string eq_substr_left = eq_substr.substr(start_idx, idx - start_idx);

					//does the equation to the left of this operator contain only matched bracket pairs? i.e. is this operator a valid node at this tree level?
					if (!count_bracket_pairs(eq_substr_left)) {

						//found a node -> add it to equation branches std::string
						//also add the operator. Don't set branch joining indexes yet, we'll make these when we make the actual branches
						eq_branches.push_back(eq_substr_left);
						eq_binops.push_back(binop[op_idx].second);

						//there must be a std::string to the right of this operator and now we want to break it down further at this tree level
						return split_equation(eq_substr, idx + 1);
					}
				}
			}

			//we get here at the last step if no further operators found
			eq_branches.push_back(eq_substr.substr(start_idx));
		};

		//if possible split the equation into a number of substrings which are joined by binary operators at this level
		split_equation(eq_str, 0);

		//if we managed to split the equation down into a number of branches then must process all the branches
		if (eq_binops.size()) {

			//first make all the branches separately, then we'll join them
			for (int idx = 0; idx < eq_branches.size(); idx++) {

				//make branch
				std::vector<EqComp::FSPEC> fspec_new_branch;
				if (!make_branch(eq_branches[idx], fspec, fspec_new_branch)) return false;

				//set branch joining indexes in operators
				if (idx < eq_branches.size() - 1) eq_binops[idx].bin_idx1 = fspec.size() - 1;

				if (idx > 0) {

					eq_binops[idx - 1].bin_idx2 = fspec.size() - 1;
				}
			}

			//now start joining branches in turn in order of operator rank (this is specified by the order in binop member data)
			//1) ^ 2) / 3) * 4) - 5) +

			//parse eq_binops and eliminate one operator at a time by making a new node; recurse until none left
			//consider them in order of precedence, thus sort them first
			std::sort(eq_binops.begin(), eq_binops.end());

			for (int idx = 0; idx < eq_binops.size(); idx++) {

				fspec.push_back(std::vector<EqComp::FSPEC>{eq_binops[idx]});

				//adjust indexes in remaining eq_binops : all branch indexes which contain an index from eq_binops[idx] must be replaced with EqComp::FSPEC.size() - 1
				for (int idx_rest = idx + 1; idx_rest < eq_binops.size(); idx_rest++) {

					if (eq_binops[idx_rest].bin_idx1 == eq_binops[idx].bin_idx2) eq_binops[idx_rest].bin_idx1 = fspec.size() - 1;
					if (eq_binops[idx_rest].bin_idx2 == eq_binops[idx].bin_idx1) eq_binops[idx_rest].bin_idx2 = fspec.size() - 1;
				}
			}

			if (fspec_branch.size()) {

				fspec.push_back(fspec_branch);
			}

			return true;
		}

		//next search for:

		//1. function name with opening bracket, e.g. "sin("
		//		-> on finding this go a level deeper by looking inside the bracket pair

		//2. opening bracket "("
		//		-> on finding this go a level deeper by looking inside the bracket pair

		//3. variable name, e.g. "x", "y", "z", "t"
		//		-> on finding this finish branch and return (basis function)

		//4. user constant
		//		-> on finding this finish branch and return (basis function)

		//5. reserved constant
		//		-> on finding this finish branch and return (basis function)

		//6. number
		//		-> on finding this finish branch and return (basis function)

		//1. functions
		for (int idx_ufunc = 0; idx_ufunc < unary_funcs.size(); idx_ufunc++) {

			//find function name
			size_t idx_found = eq_str.find(unary_funcs[idx_ufunc].first);

			if (idx_found == 0) {

				//found name: insert specifier in current branch
				fspec_branch.insert(fspec_branch.begin(), unary_funcs[idx_ufunc].second);

				//find closing bracket ")" so we can get substring between brackets
				//number of open ( brackets : when this reaches zero we found last closing bracket
				int num_open = 1;
				int idx_start = unary_funcs[idx_ufunc].first.length();

				for (int idx_str = idx_start; idx_str < eq_str.length(); idx_str++) {

					if (eq_str[idx_str] == '(') { num_open++; continue; }
					if (eq_str[idx_str] == ')') {

						num_open--;

						if (!num_open) {

							//found substring between brackets : work on this by going a level deeper to build branch
							std::string eq_substr = eq_str.substr(idx_start, idx_str - idx_start);

							//go up a step in depth in this branch
							return make_branch(eq_substr, fspec, fspec_branch);
						}
					}
				}

				//If we are here then number of brackets didn't match, so equation is wrong
				return false;
			}
		}

		//2. bracket
		if (eq_str.length() && eq_str[0] == '(') {

			int num_open = 1;

			//get substring enclosed by matching bracket pair
			for (int idx = 1; idx < eq_str.length(); idx++) {

				if (eq_str[idx] == '(') num_open++;
				if (eq_str[idx] == ')') {

					num_open--;
					if (!num_open) {

						std::string eq_substr = eq_str.substr(1, idx - 1);

						//start new branch
						return make_branch(eq_substr, fspec, fspec_branch);
					}
				}
			}
		}

		//3. variables
		for (int idx_var = 0; idx_var < variables.size(); idx_var++) {

			//since a variable forms a basis function, if this is a user variable then it must be the entire eq_str we are currently working on.
			if (eq_str == variables[idx_var]) {

				//found variable name
				//insert variable in branch : this is the end of this branch so add it to tree
				fspec_branch.insert(fspec_branch.begin(), EqComp::FSPEC(EqComp::FUNC_BVAR, idx_var));
				fspec.push_back(fspec_branch);
				return true;
			}
		}

		//4. user constants
		for (int idx_uc = 0; idx_uc < userconstants.size(); idx_uc++) {

			//since a user constant forms a basis function, if this is a user constant then it must be the entire eq_str we are currently working on.
			if (eq_str == userconstants[idx_uc]) {

				//found user constant name
				//insert constant in branch : this is the end of this branch so add it to tree
				fspec_branch.insert(fspec_branch.begin(), EqComp::FSPEC(EqComp::FUNC_CONST, userconstants_values[idx_uc]));
				fspec.push_back(fspec_branch);

				return true;
			}
		}

		//5. reserved constants
		for (int idx_ruc = 0; idx_ruc < reserved_constants.size(); idx_ruc++) {

			//since a reserved constant forms a basis function, if this is a reserved constant then it must be the entire eq_str we are currently working on.
			if (eq_str == reserved_constants[idx_ruc].first) {

				//found reserved constant name
				//insert constant in branch : this is the end of this branch so add it to tree
				fspec_branch.insert(fspec_branch.begin(), EqComp::FSPEC(EqComp::FUNC_CONST, reserved_constants[idx_ruc].second));
				fspec.push_back(fspec_branch);

				return true;
			}
		}

		//6. constants

		//since a constant forms a basis function, if this is a constant then it must be the entire eq_str we are currently working on.
		if (has_numbers_only(eq_str, "")) {

			fspec_branch.insert(fspec_branch.begin(), EqComp::FSPEC(EqComp::FUNC_CONST, (double)ToNum(eq_str)));
			fspec.push_back(fspec_branch);
			return true;
		}

		//something not right
		return false;
	}

public:

	/////////////////////////////////////////////////////////
	//
	// CONSTRUCTOR

	TEquation(void) :
		ProgramState<TEquation<BVarType...>, std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::vector<double>>, std::tuple<>>
		(this, { VINFO(text_equation), VINFO(variables), VINFO(userconstants), VINFO(userconstants_values) }, {})
	{}

	//make equation object with names of user variables
	TEquation(const std::vector<std::string>& bvar_names) :
		ProgramState<TEquation<BVarType...>, std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::vector<double>>, std::tuple<>>
		(this, { VINFO(text_equation), VINFO(variables), VINFO(userconstants), VINFO(userconstants_values) }, {})
	{
		if (sizeof...(BVarType) == bvar_names.size()) {

			variables = bvar_names;
		}
		else error_on_create = true;
	}

	TEquation(const TEquation& copy_this) :
		ProgramState<TEquation<BVarType...>, std::tuple<std::string, std::vector<std::string>, std::vector<std::string>, std::vector<double>>, std::tuple<>>
		(this, { VINFO(text_equation), VINFO(variables), VINFO(userconstants), VINFO(userconstants_values) }, {})
	{
		*this = copy_this;
	}

	~TEquation() { clear_Funcs(); }

	//---------Assignment operator

	TEquation& operator=(const TEquation& copy_this)
	{
		clear_Funcs();

		text_equation = copy_this.text_equation;
		variables = copy_this.variables;
		userconstants = copy_this.userconstants;
		userconstants_values = copy_this.userconstants_values;

		pCurieWeiss = copy_this.pCurieWeiss;
		pLongRelSus = copy_this.pLongRelSus;
		pCurieWeiss1 = copy_this.pCurieWeiss1;
		pCurieWeiss2 = copy_this.pCurieWeiss2;
		pLongRelSus1 = copy_this.pLongRelSus1;
		pLongRelSus2 = copy_this.pLongRelSus2;
		pAlpha1 = copy_this.pAlpha1;
		pAlpha2 = copy_this.pAlpha2;

		remake_equation();

		return *this;
	}

	//-------------------Implement ProgramState method

	void RepairObjectState(void)
	{
		//remake the TEquation object from reloaded text equation std::string
		clear_Funcs();

		std::vector<std::pair<std::string, double>> constants_and_values(userconstants.size());

		for (int idx = 0; idx < constants_and_values.size(); idx++) {

			constants_and_values[idx] = std::pair<std::string, double>(userconstants[idx], userconstants_values[idx]);
		}

		make_from_string(text_equation, constants_and_values);
	}

	/////////////////////////////////////////////////////////
	//
	// USER CONSTANTS

	//set value of named constant
	bool set_constant(const std::string& name, double value, bool remake_equation = true)
	{
		int idx_uc = search_vector(userconstants, name);

		if (idx_uc >= 0) {

			//constant exists : update value
			userconstants_values[idx_uc] = value;
		}
		else {

			//new constant defined: add it
			create_constants({ {name, value} });
		}

		if (remake_equation) {

			//remake equation
			std::vector<std::pair<std::string, double>> constants(userconstants.size());
			for (int idx = 0; idx < constants.size(); idx++) {

				constants[idx] = std::pair<std::string, double>(userconstants[idx], userconstants_values[idx]);
			}

			return make_from_string(text_equation, constants);
		}
		else return true;
	}

	//set value of named constants
	bool set_constants(const std::vector<std::pair<std::string, double>>& constants, bool remake_equation = true)
	{
		//if constants already defined then the following call onyl sets values
		create_constants(constants);

		if (remake_equation) {

			//remake equation
			return make_from_string(text_equation, constants);
		}
		else return true;
	}

	//get value of named constant
	double get_constant(const std::string& name)
	{
		int idx_uc = search_vector(userconstants, name);

		if (idx_uc >= 0) {

			return userconstants_values[idx_uc];
		}
		else return 0.0;
	}

	/////////////////////////////////////////////////////////
	//
	// MAKE EQUATION

	//remake equation from currently stored equation text and user constants
	bool remake_equation(void)
	{
		//remake equation
		std::vector<std::pair<std::string, double>> constants(userconstants.size());
		for (int idx = 0; idx < constants.size(); idx++) {

			constants[idx] = std::pair<std::string, double>(userconstants[idx], userconstants_values[idx]);
		}

		return make_from_string(text_equation, constants);
	}

	bool make_from_string(std::string eq_str, const std::vector<std::pair<std::string, double>>& constants = std::vector<std::pair<std::string, double>>())
	{
		if (!eq_str.length()) return false;

		clear_Funcs();

		if (error_on_create) return false;

		//From equation std::string build an equation with Function objects which can be numerically evaluated.

		//The equation can be represented as a tree with all the branches converging to a single point, the evaluated value.
		//all lowest-level branches start with a basis function, i.e. either a constant or a function variable.
		//Two branches can be combined using a binary operator (e.g. arithmetic operation), thus higher-level branches can also start with a binary operator.

		//1. Translate the text equation into a EqComp::FSPEC equation specifier. This is a logical representation of the text equation.
		//2. Translate the EqComp::FSPEC specifier structure into an actual equation with allocated Function objects.

		//prepare the equation std::string:
		//remove white spaces
		eq_str = trimspaces(eq_str);

		//- signs can be used to mean multiplication by -1. Explicitly replace these by "-1*" if "-" appears a) at the start, b) after a ( bracket, c) after a binary operator.
		//do not replace if it's followed by a digit.

		//store user equation as text also
		text_equation = eq_str;

		int idx = 0;
		while (idx < eq_str.length()) {

			if (eq_str[idx] == '-') {

				if (idx < eq_str.length() - 1 && std::string("0123456789").find(eq_str.substr(idx + 1, 1)) == std::string::npos) {

					if (idx == 0 || eq_str[idx - 1] == '(' || is_binary_operator(eq_str.substr(idx - 1, 1))) {

						eq_str.replace(idx, 1, "-1*");

						idx += 3;
						continue;
					}
				}
			}

			idx++;
		}

		auto make_equation_component = [&](std::string eq_str_component, Equation_Component<BVarType...>& eq_component) -> bool {

			//each entry in EqComp::FSPEC specifies a branch as a cascade of EqComp::FSPEC objects from bottom up.

			//expand std::string if needed according to entries in equation_expanders
			while (expand_equation_string(eq_str_component));

			auto& fspec = eq_component.get_fspec();

			std::vector<EqComp::FSPEC> fspec_branch;

			//now translate text to EqComp::FSPEC
			bool success = make_branch(eq_str_component, fspec, fspec_branch);

			//finally make the equation with Function objects
			if (success) success &= eq_component.make();

			return success;
		};

		bool success = true;

		std::vector<std::string> eq_str_components = split(eq_str, ",");

		create_constants(constants);

		success &= make_equation_component(eq_str_components[0], eq_component_1);

		if (eq_str_components.size() >= 2) success &= make_equation_component(eq_str_components[1], eq_component_2);
		if (eq_str_components.size() >= 3) success &= make_equation_component(eq_str_components[2], eq_component_3);

		//if equation std::string was not correct then need to clear Functions otherwise the structure is probably unstable and will crash the program if executed.
		if (!success) clear_Funcs();

		//Set Special functions if available
		Set_StoredSpecialFunctions();

		return success;
	}

	/////////////////////////////////////////////////////////
	//
	// DESTROY EQUATION

	void clear(void)
	{
		clear_Funcs();
		text_equation = "";
		userconstants.clear();
		userconstants_values.clear();
	}

	/////////////////////////////////////////////////////////
	//
	// EQUATION TYPE CHECKING

	bool is_set(void) const { return eq_component_1.is_set(); }
	bool is_set_dual(void) const { return eq_component_1.is_set() && eq_component_2.is_set(); }
	bool is_set_vector(void) const { return eq_component_1.is_set() && eq_component_2.is_set() && eq_component_3.is_set(); }

	/////////////////////////////////////////////////////////
	//
	// SHOW VARIOUS PROPERTIES

	std::string show_equation(void) const { return text_equation; }

	std::string show_functions(void)
	{
		std::vector<std::string> functions_list(unary_funcs.size());

		for (int idx = 0; idx < unary_funcs.size(); idx++) {

			functions_list[idx] = unary_funcs[idx].first.substr(0, unary_funcs[idx].first.length() - 1);
		}

		return combine(functions_list, ", ");
	}

	std::string show_reserved_constants(void)
	{
		std::vector<std::string> reserved_constants_list(reserved_constants.size());

		for (int idx = 0; idx < reserved_constants.size(); idx++) {

			reserved_constants_list[idx] = reserved_constants[idx].first;
		}

		return combine(reserved_constants_list, ", ");
	}

	std::string show_variables(void)
	{
		return combine(variables, ", ");
	}

	std::string show_user_constants(void)
	{
		return combine(userconstants, ", ");
	}

	/////////////////////////////////////////////////////////
	//
	// LOGICAL EQUATION GETTERS

	std::vector< std::vector<EqComp::FSPEC> >& get_scalar_fspec(void) { return eq_component_1.get_fspec(); }

	std::vector<std::vector< std::vector<EqComp::FSPEC> >> get_dual_fspec(void)
	{
		return { eq_component_1.get_fspec(), eq_component_2.get_fspec() };
	}

	std::vector<std::vector< std::vector<EqComp::FSPEC> >> get_vector_fspec(void)
	{
		return { eq_component_1.get_fspec(), eq_component_2.get_fspec(), eq_component_3.get_fspec() };
	}

	/////////////////////////////////////////////////////////
	//
	// SET SPECIAL FUNCTIONS

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

		eq_component_1.Set_SpecialFunction(type, pSpecialFunc);
		eq_component_2.Set_SpecialFunction(type, pSpecialFunc);
		eq_component_3.Set_SpecialFunction(type, pSpecialFunc);
	}

	/////////////////////////////////////////////////////////
	//
	// EVALUATE EQUATION

	//evaluate scalar equation
	double evaluate(BVarType... bvars) const
	{
		return eq_component_1.evaluate(bvars...);
	}

	//evaluate dual equation
	DBL2 evaluate_dual(BVarType... bvars)
	{
		DBL2 value = DBL2(
			eq_component_1.evaluate(bvars...),
			eq_component_2.evaluate(bvars...)
		);

		return value;
	}

	//evaluate vector equation
	DBL3 evaluate_vector(BVarType... bvars)
	{
		DBL3 value = DBL3(
			eq_component_1.evaluate(bvars...),
			eq_component_2.evaluate(bvars...),
			eq_component_3.evaluate(bvars...)
		);

		return value;
	}
};
