#pragma once

#include "TEquation_Function.h"
#include "Funcs_Strings.h"

template <typename ... BVarType>
class Equation_Component 
{

private:

	//the equation component in logical form
	std::vector< std::vector<EqComp::FSPEC> > eq_fspec;

	//the equation component with Function objects ready for run-time evaluation
	std::vector< std::vector<EqComp::Function<BVarType...>*> > Funcs;

	//text specifiers for user-defined function variables
	std::vector<double>& varvec;

private:

	//check if branch in eq_fspec at idx index is just a constant
	bool is_constant_branch(int idx)
	{
		return (eq_fspec[idx].size() == 1 && eq_fspec[idx][0].is_constant());
	}

	//erase element from eq_fspec at given index
	void fspec_erase(int idx)
	{
		if (idx >= eq_fspec.size() || idx < 0) return;

		eq_fspec.erase(eq_fspec.begin() + idx);

		//now adjust all binary operators which have indexes of idx or above
		for (int op_idx = 0; op_idx < eq_fspec.size(); op_idx++) {

			if (eq_fspec[op_idx].size() && eq_fspec[op_idx][0].is_binary_operator()) {

				if (eq_fspec[op_idx][0].bin_idx1 >= idx) eq_fspec[op_idx][0].bin_idx1--;
				if (eq_fspec[op_idx][0].bin_idx2 >= idx) eq_fspec[op_idx][0].bin_idx2--;
			}
		}
	}

	//optimize eq_fspec so we reduce the number of operations
	void optimize(void)
	{
		//1. Search for FUNC_MUL (or FUNC_DIV):

		//1a. If FUNC_MUL (or FUNC_DIV) joins two constants only then replace it by their product and erase the constants
		//1b. If FUNC_MUL (or FUNC_DIV) joins a constant with anything else (unary function) then set the constant as parameter in last function of non-constant row, erase constant, and replace FUNC_MUL by the changed non-constant row

		//2. Search for FUNC_POW operators. If exactly one of either the basis or exponent are constants then convert to a unary function : delete the constant and ^ operator, and add a FUNC_POWER_EXPCONST or FUNC_POWER_BASECONST function at the end
		//of the non-constant branch, setting its base_or_exponent value.

		auto combine_constants = [&](void) -> bool {

			bool changes_made = false;

			int idx = 0;
			while (idx < eq_fspec.size()) {

				if (eq_fspec[idx][0].type == EqComp::FUNC_MUL || eq_fspec[idx][0].type == EqComp::FUNC_DIV) {

					int bin_idx1 = eq_fspec[idx][0].bin_idx1;
					int bin_idx2 = eq_fspec[idx][0].bin_idx2;

					//1a.
					if (is_constant_branch(bin_idx1) && is_constant_branch(bin_idx2)) {

						changes_made = true;

						//replace the multiplication (or division) with a combined constant
						if (eq_fspec[idx][0].type == EqComp::FUNC_MUL) {

							eq_fspec[idx][0] = EqComp::FSPEC(EqComp::FUNC_CONST, eq_fspec[bin_idx1][0].param * eq_fspec[bin_idx2][0].param);
						}
						else {

							eq_fspec[idx][0] = EqComp::FSPEC(EqComp::FUNC_CONST, eq_fspec[bin_idx1][0].param / eq_fspec[bin_idx2][0].param);
						}

						//erase larger index first (otherwise the other becomes invalid)
						if (bin_idx1 > bin_idx2) { fspec_erase(bin_idx1); fspec_erase(bin_idx2); }
						else { fspec_erase(bin_idx2); fspec_erase(bin_idx1); }

						idx -= 2;
						if (idx < 0) idx = 0;
						continue;
					}

					//1b.
					if (is_constant_branch(bin_idx1) || is_constant_branch(bin_idx2)) {

						//if FUNC_DIV then only proceed if bin_idx2 is the constant, since only then can we transform it to a multiplicative constant
						if (eq_fspec[idx][0].type == EqComp::FUNC_MUL || is_constant_branch(bin_idx2)) {

							changes_made = true;

							int bin_idx_const, bin_idx_nconst;
							if (is_constant_branch(bin_idx1)) { bin_idx_const = bin_idx1; bin_idx_nconst = bin_idx2; }
							else { bin_idx_const = bin_idx2; bin_idx_nconst = bin_idx1; }

							//set parameter in last function of non-constant branch
							if (eq_fspec[idx][0].type == EqComp::FUNC_MUL) {

								eq_fspec[bin_idx_nconst].back().make_pmul(eq_fspec[bin_idx_const][0].param);
							}
							else {

								eq_fspec[bin_idx_nconst].back().make_pmul(1 / eq_fspec[bin_idx_const][0].param);
							}

							//replace FUNC_MUL (or FUNC_DIV) with the non-constant branch
							eq_fspec[idx].insert(eq_fspec[idx].begin() + 1, eq_fspec[bin_idx_nconst].begin(), eq_fspec[bin_idx_nconst].end());
							eq_fspec[idx].erase(eq_fspec[idx].begin());

							//erase both the constant and non-constant branches
							if (bin_idx_const > bin_idx_nconst) { fspec_erase(bin_idx_const); fspec_erase(bin_idx_nconst); }
							else { fspec_erase(bin_idx_nconst); fspec_erase(bin_idx_const); }

							idx -= 2;
							if (idx < 0) idx = 0;
							continue;
						}
					}
				}
				
				//2
				if (eq_fspec[idx][0].type == EqComp::FUNC_POW) {

					int bin_idx1 = eq_fspec[idx][0].bin_idx1;
					int bin_idx2 = eq_fspec[idx][0].bin_idx2;

					if (is_constant_branch(bin_idx1) || is_constant_branch(bin_idx2)) {

						changes_made = true;

						//both are constant - reduce to a constant
						if (is_constant_branch(bin_idx1) && is_constant_branch(bin_idx2)) {

							//replace the power operation with a constant
							eq_fspec[idx][0] = EqComp::FSPEC(EqComp::FUNC_CONST, pow(eq_fspec[bin_idx1][0].param, eq_fspec[bin_idx2][0].param));

							//erase larger index first (otherwise the other becomes invalid)
							if (bin_idx1 > bin_idx2) { fspec_erase(bin_idx1); fspec_erase(bin_idx2); }
							else { fspec_erase(bin_idx2); fspec_erase(bin_idx1); }
						}
						else {

							if (is_constant_branch(bin_idx1)) { 
								
								//first replace the ^ operator by a unary function with the correct base_or_exponent value
								eq_fspec[idx][0] = EqComp::FSPEC(EqComp::FUNC_POWER_BASECONST);
								eq_fspec[idx][0].base_or_exponent = eq_fspec[bin_idx1][0].param;

								//copy this modified branch to the end of the non-constant branch
								eq_fspec[bin_idx2].insert(eq_fspec[bin_idx2].end(), eq_fspec[idx].begin(), eq_fspec[idx].end());

								//search for operators which work on idx and change their indexes to work in bin_idx2 instead : idx will soon be deleted and this branch was moved to bin_idx2
								for (int op_idx = 0; op_idx < eq_fspec.size(); op_idx++) {

									if (eq_fspec[op_idx].size() && eq_fspec[op_idx][0].is_binary_operator()) {

										if (eq_fspec[op_idx][0].bin_idx1 == idx) eq_fspec[op_idx][0].bin_idx1 = bin_idx2;
										if (eq_fspec[op_idx][0].bin_idx2 == idx) eq_fspec[op_idx][0].bin_idx2 = bin_idx2;
									}
								}

								//delete the constant and branch previously starting with the ^ operator
								if (bin_idx1 > idx) { fspec_erase(bin_idx1); fspec_erase(idx); }
								else { fspec_erase(idx); fspec_erase(bin_idx1); }
							}
							else { 
							
								//first replace the ^ operator by a unary function with the correct base_or_exponent value
								eq_fspec[idx][0] = EqComp::FSPEC(EqComp::FUNC_POWER_EXPCONST);
								eq_fspec[idx][0].base_or_exponent = eq_fspec[bin_idx2][0].param;

								//copy this modified branch to the end of the non-constant branch
								eq_fspec[bin_idx1].insert(eq_fspec[bin_idx1].end(), eq_fspec[idx].begin(), eq_fspec[idx].end());

								//search for operators which work on idx and change their indexes to work in bin_idx1 instead : idx will soon be deleted and this branch was moved to bin_idx1
								for (int op_idx = 0; op_idx < eq_fspec.size(); op_idx++) {

									if (eq_fspec[op_idx].size() && eq_fspec[op_idx][0].is_binary_operator()) {

										if (eq_fspec[op_idx][0].bin_idx1 == idx) eq_fspec[op_idx][0].bin_idx1 = bin_idx1;
										if (eq_fspec[op_idx][0].bin_idx2 == idx) eq_fspec[op_idx][0].bin_idx2 = bin_idx1;
									}
								}

								//delete the constant and branch previously starting with the ^ operator
								if (bin_idx2 > idx) { fspec_erase(bin_idx2); fspec_erase(idx); }
								else { fspec_erase(idx); fspec_erase(bin_idx2); }
							}

							idx -= 2;
							if (idx < 0) idx = 0;
							continue;
						}
					}
				}
				
				idx++;
			}

			return changes_made;
		};

		while (combine_constants());
	}

public:

	Equation_Component(std::vector<double>& varvec_) :
		varvec(varvec_)
	{}

	/////////////////////////////////////////////////////////
	//
	// CLEAR MEMORY

	void clear_Funcs(void)
	{
		for (int idx = 0; idx < Funcs.size(); idx++) {

			for (int idx_tree = 0; idx_tree < Funcs[idx].size(); idx_tree++) {

				delete Funcs[idx][idx_tree];
			}

			Funcs[idx].clear();
		}

		Funcs.clear();
		eq_fspec.clear();
	}

	/////////////////////////////////////////////////////////
	//
	// CREATE ACTUAL EQUATION WITH FUNCTION OBJECTS (eq_fspec should be set)

	//make Equation using eq_fspec logical representation
	bool make(void)
	{
		//eq_fspec may need to be "fixed" before using it to create the equation with Function objects.
		//if eq_fspec starts with a unary operator, attach the branch to the previous branch. All branches must either start with a basis function or with a binary operator which combines 2 branches.
		int idx = 0;
		while (idx < eq_fspec.size()) {

			if (!eq_fspec[idx].size()) {

				fspec_erase(idx);
				continue;
			}

			if (eq_fspec[idx][0].is_unary_function() && idx > 0) {

				eq_fspec[idx - 1].insert(eq_fspec[idx - 1].end(), eq_fspec[idx].begin(), eq_fspec[idx].end());
				fspec_erase(idx);

				continue;
			}

			idx++;
		}

		//optimize fspec to reduce number of operations where possible
		optimize();

		//now make Function objects
		for (int idx = 0; idx < eq_fspec.size(); idx++) {

			//new branch
			Funcs.push_back(std::vector<EqComp::Function<BVarType...>*>());

			if (eq_fspec[idx][0].is_binary_operator()) {

				int bin_idx1 = eq_fspec[idx][0].bin_idx1;
				int bin_idx2 = eq_fspec[idx][0].bin_idx2;

				if (bin_idx1 >= Funcs.size() || bin_idx2 >= Funcs.size()) return false;

				Funcs.back().push_back(new EqComp::Function<BVarType...>(eq_fspec[idx][0], *Funcs[bin_idx1].back(), *Funcs[bin_idx2].back(), varvec));
			}
			else {

				Funcs.back().push_back(new EqComp::Function<BVarType...>(eq_fspec[idx][0], varvec));
			}

			for (int idx_tree = 1; idx_tree < eq_fspec[idx].size(); idx_tree++) {

				Funcs.back().push_back(new EqComp::Function<BVarType...>(eq_fspec[idx][idx_tree], *Funcs.back()[idx_tree - 1], varvec));
			}
		}

		return true;
	}

	/////////////////////////////////////////////////////////
	//
	// IS THIS EQUATION SET?

	bool is_set(void) const { return Funcs.size(); }

	/////////////////////////////////////////////////////////
	//
	// GET FSPEC

	std::vector< std::vector<EqComp::FSPEC> >& get_fspec(void) { return eq_fspec; }

	/////////////////////////////////////////////////////////
	//
	// SET SPECIAL FUNCTIONS

	void Set_SpecialFunction(EqComp::FUNC_ type, std::shared_ptr<Funcs_Special> pSpecialFunc)
	{
		for (int idx = 0; idx < Funcs.size(); idx++) {

			for (int idx_tree = 0; idx_tree < Funcs[idx].size(); idx_tree++) {

				Funcs[idx][idx_tree]->Set_SpecialFunction(type, pSpecialFunc);
			}
		}
	}

	/////////////////////////////////////////////////////////
	//
	// EVALUATE EQUATION

	//evaluate scalar equation
	double evaluate(BVarType... bvars) const
	{
		if (Funcs.size()) return Funcs.back().back()->evaluate(bvars...);
		else return 0.0;
	}
};