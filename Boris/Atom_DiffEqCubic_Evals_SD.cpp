#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_SD

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//Steepest-Descent solver with Barzilai-Borwein stepsizes (Journal of Numerical Analysis (1988) 8, 141-148).
//Use alternating stepsize selection. 
//Also tested variable threshold stepsize selection rules as in doi:10.1016/j.acha.2009.02.003, but decided against it as less reliable (although can be faster) and requiring parameters to tweak.
//SD solver for micromagnetics described in https://doi.org/10.1063/1.4862839 - the implementation here is similar but not identical.

//--------------------------------------------- Steepest Descent Solver

void Atom_DifferentialEquationCubic::RunSD_Start(void)
{
	//set new magnetization vectors
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				/////////////////////////

				double mu_s = paMesh->mu_s;
				double grel = paMesh->grel;
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->grel, grel);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				/////////////////////////

				//calculate m cross Heff1 (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				DBL3 mxHeff = (GAMMA / 2) * (m ^ H);

				/////////////////////////

				//current torque value G = m x (m x H)
				DBL3 G = m ^ mxHeff;

				//save calculated torque for next time
				sEval0[idx] = G;

				//save current m for next time
				sM1[idx] = m;

				/////////////////////////

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff1)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M1 value, Heff1 is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				//multiplication by grel important : can set grel to zero in entire mesh, or parts of mesh, which results in spins freezing.
				double s = dT * GAMMA * grel / 4.0;

				DBL3 s_mxH = m ^ (s * H);
				m = ((1 - (s_mxH*s_mxH)) * m - 2 * (m ^ s_mxH)) / (1 + (s_mxH*s_mxH));

				//set new M1
				paMesh->M1[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				paMesh->M1[idx].renormalize(mu_s);

				/////////////////////////
			}
		}
	}
}

//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
//must reset the static delta_... quantities before running these across all meshes
void Atom_DifferentialEquationCubic::RunSD_BB(void)
{
	double _delta_m_sq = 0.0;
	double _delta_G_sq = 0.0;
	double _delta_m_dot_delta_G = 0.0;

	//set new magnetization vectors
#pragma omp parallel for reduction(+:_delta_m_sq, _delta_G_sq, _delta_m_dot_delta_G)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				/////////////////////////

				double mu_s = paMesh->mu_s;
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				//calculate m cross Heff1 (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				DBL3 mxHeff = (GAMMA / 2) * (m ^ H);

				/////////////////////////

				//current torque value G = m x (m x H)
				DBL3 G = m ^ mxHeff;

				//change in torque
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur
				DBL3 delta_G = (G - sEval0[idx]) / 1e6;

				//save calculated torque for next time
				sEval0[idx] = G;

				/////////////////////////

				//change in M1
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				DBL3 delta_m = (m - sM1[idx]) / 1e6;

				//save current m for next time
				sM1[idx] = m;

				/////////////////////////

				//calculate num and denom for the two Barzilai-Borwein stepsize solutions (see Journal of Numerical Analysis (1988) 8, 141-148) so we can find new stepsize
				_delta_m_sq += delta_m * delta_m;
				_delta_G_sq += delta_G * delta_G;
				_delta_m_dot_delta_G += delta_m * delta_G;
			}
		}
	}

	//accumulate across all meshes -> remember these should have been set to zero before starting a run across all meshes
	delta_m_sq += _delta_m_sq;
	delta_G_sq += _delta_G_sq;
	delta_m_dot_delta_G += _delta_m_dot_delta_G;
}

//3. set new moment vectors
void Atom_DifferentialEquationCubic::RunSD_Advance_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();
	if (calculate_dmdt) dmdt_reduction.new_minmax_reduction();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

	//now we have the stepsize, set new M1 values
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				/////////////////////////

				double mu_s = paMesh->mu_s;
				double grel = paMesh->grel;
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->grel, grel);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				//obtained maximum normalized torque term
				if (IsNZ(grel)) {

					double _mxh = GetMagnitude(m ^ H) / (conversion * paMesh->M1[idx].norm());
					mxh_reduction.reduce_max(_mxh);
				}

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff1)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M1 value, Heff1 is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				//multiplication by grel important : can set grel to zero in entire mesh, or parts of mesh, which results in spins freezing.
				double s = dT * GAMMA * grel / 4.0;

				DBL3 s_mxH = m ^ (s * H);
				m = ((1 - (s_mxH*s_mxH)) * m - 2 * (m ^ s_mxH)) / (1 + (s_mxH*s_mxH));

				//set new M1
				paMesh->M1[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				paMesh->M1[idx].renormalize(mu_s);

				//use the flag check here to avoid doing dmdt reduction if not enabled (mxh and dmdt conditions are equivalent here, mxh is more likely to be used - at least by me!)
				if (calculate_dmdt && IsNZ(grel)) {

					//obtained maximum dmdt term
					double Mnorm = paMesh->M1[idx].norm();
					double _dmdt = GetMagnitude(paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * grel * Mnorm * conversion * Mnorm);
					dmdt_reduction.reduce_max(_dmdt);
				}
			}
		}
	}

	mxh_reduction.maximum();
	if (calculate_dmdt) dmdt_reduction.maximum();
}

void Atom_DifferentialEquationCubic::RunSD_Advance(void)
{
	//now we have the stepsize, set new M1 values
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				/////////////////////////

				double mu_s = paMesh->mu_s;
				double grel = paMesh->grel;
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s, paMesh->grel, grel);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff1)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M1 value, Heff1 is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly.

				//multiplication by grel important : can set grel to zero in entire mesh, or parts of mesh, which results in spins freezing.
				double s = dT * GAMMA * grel / 4.0;

				DBL3 s_mxH = m ^ (s * H);
				m = ((1 - (s_mxH*s_mxH)) * m - 2 * (m ^ s_mxH)) / (1 + (s_mxH*s_mxH));

				//set new M1
				paMesh->M1[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				paMesh->M1[idx].renormalize(mu_s);
			}
		}
	}
}

#endif
#endif