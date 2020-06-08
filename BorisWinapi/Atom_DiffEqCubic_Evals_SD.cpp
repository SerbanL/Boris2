#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_SD

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
	//set new magnetisation vectors
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				/////////////////////////

				double mu_s = paMesh->mu_s;
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				/////////////////////////

				//calculate M1 cross Heff1 (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				DBL3 MxHeff = (GAMMA / 2) * (paMesh->M1[idx] ^ H);

				/////////////////////////

				//current torque value G = m x (M1 x H)
				DBL3 G = m ^ MxHeff;

				//save calculated torque for next time
				sEval0[idx] = G;

				//save current M1 for next time
				sM1[idx] = paMesh->M1[idx];

				/////////////////////////

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff1)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M1 value, Heff1 is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly (m_next . m = 1, i.e. norm conserved) as:
				//
				//Let s = dT * GAMMA / 4, the stepsize; H = Heff1. Then:
				//m = (m - s(m.H * (m - 2sH) + 2sm(H.H) - 2H)) / (1 + sm.H)

				double mH = m * H;

				double s = dT * GAMMA / 4.0;

				m = (m - s*(mH * (m - 2*s*H) + 2*s*m*(H*H) - 2*H)) / (1 + s*mH);

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
	double _delta_M_sq = 0.0;
	double _delta_G_sq = 0.0;
	double _delta_M_dot_delta_G = 0.0;

	//set new magnetisation vectors
#pragma omp parallel for reduction(+:_delta_M_sq, _delta_G_sq, _delta_M_dot_delta_G)
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				/////////////////////////

				double mu_s = paMesh->mu_s;
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				//calculate M1 cross Heff1 (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				DBL3 MxHeff = (GAMMA / 2) * (paMesh->M1[idx] ^ H);

				/////////////////////////

				//current torque value G = m x (M1 x H)
				DBL3 G = m ^ MxHeff;

				//change in torque
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				DBL3 delta_G = (G - sEval0[idx]) / 1e6;

				//save calculated torque for next time
				sEval0[idx] = G;

				/////////////////////////

				//change in M1
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				DBL3 delta_M = (paMesh->M1[idx] - sM1[idx]) / 1e6;

				//save current M1 for next time
				sM1[idx] = paMesh->M1[idx];

				/////////////////////////

				//calculate num and denom for the two Barzilai-Borwein stepsize solutions (see Journal of Numerical Analysis (1988) 8, 141-148) so we can find new stepsize
				_delta_M_sq += delta_M * delta_M;
				_delta_G_sq += delta_G * delta_G;
				_delta_M_dot_delta_G += delta_M * delta_G;
			}
		}
	}

	//accumulate across all meshes -> remember these should have been set to zero before starting a run across all meshes
	delta_M_sq += _delta_M_sq;
	delta_G_sq += _delta_G_sq;
	delta_M_dot_delta_G += _delta_M_dot_delta_G;
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
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				//obtained maximum normalized torque term
				double _mxh = GetMagnitude(m ^ H) / (conversion * paMesh->M1[idx].norm());
				mxh_reduction.reduce_max(_mxh);

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff1)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M1 value, Heff1 is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly (m_next . m = 1, i.e. norm conserved) as:
				//
				//Let s = dT * GAMMA / 4, the stepsize; H = Heff1. Then:
				//m = (m - s(m.H * (m - 2sH) + 2sm(H.H) - 2H)) / (1 + sm.H)

				double mH = m * H;

				double s = dT * GAMMA / 4.0;

				m = (m - s * (mH * (m - 2 * s*H) + 2 * s*m*(H*H) - 2 * H)) / (1 + s * mH);

				//set new M1
				paMesh->M1[idx] = m * mu_s;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				paMesh->M1[idx].renormalize(mu_s);

				//use the flag check here to avoid doing dmdt reduction if not enabled (mxh and dmdt conditions are equivalent here, mxh is more likely to be used - at least by me!)
				if (calculate_dmdt) {

					//obtained maximum dmdt term
					double Mnorm = paMesh->M1[idx].norm();
					double _dmdt = GetMagnitude(paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm);
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
				paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);

				DBL3 m = paMesh->M1[idx] / mu_s;
				DBL3 H = paMesh->Heff1[idx];

				//The updating equation is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff1)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M1 value, Heff1 is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above equation can be solved for m_next explicitly (m_next . m = 1, i.e. norm conserved) as:
				//
				//Let s = dT * GAMMA / 4, the stepsize; H = Heff1. Then:
				//m = (m - s(m.H * (m - 2sH) + 2sm(H.H) - 2H)) / (1 + sm.H)

				double mH = m * H;

				double s = dT * GAMMA / 4.0;

				m = (m - s * (mH * (m - 2 * s*H) + 2 * s*m*(H*H) - 2 * H)) / (1 + s * mH);

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