#include "stdafx.h"
#include "DiffEqAFM.h"
#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_SD

//Steepest-Descent solver with Barzilai-Borwein stepsizes (Journal of Numerical Analysis (1988) 8, 141-148).
//Use alternating stepsize selection. 
//Also tested variable threshold stepsize selection rules as in doi:10.1016/j.acha.2009.02.003, but decided against it as less reliable (although can be faster) and requiring parameters to tweak.
//SD solver for micromagnetics described in https://doi.org/10.1063/1.4862839 - the implementation here is similar but not identical.

//--------------------------------------------- Steepest Descent Solver

void DifferentialEquationAFM::RunSD_Start(void)
{
	//set new magnetisation vectors
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				/////////////////////////

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);

				DBL3 m = pMesh->M[idx] / Ms_AFM.i;
				DBL3 H = pMesh->Heff[idx];

				DBL3 m2 = pMesh->M2[idx] / Ms_AFM.j;
				DBL3 H2 = pMesh->Heff2[idx];

				/////////////////////////

				//calculate M cross Heff (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				DBL3 MxHeff = (GAMMA / 2) * (pMesh->M[idx] ^ H);
				DBL3 MxHeff2 = (GAMMA / 2) * (pMesh->M2[idx] ^ H2);

				/////////////////////////

				//current torque value G = m x (M x H)
				DBL3 G = m ^ MxHeff;
				DBL3 G2 = m2 ^ MxHeff2;

				//save calculated torque for next time
				sEval0[idx] = G;
				sEval0_2[idx] = G2;

				//save current M for next time
				sM1[idx] = pMesh->M[idx];
				sM1_2[idx] = pMesh->M2[idx];

				/////////////////////////

				//The updating formula is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above formula can be solved for m_next explicitly (m_next . m = 1, i.e. norm conserved) as:
				//
				//Let s = dT * GAMMA / 4, the stepsize; H = Heff. Then:
				//m = (m - s(m.H * (m - 2sH) + 2sm(H.H) - 2H)) / (1 + sm.H)

				double mH = m * H;
				double mH2 = m2 * H2;

				double s = dT * GAMMA / 4.0;

				m = (m - s*(mH * (m - 2*s*H) + 2*s*m*(H*H) - 2*H)) / (1 + s*mH);
				m2 = (m2 - s*(mH2 * (m2 - 2*s*H2) + 2*s*m2*(H2*H2) - 2*H2)) / (1 + s*mH2);

				//set new M
				pMesh->M[idx] = m * Ms_AFM.i;
				pMesh->M2[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);

				/////////////////////////
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
//must reset the static delta_... quantities before running these across all meshes
void DifferentialEquationAFM::RunSD_BB(void)
{
	double _delta_M_sq = 0.0;
	double _delta_G_sq = 0.0;
	double _delta_M_dot_delta_G = 0.0;

	double _delta_M2_sq = 0.0;
	double _delta_G2_sq = 0.0;
	double _delta_M2_dot_delta_G2 = 0.0;

	//set new magnetisation vectors
#pragma omp parallel for reduction(+:_delta_M_sq, _delta_G_sq, _delta_M_dot_delta_G, _delta_M2_sq, _delta_G2_sq, _delta_M2_dot_delta_G2)
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				/////////////////////////

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);

				DBL3 m = pMesh->M[idx] / Ms_AFM.i;
				DBL3 H = pMesh->Heff[idx];

				DBL3 m2 = pMesh->M2[idx] / Ms_AFM.j;
				DBL3 H2 = pMesh->Heff2[idx];

				//calculate M cross Heff (multiplication by GAMMA/2 not necessary as this could be absorbed in the stepsize, but keep it for a more natural step size value from the user point of view - i.e. a time step).
				DBL3 MxHeff = (GAMMA / 2) * (pMesh->M[idx] ^ H);
				DBL3 MxHeff2 = (GAMMA / 2) * (pMesh->M2[idx] ^ H2);

				/////////////////////////

				//current torque value G = m x (M x H)
				DBL3 G = m ^ MxHeff;
				DBL3 G2 = m2 ^ MxHeff2;

				//change in torque
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				DBL3 delta_G = (G - sEval0[idx]) / 1e6;
				DBL3 delta_G2 = (G2 - sEval0_2[idx]) / 1e6;

				//save calculated torque for next time
				sEval0[idx] = G;
				sEval0_2[idx] = G2;

				/////////////////////////

				//change in M
				//divide by 1e6 to stop the accumulated value having a large exponent -> both num and denom are divided by same value; if exponent too large when dividing num by denom significant loss of precision can occur.
				//Also you don't want to normalize to Ms since Ms can vary between different meshes, or even in this same mesh.
				DBL3 delta_M = (pMesh->M[idx] - sM1[idx]) / 1e6;
				DBL3 delta_M2 = (pMesh->M2[idx] - sM1_2[idx]) / 1e6;

				//save current M for next time
				sM1[idx] = pMesh->M[idx];
				sM1_2[idx] = pMesh->M2[idx];

				/////////////////////////

				//calculate num and denom for the two Barzilai-Borwein stepsize solutions (see Journal of Numerical Analysis (1988) 8, 141-148) so we can find new stepsize
				_delta_M_sq += delta_M * delta_M;
				_delta_G_sq += delta_G * delta_G;
				_delta_M_dot_delta_G += delta_M * delta_G;

				_delta_M2_sq += delta_M2 * delta_M2;
				_delta_G2_sq += delta_G2 * delta_G2;
				_delta_M2_dot_delta_G2 += delta_M2 * delta_G2;
			}
		}
	}

	//accumulate across all meshes -> remember these should have been set to zero before starting a run across all meshes
	delta_M_sq += _delta_M_sq;
	delta_G_sq += _delta_G_sq;
	delta_M_dot_delta_G += _delta_M_dot_delta_G;

	delta_M2_sq += _delta_M2_sq;
	delta_G2_sq += _delta_G2_sq;
	delta_M2_dot_delta_G2 += _delta_M2_dot_delta_G2;
}

//3. set new magnetization vectors
void DifferentialEquationAFM::RunSD_Advance_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();
	if (calculate_dmdt) dmdt_reduction.new_minmax_reduction();

	//now we have the stepsize, set new M values
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				/////////////////////////

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);

				DBL3 m = pMesh->M[idx] / Ms_AFM.i;
				DBL3 H = pMesh->Heff[idx];

				DBL3 m2 = pMesh->M2[idx] / Ms_AFM.j;
				DBL3 H2 = pMesh->Heff2[idx];

				//obtained maximum normalized torque term
				double _mxh = GetMagnitude(m ^ H) / pMesh->M[idx].norm();
				mxh_reduction.reduce_max(_mxh);

				//The updating formula is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above formula can be solved for m_next explicitly (m_next . m = 1, i.e. norm conserved) as:
				//
				//Let s = dT * GAMMA / 4, the stepsize; H = Heff. Then:
				//m = (m - s(m.H * (m - 2sH) + 2sm(H.H) - 2H)) / (1 + sm.H)

				double mH = m * H;
				double mH2 = m2 * H2;

				double s = dT * GAMMA / 4.0;

				m = (m - s * (mH * (m - 2 * s*H) + 2 * s*m*(H*H) - 2 * H)) / (1 + s * mH);
				m2 = (m2 - s * (mH2 * (m2 - 2 * s*H2) + 2 * s*m2*(H2*H2) - 2 * H2)) / (1 + s * mH2);

				//set new M
				pMesh->M[idx] = m * Ms_AFM.i;
				pMesh->M2[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);

				//use the flag check here to avoid doing dmdt reduction if not enabled (mxh and dmdt conditions are equivalent here, mxh is more likely to be used - at least by me!)
				if (calculate_dmdt) {

					//obtained maximum dmdt term
					double Mnorm = pMesh->M[idx].norm();
					double _dmdt = GetMagnitude(pMesh->M[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * Mnorm);
					dmdt_reduction.reduce_max(_dmdt);
				}
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics are disabled in this mesh)
		mxh_reduction.maximum();
		if (calculate_dmdt) dmdt_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
		if (calculate_dmdt) dmdt_reduction.max = 0.0;
	}
}

void DifferentialEquationAFM::RunSD_Advance(void)
{
	//now we have the stepsize, set new M values
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				/////////////////////////

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);

				DBL3 m = pMesh->M[idx] / Ms_AFM.i;
				DBL3 H = pMesh->Heff[idx];

				DBL3 m2 = pMesh->M2[idx] / Ms_AFM.j;
				DBL3 H2 = pMesh->Heff2[idx];

				//The updating formula is (see https://doi.org/10.1063/1.4862839):

				//m_next = m - (dT/2) * (m_next + m) x ((gamma/2)m x Heff)
				//Here gamma = mu0 * |gamma_e| as usual., m is the current normalized M value, Heff is the current effective field, and we need to find m_next.
				//This is applicable to the LLGStatic approach, i.e. no precession term and damping set to 1.
				//M_next = m_next * Ms

				//The above formula can be solved for m_next explicitly (m_next . m = 1, i.e. norm conserved) as:
				//
				//Let s = dT * GAMMA / 4, the stepsize; H = Heff. Then:
				//m = (m - s(m.H * (m - 2sH) + 2sm(H.H) - 2H)) / (1 + sm.H)

				double mH = m * H;
				double mH2 = m2 * H2;

				double s = dT * GAMMA / 4.0;

				m = (m - s * (mH * (m - 2 * s*H) + 2 * s*m*(H*H) - 2 * H)) / (1 + s * mH);
				m2 = (m2 - s * (mH2 * (m2 - 2 * s*H2) + 2 * s*m2*(H2*H2) - 2 * H2)) / (1 + s * mH2);

				//set new M
				pMesh->M[idx] = m * Ms_AFM.i;
				pMesh->M2[idx] = m2 * Ms_AFM.j;

				//renormalize - method is supposed to conserve norm, but best to renormalize anyway.
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

#endif