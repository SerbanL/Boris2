#include "stdafx.h"
#include "Mesh_AntiFerromagnetic.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "MeshParamsControl.h"
#include "SuperMesh.h"

//Take a Monte Carlo step in this mesh
void AFMesh::Iterate_MonteCarlo(double acceptance_rate)
{
	if (mc_disabled) return;

	Iterate_MonteCarlo_Parallel_Classic();

	///////////////////////////////////////////////////////////////
	// ADAPTIVE CONE ANGLE

	MonteCarlo_AdaptiveAngle(mc_cone_angledeg, acceptance_rate);
}

#if COMPILECUDA == 1
//Take a Monte Carlo step in this mesh
void AFMesh::Iterate_MonteCarloCUDA(double acceptance_rate)
{
	if (mc_disabled) return;

	if (pMeshCUDA) {

		mc_acceptance_rate = reinterpret_cast<AFMeshCUDA*>(pMeshCUDA)->Iterate_MonteCarloCUDA_Classic(mc_cone_angledeg, acceptance_rate);
	}

	///////////////////////////////////////////////////////////////
	// ADAPTIVE CONE ANGLE

	MonteCarlo_AdaptiveAngle(mc_cone_angledeg, acceptance_rate);
}
#endif

void AFMesh::Iterate_MonteCarlo_Parallel_Classic(void)
{
	//number of moves in this step : one per each spin
	int num_moves = M.get_nonempty_cells();

	//number of cells in this mesh
	unsigned N = n.dim();

	//recalculate it
	mc_acceptance_rate = 0.0;

	///////////////////////////////////////////////////////////////
	// PARALLEL MONTE-CARLO METROPOLIS

	//make sure number of calls to rand doesn't match the prng period or we're asking for trouble with MC.
	prng.check_periodicity();

	//red-black : two passes will be taken
	int rb = 0;
	while (rb < 2) {

		double acceptance_rate = 0.0;

#pragma omp parallel for reduction(+:acceptance_rate)
		for (int idx_jk = 0; idx_jk < M.n.y * M.n.z; idx_jk++) {

			int j = idx_jk % M.n.y;
			int k = (idx_jk / M.n.y) % M.n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < M.n.x; i += 2) {

				int spin_idx = i + j * M.n.x + k * M.n.x*M.n.y;

				//only consider non-empty and non-frozen cells
				if (M.is_not_empty(spin_idx) && !M.is_skipcell(spin_idx)) {

					DBL2 MsAFM_val = Ms_AFM;
					DBL2 susrelAFM_val = susrel_AFM;
					update_parameters_mcoarse(spin_idx, Ms_AFM, MsAFM_val, susrel_AFM, susrelAFM_val);

					double Temperature;
					if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(spin_idx)];
					else Temperature = base_temperature;

					DBL2 Ms0_AFM = Ms_AFM.get0();
					double me_A = MsAFM_val.i / Ms0_AFM.i;
					double me_B = MsAFM_val.j / Ms0_AFM.j;

					DBL2 tau_ii_val = tau_ii;
					DBL2 tau_ij_val = tau_ij;
					DBL2 mu_val = atomic_moment_AFM;

					////////////////////
					//
					// 1. Generate moved spins

					//Sub-lattice A : move spin

					//Picked spin is M[spin_idx]
					DBL3 M_old_A = M[spin_idx];

					//obtain rotated spin in a cone around the picked spin
					double theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
					double phi_rot = prng.rand() * 2 * PI;
					//Move spin in cone with uniform random probability distribution.
					DBL3 M_new_A = relrotate_polar(M_old_A, theta_rot, phi_rot);

					//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
					if (Temperature > 0.0) {

						double sigma = 2 * me_A*sqrt(susrelAFM_val.i*BOLTZMANN*Temperature / (h.dim() * Ms0_AFM.i));
						if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
						M_new_A *= 1 + (prng.rand() * 2 * sigma - sigma);
					}

					//Sub-lattice B : move spin

					//Picked spin is M2[spin_idx]
					DBL3 M_old_B = M2[spin_idx];

					//obtain rotated spin in a cone around the picked spin
					theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
					phi_rot = prng.rand() * 2 * PI;
					//Move spin in cone with uniform random probability distribution.
					DBL3 M_new_B = relrotate_polar(M_old_B, theta_rot, phi_rot);

					//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
					if (Temperature > 0.0) {

						double sigma = 2 * me_B*sqrt(susrelAFM_val.j*BOLTZMANN*Temperature / (h.dim() * Ms0_AFM.j));
						if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
						M_new_B *= 1 + (prng.rand() * 2 * sigma - sigma);
					}

					////////////////////
					//
					// 2. Energy change from sub-lattice A spin move

					//1. Find energy change
					DBL2 energy_delta = DBL2();
					for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

						energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx, M_new_A, M2[spin_idx]);
					}

					//2. Find contribution to free energy change from longitudinal susceptibility
					DBL3 m_A = M_old_A / Ms0_AFM.i;
					DBL3 m_B = M_old_B / Ms0_AFM.j;
					DBL3 m_new_A = M_new_A / Ms0_AFM.i;

					if (Temperature > 0.0 && Temperature <= T_Curie) {

						double m_Asq = m_A * m_A;
						double diff_A = m_Asq - me_A * me_A;

						double m_new_Asq = m_new_A * m_new_A;
						double diff_A_new = m_new_Asq - me_A * me_A;

						double m_Bsq = m_B * m_B;
						double diff_B = m_Bsq - me_B * me_B;

						double delta_elong_A = h.dim() * (Ms0_AFM.i / (8 * susrelAFM_val.i * me_A*me_A)) * 
							((diff_A_new * diff_A_new - diff_A * diff_A) * (1 + 3*tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.j/(mu_val.i*MUB)) +
							diff_B * 6 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.i*(me_A / me_B) * (m_new_Asq - m_Asq) / (mu_val.i*MUB));

						energy_delta += DBL2(delta_elong_A, 0.0);
					}
					else if (Temperature > 0.0) {

						double m_Asq = m_A * m_A;
						double m_new_Asq = m_new_A * m_new_A;

						double delta_elong_A = h.dim() * (Ms0_AFM.i / (2 * susrelAFM_val.i)) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.i*MUB)) * (m_new_Asq - m_Asq);

						energy_delta += DBL2(delta_elong_A, 0.0);
					}

					////////////////////
					//
					// 3. Accept or reject sub-lattice A spin move

					//Compute acceptance probability
					double P_accept = 0.0, P = 1.0;
					if (Temperature > 0.0) {
						
						double Mratio = (M_new_A*M_new_A) / (M_old_A*M_old_A);
						P_accept = Mratio * Mratio * exp(-energy_delta.i / (BOLTZMANN * Temperature));
						//uniform random number between 0 and 1
						P = prng.rand();
					}
					else if (energy_delta < 0) P_accept = 1.0;

					if (P <= P_accept) {

						acceptance_rate += 1.0 / num_moves;

						//set new spin
						M[spin_idx] = M_new_A;
					}

					////////////////////
					//
					// 4. Energy change from sub-lattice B spin move

					//1. Find energy change
					energy_delta = DBL2(0);
					for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

						energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx, M[spin_idx], M_new_B);
					}

					//2. Find contribution to free energy change from longitudinal susceptibility
					m_A = M[spin_idx] / Ms0_AFM.i;
					DBL3 m_new_B = M_new_B / Ms0_AFM.j;

					if (Temperature > 0.0 && Temperature <= T_Curie) {

						double m_Asq = m_A * m_A;
						double diff_A = m_Asq - me_A * me_A;

						double m_Bsq = m_B * m_B;
						double diff_B = m_Bsq - me_B * me_B;

						double m_new_Bsq = m_new_B * m_new_B;
						double diff_B_new = m_new_Bsq - me_B * me_B;

						double delta_elong_B = h.dim() * (Ms0_AFM.j / (8 * susrelAFM_val.j * me_B*me_B)) * 
							((diff_B_new * diff_B_new - diff_B * diff_B) * (1 + 3*tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.i / (mu_val.j*MUB)) +
							diff_A * 6 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.j*(me_B / me_A) * (m_new_Bsq - m_Bsq) / (mu_val.j*MUB));

						energy_delta += DBL2(0.0, delta_elong_B);
					}
					else if (Temperature > 0.0) {

						double m_Bsq = m_B * m_B;
						double m_new_Bsq = m_new_B * m_new_B;

						double delta_elong_B = h.dim() * (Ms0_AFM.j / (2 * susrelAFM_val.j)) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.j*MUB)) * (m_new_Bsq - m_Bsq);

						energy_delta += DBL2(0.0, delta_elong_B);
					}

					////////////////////
					//
					// 5. Accept or reject sub-lattice B spin move

					//Compute acceptance probability
					P_accept = 0.0, P = 1.0;
					if (Temperature > 0.0) {

						double Mratio = (M_new_B*M_new_B) / (M_old_B*M_old_B);
						P_accept = Mratio * Mratio * exp(-energy_delta.j / (BOLTZMANN * Temperature));
						//uniform random number between 0 and 1
						P = prng.rand();
					}
					else if (energy_delta < 0) P_accept = 1.0;

					if (P <= P_accept) {

						//set new spin
						M2[spin_idx] = M_new_B;
					}
				}
			}
		}

		mc_acceptance_rate += acceptance_rate;
		rb++;
	}
}

#endif
