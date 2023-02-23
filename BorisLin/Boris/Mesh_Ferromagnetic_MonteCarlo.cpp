#include "stdafx.h"
#include "Mesh_Ferromagnetic.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "MeshParamsControl.h"
#include "SuperMesh.h"

//Take a Monte Carlo step in this mesh
void FMesh::Iterate_MonteCarlo(double acceptance_rate)
{
	if (mc_disabled) return;

	Iterate_MonteCarlo_Parallel_Classic();

	///////////////////////////////////////////////////////////////
	// ADAPTIVE CONE ANGLE

	MonteCarlo_AdaptiveAngle(mc_cone_angledeg, acceptance_rate);
}

#if COMPILECUDA == 1
//Take a Monte Carlo step in this mesh
void FMesh::Iterate_MonteCarloCUDA(double acceptance_rate)
{
	if (mc_disabled) return;

	if (pMeshCUDA) {

		mc_acceptance_rate = reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->Iterate_MonteCarloCUDA_Classic(mc_cone_angledeg, acceptance_rate);
	}

	///////////////////////////////////////////////////////////////
	// ADAPTIVE CONE ANGLE

	MonteCarlo_AdaptiveAngle(mc_cone_angledeg, acceptance_rate);
}
#endif

void FMesh::Iterate_MonteCarlo_Parallel_Classic(void)
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

			if (rb == 0) prng.rand();

			int j = idx_jk % M.n.y;
			int k = (idx_jk / M.n.y) % M.n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < M.n.x; i += 2) {

				int spin_idx = i + j * M.n.x + k * M.n.x*M.n.y;
				
				//only consider non-empty and non-frozen cells
				if (M.is_not_empty(spin_idx) && !M.is_skipcell(spin_idx)) {

					double Ms_val = Ms;
					double susrel_val = susrel;
					update_parameters_mcoarse(spin_idx, Ms, Ms_val, susrel, susrel_val);

					double Temperature;
					if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(spin_idx)];
					else Temperature = base_temperature;

					double Ms0 = Ms.get0();
					double me = Ms_val / Ms0;

					//Picked spin is M[spin_idx]
					DBL3 M_old = M[spin_idx];

					//obtain rotated spin in a cone around the picked spin
					double theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
					double phi_rot = prng.rand() * 2 * PI;
					//Move spin in cone with uniform random probability distribution.
					DBL3 M_new = relrotate_polar(M_old, theta_rot, phi_rot);

					//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
					if (Temperature > 0.0) {

						double sigma = 2 * me*sqrt(susrel_val*BOLTZMANN*Temperature / (h.dim() * Ms0));
						if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
						M_new *= 1 + (prng.rand() * 2 * sigma - sigma);
					}

					//1. Find energy change
					double energy_delta = 0.0;
					for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

						energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx, M_new);
					}

					//2. Find contribution to free energy change from longitudinal susceptibility
					DBL3 m = M_old / Ms0;
					DBL3 m_new = M_new / Ms0;

					if (Temperature > 0.0 && Temperature <= T_Curie) {

						double diff = m * m - me * me;
						double diff_new = m_new * m_new - me * me;

						energy_delta += h.dim() * (Ms0 / (8 * susrel_val * me*me)) * (diff_new * diff_new - diff * diff);
					}
					else if (Temperature > 0.0) {

						double r = 3 * T_Curie / (10 * (Temperature - T_Curie));
						double m_new_sq = m_new * m_new;
						double m_sq = m * m;
						energy_delta += h.dim() * (Ms0 / (2 * susrel_val)) * (m_new_sq * (1 + r * m_new_sq)- m_sq * (1 + r*m_sq));
					}

					//Compute acceptance probability
					//Target pdf is proportional to M^2 * exp(-E/kBT), however spin picking probability is not uniform, but proportional to M^2. Thus acceptance probability required to satisfy detailed balance is min{1, (M_new^4 / M_old^4) * exp(-dE/kBT)}
					double P_accept = 0.0, P = 1.0;
					if (Temperature > 0.0) {

						//Target pdf is proportional to M^2 * exp(-E/kBT), however spin picking probability is not uniform, but proportional to M^2. Thus acceptance probability required to satisfy detailed balance is min{1, (M_new^4 / M_old^4) * exp(-dE/kBT)}
						double Mratio = (M_new*M_new) / (M_old*M_old);
						P_accept = Mratio * Mratio * exp(-energy_delta / (BOLTZMANN * Temperature));
						//uniform random number between 0 and 1
						P = prng.rand();
					}
					else if (energy_delta < 0) P_accept = 1.0;

					if (P <= P_accept) {

						acceptance_rate += 1.0 / num_moves;

						//set new spin
						M[spin_idx] = M_new;
					}
				}	
			}
		}

		mc_acceptance_rate += acceptance_rate;
		rb++;
	}
}

#endif