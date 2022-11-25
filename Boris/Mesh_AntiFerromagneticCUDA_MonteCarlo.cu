#include "Mesh_AntiFerromagneticCUDA.h"

#if COMPILECUDA == 1

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void Zero_MCAux_Atom_Mesh_AFMCUDA(cuBReal& aux_real)
{
	if (threadIdx.x == 0) aux_real = 0.0;
}

///////////////////////////////////////////////////////////////
// PARALLEL MONTE-CARLO METROPOLIS - WITH REDUCTION

__global__ void Iterate_MonteCarloCUDA_Classic_AFM_red_kernel(ManagedMeshCUDA& cuMesh, int* cuModules, int& numModules, cuReal3& Ha, cuBorisRand<>& prng, cuBReal mc_cone_angledeg, cuBReal& mc_acceptance_rate)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	cuBReal T_Curie = *cuMesh.pT_Curie;

	cuSZ3& n = M.n;

	int num_moves = M.get_nonempty_cells();

	cuBReal acceptance_rate = 0.0;

	//this method must be called with half-size : n.dim() / 2, i.e. <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int spin_idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(spin_idx % n.x, (spin_idx / n.x) % n.y, spin_idx / (n.x*n.y));

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for red squares nudge on odd planes only
			spin_idx += (int)(n.z % 2);
		}
		//else : nx is odd and ny is odd, no nudge is needed
	}
	else {

		//nx is even : for red squares nudge on a) odd rows and even planes, and b) even rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool red_nudge = (((ijk.j % 2) == 1 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 0 && (ijk.k % 2) == 1)));

		spin_idx += (int)red_nudge;
	}

	//calculate only in non-empty and non-frozen cells
	if (spin_idx < n.dim() && M.is_not_empty(spin_idx) && !M.is_skipcell(spin_idx)) {

		cuReal2 MsAFM_val = *cuMesh.pMs_AFM;
		cuReal2 susrelAFM_val = *cuMesh.psusrel_AFM;
		cuMesh.update_parameters_mcoarse(spin_idx, *cuMesh.pMs_AFM, MsAFM_val, *cuMesh.psusrel_AFM, susrelAFM_val);

		cuBReal Temperature;
		if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(spin_idx)];
		else Temperature = *cuMesh.pbase_temperature;

		cuReal2 Ms0_AFM = cuMesh.pMs_AFM->get0();
		cuBReal me_A = MsAFM_val.i / Ms0_AFM.i;
		cuBReal me_B = MsAFM_val.j / Ms0_AFM.j;

		cuReal2 tau_ii_val = *cuMesh.ptau_ii;
		cuReal2 tau_ij_val = *cuMesh.ptau_ij;
		cuReal2 mu_val = *cuMesh.patomic_moment_AFM;

		////////////////////
		//
		// 1. Generate moved spins

		//Sub-lattice A : move spin

		//Picked spin is M[spin_idx]
		cuReal3 M_old_A = M[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		cuBReal theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		cuBReal phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_A = relrotate_polar(M_old_A, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_A*sqrt(susrelAFM_val.i*BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.i));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_A *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		//Sub-lattice B : move spin

		//Picked spin is M2[spin_idx]
		cuReal3 M_old_B = M2[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_B = relrotate_polar(M_old_B, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_B*sqrt(susrelAFM_val.j*BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.j));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_B *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		////////////////////
		//
		// 2. Energy change from sub-lattice A spin move

		//1. Find energy change
		cuReal2 energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M_new_A, M2[spin_idx], cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		cuReal3 m_A = M_old_A / Ms0_AFM.i;
		cuReal3 m_B = M_old_B / Ms0_AFM.j;
		cuReal3 m_new_A = M_new_A / Ms0_AFM.i;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_new_Asq = m_new_A * m_new_A;
			cuBReal diff_A_new = m_new_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (8 * susrelAFM_val.i * me_A*me_A)) *
				((diff_A_new * diff_A_new - diff_A * diff_A) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.j / (mu_val.i*MUB)) +
					diff_B * 6 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.i*(me_A / me_B) * (m_new_Asq - m_Asq) / (mu_val.i*MUB));

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal m_new_Asq = m_new_A * m_new_A;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (2 * susrelAFM_val.i)) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.i*MUB)) * (m_new_Asq - m_Asq);

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}

		////////////////////
		//
		// 3. Accept or reject sub-lattice A spin move

		//Compute acceptance probability
		cuBReal P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_A*M_new_A) / (M_old_A*M_old_A);
			P_accept = Mratio * Mratio * exp(-energy_delta.i / (BOLTZMANN * Temperature));
			//uniform random number between 0 and 1
			P = prng.rand();
		}
		else if (energy_delta < 0) P_accept = 1.0;

		if (P <= P_accept) {

			acceptance_rate = 1.0 / num_moves;

			//set new spin
			M[spin_idx] = M_new_A;
		}

		////////////////////
		//
		// 4. Energy change from sub-lattice B spin move

		//1. Find energy change
		energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M[spin_idx], M_new_B, cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		m_A = M[spin_idx] / Ms0_AFM.i;
		cuReal3 m_new_B = M_new_B / Ms0_AFM.j;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal m_new_Bsq = m_new_B * m_new_B;
			cuBReal diff_B_new = m_new_Bsq - me_B * me_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (8 * susrelAFM_val.j * me_B*me_B)) *
				((diff_B_new * diff_B_new - diff_B * diff_B) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.i / (mu_val.j*MUB)) +
					diff_A * 6 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.j*(me_B / me_A) * (m_new_Bsq - m_Bsq) / (mu_val.j*MUB));

			energy_delta += cuReal2(0.0, delta_elong_B);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Bsq = m_B * m_B;
			cuBReal m_new_Bsq = m_new_B * m_new_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (2 * susrelAFM_val.j)) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.j*MUB)) * (m_new_Bsq - m_Bsq);

			energy_delta += cuReal2(0.0, delta_elong_B);
		}

		////////////////////
		//
		// 5. Accept or reject sub-lattice B spin move

		//Compute acceptance probability
		P_accept = 0.0; P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_B*M_new_B) / (M_old_B*M_old_B);
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

	reduction_sum(0, 1, &acceptance_rate, mc_acceptance_rate);
}

__global__ void Iterate_MonteCarloCUDA_Classic_AFM_black_kernel(ManagedMeshCUDA& cuMesh, int* cuModules, int& numModules, cuReal3& Ha, cuBorisRand<>& prng, cuBReal mc_cone_angledeg, cuBReal& mc_acceptance_rate)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	cuBReal T_Curie = *cuMesh.pT_Curie;

	cuSZ3& n = M.n;

	int num_moves = M.get_nonempty_cells();

	cuBReal acceptance_rate = 0.0;

	//this method must be called with half-size : n.dim() / 2, i.e. <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int spin_idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(spin_idx % n.x, (spin_idx / n.x) % n.y, spin_idx / (n.x*n.y));

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for black squares nudge on even planes only
			spin_idx += (int)(n.z % 2 == 0);
		}
		else {

			//nx is odd and ny is odd, nudge everything by 1 for black squares
			spin_idx++;
		}
	}
	else {

		//nx is even : for black squares nudge on a) even rows and even planes, and b) odd rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool black_nudge = (((ijk.j % 2) == 0 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 1 && (ijk.k % 2) == 1)));

		spin_idx += (int)black_nudge;
	}

	//calculate only in non-empty and non-frozen cells
	if (spin_idx < n.dim() && M.is_not_empty(spin_idx) && !M.is_skipcell(spin_idx)) {

		cuReal2 MsAFM_val = *cuMesh.pMs_AFM;
		cuReal2 susrelAFM_val = *cuMesh.psusrel_AFM;
		cuMesh.update_parameters_mcoarse(spin_idx, *cuMesh.pMs_AFM, MsAFM_val, *cuMesh.psusrel_AFM, susrelAFM_val);

		cuBReal Temperature;
		if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(spin_idx)];
		else Temperature = *cuMesh.pbase_temperature;

		cuReal2 Ms0_AFM = cuMesh.pMs_AFM->get0();
		cuBReal me_A = MsAFM_val.i / Ms0_AFM.i;
		cuBReal me_B = MsAFM_val.j / Ms0_AFM.j;

		cuReal2 tau_ii_val = *cuMesh.ptau_ii;
		cuReal2 tau_ij_val = *cuMesh.ptau_ij;
		cuReal2 mu_val = *cuMesh.patomic_moment_AFM;

		////////////////////
		//
		// 1. Generate moved spins

		//Sub-lattice A : move spin

		//Picked spin is M[spin_idx]
		cuReal3 M_old_A = M[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		cuBReal theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		cuBReal phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_A = relrotate_polar(M_old_A, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_A*sqrt(susrelAFM_val.i*(cuBReal)BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.i));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_A *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		//Sub-lattice B : move spin

		//Picked spin is M2[spin_idx]
		cuReal3 M_old_B = M2[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_B = relrotate_polar(M_old_B, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {
			
			cuBReal sigma = 2 * me_B*sqrt(susrelAFM_val.j*BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.j));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_B *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		////////////////////
		//
		// 2. Energy change from sub-lattice A spin move

		//1. Find energy change
		cuReal2 energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M_new_A, M2[spin_idx], cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		cuReal3 m_A = M_old_A / Ms0_AFM.i;
		cuReal3 m_B = M_old_B / Ms0_AFM.j;
		cuReal3 m_new_A = M_new_A / Ms0_AFM.i;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_new_Asq = m_new_A * m_new_A;
			cuBReal diff_A_new = m_new_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (8 * susrelAFM_val.i * me_A*me_A)) *
				((diff_A_new * diff_A_new - diff_A * diff_A) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.j / (mu_val.i*MUB)) +
					diff_B * 6 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.i*(me_A / me_B) * (m_new_Asq - m_Asq) / (mu_val.i*MUB));

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal m_new_Asq = m_new_A * m_new_A;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (2 * susrelAFM_val.i)) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.i*MUB)) * (m_new_Asq - m_Asq);

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}

		////////////////////
		//
		// 3. Accept or reject sub-lattice A spin move

		//Compute acceptance probability
		cuBReal P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_A*M_new_A) / (M_old_A*M_old_A);
			P_accept = Mratio * Mratio * exp(-energy_delta.i / (BOLTZMANN * Temperature));
			//uniform random number between 0 and 1
			P = prng.rand();
		}
		else if (energy_delta < 0) P_accept = 1.0;

		if (P <= P_accept) {

			acceptance_rate = 1.0 / num_moves;

			//set new spin
			M[spin_idx] = M_new_A;
		}

		////////////////////
		//
		// 4. Energy change from sub-lattice B spin move

		//1. Find energy change
		energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M[spin_idx], M_new_B, cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		m_A = M[spin_idx] / Ms0_AFM.i;
		cuReal3 m_new_B = M_new_B / Ms0_AFM.j;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal m_new_Bsq = m_new_B * m_new_B;
			cuBReal diff_B_new = m_new_Bsq - me_B * me_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (8 * susrelAFM_val.j * me_B*me_B)) *
				((diff_B_new * diff_B_new - diff_B * diff_B) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.i / (mu_val.j*MUB)) +
					diff_A * 6 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.j*(me_B / me_A) * (m_new_Bsq - m_Bsq) / (mu_val.j*MUB));

			energy_delta += cuReal2(0.0, delta_elong_B);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Bsq = m_B * m_B;
			cuBReal m_new_Bsq = m_new_B * m_new_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (2 * susrelAFM_val.j)) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.j*MUB)) * (m_new_Bsq - m_Bsq);

			energy_delta += cuReal2(0.0, delta_elong_B);
		}

		////////////////////
		//
		// 5. Accept or reject sub-lattice B spin move

		//Compute acceptance probability
		P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_B*M_new_B) / (M_old_B*M_old_B);
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

	reduction_sum(0, 1, &acceptance_rate, mc_acceptance_rate);
}

///////////////////////////////////////////////////////////////
// PARALLEL MONTE-CARLO METROPOLIS - WITHOUT REDUCTION

__global__ void Iterate_MonteCarloCUDA_Classic_AFM_red_kernel(ManagedMeshCUDA& cuMesh, int* cuModules, int& numModules, cuReal3& Ha, cuBorisRand<>& prng, cuBReal mc_cone_angledeg)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	cuBReal T_Curie = *cuMesh.pT_Curie;

	cuSZ3& n = M.n;

	int num_moves = M.get_nonempty_cells();

	//this method must be called with half-size : n.dim() / 2, i.e. <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int spin_idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(spin_idx % n.x, (spin_idx / n.x) % n.y, spin_idx / (n.x*n.y));

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for red squares nudge on odd planes only
			spin_idx += (int)(n.z % 2);
		}
		//else : nx is odd and ny is odd, no nudge is needed
	}
	else {

		//nx is even : for red squares nudge on a) odd rows and even planes, and b) even rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool red_nudge = (((ijk.j % 2) == 1 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 0 && (ijk.k % 2) == 1)));

		spin_idx += (int)red_nudge;
	}

	//calculate only in non-empty and non-frozen cells
	if (spin_idx < n.dim() && M.is_not_empty(spin_idx) && !M.is_skipcell(spin_idx)) {

		cuReal2 MsAFM_val = *cuMesh.pMs_AFM;
		cuReal2 susrelAFM_val = *cuMesh.psusrel_AFM;
		cuMesh.update_parameters_mcoarse(spin_idx, *cuMesh.pMs_AFM, MsAFM_val, *cuMesh.psusrel_AFM, susrelAFM_val);

		cuBReal Temperature;
		if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(spin_idx)];
		else Temperature = *cuMesh.pbase_temperature;

		cuReal2 Ms0_AFM = cuMesh.pMs_AFM->get0();
		cuBReal me_A = MsAFM_val.i / Ms0_AFM.i;
		cuBReal me_B = MsAFM_val.j / Ms0_AFM.j;

		cuReal2 tau_ii_val = *cuMesh.ptau_ii;
		cuReal2 tau_ij_val = *cuMesh.ptau_ij;
		cuReal2 mu_val = *cuMesh.patomic_moment_AFM;

		////////////////////
		//
		// 1. Generate moved spins

		//Sub-lattice A : move spin

		//Picked spin is M[spin_idx]
		cuReal3 M_old_A = M[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		cuBReal theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		cuBReal phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_A = relrotate_polar(M_old_A, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_A*sqrt(susrelAFM_val.i*(cuBReal)BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.i));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_A *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		//Sub-lattice B : move spin

		//Picked spin is M2[spin_idx]
		cuReal3 M_old_B = M2[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_B = relrotate_polar(M_old_B, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_B*sqrt(susrelAFM_val.j*BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.j));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_B *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		////////////////////
		//
		// 2. Energy change from sub-lattice A spin move

		//1. Find energy change
		cuReal2 energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M_new_A, M2[spin_idx], cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		cuReal3 m_A = M_old_A / Ms0_AFM.i;
		cuReal3 m_B = M_old_B / Ms0_AFM.j;
		cuReal3 m_new_A = M_new_A / Ms0_AFM.i;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_new_Asq = m_new_A * m_new_A;
			cuBReal diff_A_new = m_new_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (8 * susrelAFM_val.i * me_A*me_A)) *
				((diff_A_new * diff_A_new - diff_A * diff_A) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.j / (mu_val.i*MUB)) +
					diff_B * 6 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.i*(me_A / me_B) * (m_new_Asq - m_Asq) / (mu_val.i*MUB));

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal m_new_Asq = m_new_A * m_new_A;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (2 * susrelAFM_val.i)) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.i*MUB)) * (m_new_Asq - m_Asq);

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}

		////////////////////
		//
		// 3. Accept or reject sub-lattice A spin move

		//Compute acceptance probability
		cuBReal P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_A*M_new_A) / (M_old_A*M_old_A);
			P_accept = Mratio * Mratio * exp(-energy_delta.i / (BOLTZMANN * Temperature));
			//uniform random number between 0 and 1
			P = prng.rand();
		}
		else if (energy_delta < 0) P_accept = 1.0;

		if (P <= P_accept) {

			//set new spin
			M[spin_idx] = M_new_A;
		}

		////////////////////
		//
		// 4. Energy change from sub-lattice B spin move

		//1. Find energy change
		energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M[spin_idx], M_new_B, cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		m_A = M[spin_idx] / Ms0_AFM.i;
		cuReal3 m_new_B = M_new_B / Ms0_AFM.j;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal m_new_Bsq = m_new_B * m_new_B;
			cuBReal diff_B_new = m_new_Bsq - me_B * me_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (8 * susrelAFM_val.j * me_B*me_B)) *
				((diff_B_new * diff_B_new - diff_B * diff_B) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.i / (mu_val.j*MUB)) +
					diff_A * 6 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.j*(me_B / me_A) * (m_new_Bsq - m_Bsq) / (mu_val.j*MUB));

			energy_delta += cuReal2(0.0, delta_elong_B);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Bsq = m_B * m_B;
			cuBReal m_new_Bsq = m_new_B * m_new_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (2 * susrelAFM_val.j)) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.j*MUB)) * (m_new_Bsq - m_Bsq);

			energy_delta += cuReal2(0.0, delta_elong_B);
		}

		////////////////////
		//
		// 5. Accept or reject sub-lattice B spin move

		//Compute acceptance probability
		P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_B*M_new_B) / (M_old_B*M_old_B);
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

__global__ void Iterate_MonteCarloCUDA_Classic_AFM_black_kernel(ManagedMeshCUDA& cuMesh, int* cuModules, int& numModules, cuReal3& Ha, cuBorisRand<>& prng, cuBReal mc_cone_angledeg)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC_VC<cuBReal>& Temp = *cuMesh.pTemp;

	cuBReal T_Curie = *cuMesh.pT_Curie;

	cuSZ3& n = M.n;

	int num_moves = M.get_nonempty_cells();

	//this method must be called with half-size : n.dim() / 2, i.e. <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int spin_idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(spin_idx % n.x, (spin_idx / n.x) % n.y, spin_idx / (n.x*n.y));

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for black squares nudge on even planes only
			spin_idx += (int)(n.z % 2 == 0);
		}
		else {

			//nx is odd and ny is odd, nudge everything by 1 for black squares
			spin_idx++;
		}
	}
	else {

		//nx is even : for black squares nudge on a) even rows and even planes, and b) odd rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool black_nudge = (((ijk.j % 2) == 0 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 1 && (ijk.k % 2) == 1)));

		spin_idx += (int)black_nudge;
	}

	//calculate only in non-empty and non-frozen cells
	if (spin_idx < n.dim() && M.is_not_empty(spin_idx) && !M.is_skipcell(spin_idx)) {

		cuReal2 MsAFM_val = *cuMesh.pMs_AFM;
		cuReal2 susrelAFM_val = *cuMesh.psusrel_AFM;
		cuMesh.update_parameters_mcoarse(spin_idx, *cuMesh.pMs_AFM, MsAFM_val, *cuMesh.psusrel_AFM, susrelAFM_val);

		cuBReal Temperature;
		if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(spin_idx)];
		else Temperature = *cuMesh.pbase_temperature;

		cuReal2 Ms0_AFM = cuMesh.pMs_AFM->get0();
		cuBReal me_A = MsAFM_val.i / Ms0_AFM.i;
		cuBReal me_B = MsAFM_val.j / Ms0_AFM.j;

		cuReal2 tau_ii_val = *cuMesh.ptau_ii;
		cuReal2 tau_ij_val = *cuMesh.ptau_ij;
		cuReal2 mu_val = *cuMesh.patomic_moment_AFM;

		////////////////////
		//
		// 1. Generate moved spins

		//Sub-lattice A : move spin

		//Picked spin is M[spin_idx]
		cuReal3 M_old_A = M[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		cuBReal theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		cuBReal phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_A = relrotate_polar(M_old_A, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_A*sqrt(susrelAFM_val.i*(cuBReal)BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.i));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_A *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		//Sub-lattice B : move spin

		//Picked spin is M2[spin_idx]
		cuReal3 M_old_B = M2[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
		phi_rot = prng.rand() * 2 * PI;
		//Move spin in cone with uniform random probability distribution.
		cuReal3 M_new_B = relrotate_polar(M_old_B, theta_rot, phi_rot);

		//now allow magnetization length to change slightly with a Gaussian pdf around current value with sigma value from the normal distribution of P(m^2).
		if (Temperature > 0.0) {

			cuBReal sigma = 2 * me_B*sqrt(susrelAFM_val.j*BOLTZMANN*Temperature / (M.h.dim() * Ms0_AFM.j));
			if (Temperature >= T_Curie || sigma > 0.03) sigma = 0.03;
			M_new_B *= 1 + (prng.rand() * 2 * sigma - sigma);
		}

		////////////////////
		//
		// 2. Energy change from sub-lattice A spin move

		//1. Find energy change
		cuReal2 energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M_new_A, M2[spin_idx], cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		cuReal3 m_A = M_old_A / Ms0_AFM.i;
		cuReal3 m_B = M_old_B / Ms0_AFM.j;
		cuReal3 m_new_A = M_new_A / Ms0_AFM.i;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_new_Asq = m_new_A * m_new_A;
			cuBReal diff_A_new = m_new_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (8 * susrelAFM_val.i * me_A*me_A)) *
				((diff_A_new * diff_A_new - diff_A * diff_A) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.j / (mu_val.i*MUB)) +
					diff_B * 6 * tau_ij_val.i*BOLTZMANN*T_Curie*susrelAFM_val.i*(me_A / me_B) * (m_new_Asq - m_Asq) / (mu_val.i*MUB));

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal m_new_Asq = m_new_A * m_new_A;

			cuBReal delta_elong_A = M.h.dim() * (Ms0_AFM.i / (2 * susrelAFM_val.i)) * (1 + 3 * tau_ij_val.i*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.i*MUB)) * (m_new_Asq - m_Asq);

			energy_delta += cuReal2(delta_elong_A, 0.0);
		}

		////////////////////
		//
		// 3. Accept or reject sub-lattice A spin move

		//Compute acceptance probability
		cuBReal P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_A*M_new_A) / (M_old_A*M_old_A);
			P_accept = Mratio * Mratio * exp(-energy_delta.i / (BOLTZMANN * Temperature));
			//uniform random number between 0 and 1
			P = prng.rand();
		}
		else if (energy_delta < 0) P_accept = 1.0;

		if (P <= P_accept) {

			//set new spin
			M[spin_idx] = M_new_A;
		}

		////////////////////
		//
		// 4. Energy change from sub-lattice B spin move

		//1. Find energy change
		energy_delta = cuMesh.Get_EnergyChange_AFM(spin_idx, M[spin_idx], M_new_B, cuModules, numModules, Ha);

		//2. Find contribution to free energy change from longitudinal susceptibility
		m_A = M[spin_idx] / Ms0_AFM.i;
		cuReal3 m_new_B = M_new_B / Ms0_AFM.j;

		if (Temperature > 0.0 && Temperature <= T_Curie) {

			cuBReal m_Asq = m_A * m_A;
			cuBReal diff_A = m_Asq - me_A * me_A;

			cuBReal m_Bsq = m_B * m_B;
			cuBReal diff_B = m_Bsq - me_B * me_B;

			cuBReal m_new_Bsq = m_new_B * m_new_B;
			cuBReal diff_B_new = m_new_Bsq - me_B * me_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (8 * susrelAFM_val.j * me_B*me_B)) *
				((diff_B_new * diff_B_new - diff_B * diff_B) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.i / (mu_val.j*MUB)) +
					diff_A * 6 * tau_ij_val.j*BOLTZMANN*T_Curie*susrelAFM_val.j*(me_B / me_A) * (m_new_Bsq - m_Bsq) / (mu_val.j*MUB));

			energy_delta += cuReal2(0.0, delta_elong_B);
		}
		else if (Temperature > 0.0) {

			cuBReal m_Bsq = m_B * m_B;
			cuBReal m_new_Bsq = m_new_B * m_new_B;

			cuBReal delta_elong_B = M.h.dim() * (Ms0_AFM.j / (2 * susrelAFM_val.j)) * (1 + 3 * tau_ij_val.j*BOLTZMANN*T_Curie*(susrelAFM_val.i + susrelAFM_val.j) / (mu_val.j*MUB)) * (m_new_Bsq - m_Bsq);

			energy_delta += cuReal2(0.0, delta_elong_B);
		}

		////////////////////
		//
		// 5. Accept or reject sub-lattice B spin move

		//Compute acceptance probability
		P_accept = 0.0, P = 1.0;
		if (Temperature > 0.0) {

			cuBReal Mratio = (M_new_B*M_new_B) / (M_old_B*M_old_B);
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

cuBReal AFMeshCUDA::Iterate_MonteCarloCUDA_Classic(cuBReal mc_cone_angledeg, double target_acceptance_rate)
{
	if (mc_acceptance_reduction_counter == 0) Zero_MCAux_Atom_Mesh_AFMCUDA <<< 1, CUDATHREADS >>> (mc_acceptance_rate);

	//Field set
	if (pHa) {

		if (mc_acceptance_reduction_counter == 0) {

			//with acceptance rate reduction
			Iterate_MonteCarloCUDA_Classic_AFM_red_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, *pHa, prng, mc_cone_angledeg, mc_acceptance_rate);
			Iterate_MonteCarloCUDA_Classic_AFM_black_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, *pHa, prng, mc_cone_angledeg, mc_acceptance_rate);
		}
		else {

			//without acceptance rate reduction
			Iterate_MonteCarloCUDA_Classic_AFM_red_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, *pHa, prng, mc_cone_angledeg);
			Iterate_MonteCarloCUDA_Classic_AFM_black_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, *pHa, prng, mc_cone_angledeg);
		}
	}
	//No field (or rather ZeemanCUDA module not added)
	else {

		cu_obj<cuReal3> Ha;
		Ha.from_cpu(cuReal3());

		if (mc_acceptance_reduction_counter == 0) {

			//with acceptance rate reduction
			Iterate_MonteCarloCUDA_Classic_AFM_red_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, Ha, prng, mc_cone_angledeg, mc_acceptance_rate);
			Iterate_MonteCarloCUDA_Classic_AFM_black_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, Ha, prng, mc_cone_angledeg, mc_acceptance_rate);
		}
		else {

			//without acceptance rate reduction
			Iterate_MonteCarloCUDA_Classic_AFM_red_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, Ha, prng, mc_cone_angledeg);
			Iterate_MonteCarloCUDA_Classic_AFM_black_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuMesh, cuModules, cuNumModules, Ha, prng, mc_cone_angledeg);
		}
	}

	if (mc_acceptance_reduction_counter == 0) {

		mc_acceptance_rate_last = mc_acceptance_rate.to_cpu();

		//is acceptance rate close enough to target acceptance? If yes don't do reduction next time.
		if (abs(target_acceptance_rate - mc_acceptance_rate_last) < MONTECARLO_ACCEPTANCETOLERANCE) mc_acceptance_reduction_counter = 1;

		return mc_acceptance_rate_last;
	}
	else {

		//increase counter until we come full circle : when mc_acceptance_reduction_counter becomes zero again we'll have to do reduction just to check.
		mc_acceptance_reduction_counter = (mc_acceptance_reduction_counter + 1) % MONTECARLO_REDUCTIONITERS;

		//return exact terget acceptance rate means cone angle will not be adjusted
		return target_acceptance_rate;
	}
}

#endif

#endif