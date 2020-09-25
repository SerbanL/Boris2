#include "Atom_Mesh_CubicCUDA.h"

#if COMPILECUDA == 1

#ifdef MESH_COMPILATION_ATOM_CUBIC && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"

__global__ void Zero_MCAux_Atom_Mesh_CubicCUDA(cuBReal& aux_real)
{
	if (threadIdx.x == 0) aux_real = 0.0;
}

///////////////////////////////////////////////////////////////
// PARALLEL MONTE-CARLO METROPOLIS

__global__ void Iterate_MonteCarloCUDA_Classic_Cubic_red_kernel(ManagedAtom_MeshCUDA& cuaMesh, int* cuaModules, int& numModules, cuReal3& Ha, cuBorisRand& prng, cuBReal mc_cone_angledeg, cuBReal& mc_acceptance_rate)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	cuSZ3& n = M1.n;

	cuBReal base_temperature = *cuaMesh.pbase_temperature;

	int num_moves = M1.get_nonempty_cells();

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
	
	//calculate only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	if (spin_idx < n.dim() && M1.is_not_empty(spin_idx)) {

		//Picked spin is M1[spin_idx]
		cuReal3 M1_old = M1[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		cuBReal theta_rot = prng.rand() * mc_cone_angledeg * (cuBReal)PI / 180.0;
		cuBReal phi_rot = prng.rand() * 2 * (cuBReal)PI;
		//Move spin in cone with uniform random probability distribution. This approach only requires 2 random numbers to be generated. 
		//Also using a Gaussian distribution to move spin around the initial spin is less efficient, requiring more steps to thermalize.
		cuReal3 M1_new = relrotate_polar(M1_old, theta_rot, phi_rot);

		//Find energy for current spin by summing all contributing modules energies at given spin only
		cuBReal old_energy = cuaMesh.Get_Atomistic_Energy_SC(spin_idx, cuaModules, numModules, Ha);

		//Set new spin; M1_new has same length as M1[spin_idx], but reset length anyway to avoid floating point error creep
		M1[spin_idx] = M1_new.normalized() * M1_old.norm();
		cuBReal new_energy = cuaMesh.Get_Atomistic_Energy_SC(spin_idx, cuaModules, numModules, Ha);

		//Compute acceptance probability
		cuBReal P_accept = exp(-(new_energy - old_energy) / ((cuBReal)BOLTZMANN * base_temperature));

		//uniform random number between 0 and 1
		cuBReal P = prng.rand();

		if (P > P_accept) {

			//reject move
			M1[spin_idx] = M1_old;
		}
		else acceptance_rate = 1.0 / num_moves;
	}

	reduction_sum(0, 1, &acceptance_rate, mc_acceptance_rate);
}

__global__ void Iterate_MonteCarloCUDA_Classic_Cubic_black_kernel(ManagedAtom_MeshCUDA& cuaMesh, int* cuaModules, int& numModules, cuReal3& Ha, cuBorisRand& prng, cuBReal mc_cone_angledeg, cuBReal& mc_acceptance_rate)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	cuSZ3& n = M1.n;

	cuBReal base_temperature = *cuaMesh.pbase_temperature;

	int num_moves = M1.get_nonempty_cells();

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

	//calculate only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	if (spin_idx < n.dim() && M1.is_not_empty(spin_idx)) {

		//Picked spin is M1[spin_idx]
		cuReal3 M1_old = M1[spin_idx];

		//obtain rotated spin in a cone around the picked spin
		cuBReal theta_rot = prng.rand() * mc_cone_angledeg * (cuBReal)PI / 180.0;
		cuBReal phi_rot = prng.rand() * 2 * (cuBReal)PI;
		//Move spin in cone with uniform random probability distribution. This approach only requires 2 random numbers to be generated. 
		//Also using a Gaussian distribution to move spin around the initial spin is less efficient, requiring more steps to thermalize.
		cuReal3 M1_new = relrotate_polar(M1_old, theta_rot, phi_rot);

		//Find energy for current spin by summing all contributing modules energies at given spin only
		cuBReal old_energy = cuaMesh.Get_Atomistic_Energy_SC(spin_idx, cuaModules, numModules, Ha);

		//Set new spin; M1_new has same length as M1[spin_idx], but reset length anyway to avoid floating point error creep
		M1[spin_idx] = M1_new.normalized() * M1_old.norm();
		cuBReal new_energy = cuaMesh.Get_Atomistic_Energy_SC(spin_idx, cuaModules, numModules, Ha);

		//Compute acceptance probability
		cuBReal P_accept = exp(-(new_energy - old_energy) / ((cuBReal)BOLTZMANN * base_temperature));

		//uniform random number between 0 and 1
		cuBReal P = prng.rand();

		if (P > P_accept) {

			//reject move
			M1[spin_idx] = M1_old;
		}
		else acceptance_rate = 1.0 / num_moves;
	}

	reduction_sum(0, 1, &acceptance_rate, mc_acceptance_rate);
}

cuBReal Atom_Mesh_CubicCUDA::Iterate_MonteCarloCUDA_Classic(cuBReal mc_cone_angledeg)
{
	Zero_MCAux_Atom_Mesh_CubicCUDA <<< 1, CUDATHREADS >>> (mc_acceptance_rate);

	//Field set
	if (pHa) {

		Iterate_MonteCarloCUDA_Classic_Cubic_red_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaMesh, cuaModules, cuaNumModules, *pHa, prng, mc_cone_angledeg, mc_acceptance_rate);
		Iterate_MonteCarloCUDA_Classic_Cubic_black_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaMesh, cuaModules, cuaNumModules, *pHa, prng, mc_cone_angledeg, mc_acceptance_rate);
	}
	//No field (or rather Atom_ZeemanCUDA module not added)
	else {

		cu_obj<cuReal3> Ha;
		Ha.from_cpu(cuReal3());

		Iterate_MonteCarloCUDA_Classic_Cubic_red_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaMesh, cuaModules, cuaNumModules, Ha, prng, mc_cone_angledeg, mc_acceptance_rate);
		Iterate_MonteCarloCUDA_Classic_Cubic_black_kernel <<< (n.dim() / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuaMesh, cuaModules, cuaNumModules, Ha, prng, mc_cone_angledeg, mc_acceptance_rate);
	}

	return mc_acceptance_rate.to_cpu();
}

///////////////////////////////////////////////////////////////
// PARALLEL CONSTRAINED MONTE-CARLO METROPOLIS

//Take a constrained Monte Carlo Metropolis step in this atomistic mesh
cuBReal Atom_Mesh_CubicCUDA::Iterate_MonteCarloCUDA_Constrained(cuBReal mc_cone_angledeg)
{
	Zero_MCAux_Atom_Mesh_CubicCUDA <<< 1, CUDATHREADS >>> (mc_acceptance_rate);

	//TO DO

	return mc_acceptance_rate.to_cpu();
}


#endif

#endif