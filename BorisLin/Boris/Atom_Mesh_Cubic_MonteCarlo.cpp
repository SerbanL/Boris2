#include "stdafx.h"
#include "Atom_Mesh_Cubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

//defines std::execution::par_unseq used for parallel std::sort
//#include <execution>

#include "Atom_MeshParamsControl.h"
#include "SuperMesh.h"

//Take a Monte Carlo step in this atomistic mesh
void Atom_Mesh_Cubic::Iterate_MonteCarlo(double acceptance_rate)
{
	if (mc_disabled) return;

	if (mc_constrain) {

		if (mc_parallel) Iterate_MonteCarlo_Parallel_Constrained();
		else Iterate_MonteCarlo_Serial_Constrained();
	}
	else {

		if (mc_parallel) Iterate_MonteCarlo_Parallel_Classic();
		else Iterate_MonteCarlo_Serial_Classic();
	}

	///////////////////////////////////////////////////////////////
	// ADAPTIVE CONE ANGLE

	MonteCarlo_AdaptiveAngle(mc_cone_angledeg, acceptance_rate);
}

#if COMPILECUDA == 1
//Take a Monte Carlo step in this atomistic mesh
void Atom_Mesh_Cubic::Iterate_MonteCarloCUDA(double acceptance_rate)
{ 
	if (mc_disabled) return;

	if (paMeshCUDA && mc_parallel) {

		if (mc_constrain) mc_acceptance_rate = reinterpret_cast<Atom_Mesh_CubicCUDA*>(paMeshCUDA)->Iterate_MonteCarloCUDA_Constrained(mc_cone_angledeg, acceptance_rate);
		else mc_acceptance_rate = reinterpret_cast<Atom_Mesh_CubicCUDA*>(paMeshCUDA)->Iterate_MonteCarloCUDA_Classic(mc_cone_angledeg, acceptance_rate);
	}

	///////////////////////////////////////////////////////////////
	// ADAPTIVE CONE ANGLE

	MonteCarlo_AdaptiveAngle(mc_cone_angledeg, acceptance_rate);
}
#endif

//Take a Monte Carlo step in this atomistic mesh : these functions implement the actual algorithms
void Atom_Mesh_Cubic::Iterate_MonteCarlo_Serial_Classic(void)
{
	//number of moves in this step : one per each spin
	int num_moves = M1.get_nonempty_cells();

	//number of cells in this mesh
	unsigned N = n.dim();

	//make sure indices array has correct memory allocated
	if (mc_indices.size() != N) if (!malloc_vector(mc_indices, N)) return;

	//recalculate it
	mc_acceptance_rate = 0.0;

	//reset indices array so we can shuffle them
#pragma omp parallel for
	for (int idx = 0; idx < N; idx++)
		mc_indices[idx] = idx;

	///////////////////////////////////////////////////////////////
	// SERIAL MONTE-CARLO METROPOLIS

	//make sure number of calls to rand doesn't match the prng period or we're asking for trouble with MC.
	prng.check_periodicity();

	for (int idx = N - 1; idx >= 0; idx--) {

		//pick spin index in a random order, but guarantee each spin is given the chance to move exactly once per MCM step.
		//Note this approach is more efficient than picking spins completely at random, where each spin is moved once per step only on average (tests I've ran show this approach requires more steps to thermalize).
		//Moreover the Metropolis et al. paper also moves each particle exactly once (if accepted) per step. The random picking approach is a little easier to write and requires less memory but overall is not as good.
		//Technical note: prefer using the method below to shuffle spins rather than the std::shuffle algorithm as the latter is slower.
		int rand_idx = floor(prng.rand() * (idx + 1));
		//idx < N check to be sure - it is theoretically possible idx = N can be generated but the probability is small.
		if (rand_idx >= N) rand_idx = N - 1;
		int spin_idx = mc_indices[rand_idx];
		mc_indices[rand_idx] = mc_indices[idx];

		//only consider non-empty and non-frozen cells
		if (M1.is_not_empty(spin_idx) && !M1.is_skipcell(spin_idx)) {

			//Picked spin is M1[spin_idx]
			DBL3 M1_old = M1[spin_idx];

			//obtain rotated spin in a cone around the picked spin
			double theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
			double phi_rot = prng.rand() * 2 * PI;
			//Move spin in cone with uniform random probability distribution. This approach only requires 2 random numbers to be generated. 
			//Also using a Gaussian distribution to move spin around the initial spin is less efficient, requiring more steps to thermalize.
			DBL3 M1_new = relrotate_polar(M1_old, theta_rot, phi_rot);

			//find energy change : new - old
			double energy_delta = 0.0;
			for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

				energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx, M1_new);
			}

			//Compute acceptance probability
			double P_accept = 0.0, P = 1.0;
			if (base_temperature > 0.0) {

				P_accept = exp(-energy_delta / (BOLTZMANN * base_temperature));
				//uniform random number between 0 and 1
				P = prng.rand();
			}
			else if (energy_delta < 0) P_accept = 1.0;

			if (P <= P_accept) {

				mc_acceptance_rate += 1.0 / num_moves;

				//renormalize spin to mu_s to avoid floating point error creep
				double mu_s_val = mu_s;
				update_parameters_mcoarse(spin_idx, mu_s, mu_s_val);
				M1_new.renormalize(mu_s_val);

				//set new spin
				M1[spin_idx] = M1_new;
			}
		}
	}
}

void Atom_Mesh_Cubic::Iterate_MonteCarlo_Serial_Constrained(void)
{
	//number of moves in this step : one per each spin
	int num_moves = M1.get_nonempty_cells();

	//number of cells in this mesh
	unsigned N = n.dim();

	//make sure indices array has correct memory allocated
	if (mc_indices.size() != N) if (!malloc_vector(mc_indices, N)) return;

	//recalculate it
	mc_acceptance_rate = 0.0;

	//Total moment along CMC direction
	double cmc_M = 0.0;

	//reset indices array so we can shuffle them
#pragma omp parallel for reduction(+:cmc_M)
	for (int idx = 0; idx < N; idx++) {

		mc_indices[idx] = idx;
		if (M1.is_not_empty(idx)) cmc_M += M1[idx] * cmc_n;
	}

	///////////////////////////////////////////////////////////////
	// SERIAL CONSTRAINED MONTE-CARLO METROPOLIS

	// Following Asselin et al., PRB 82, 054415 (2010) but with some differences: 
	// 1) adaptive cone angle for target acceptance of 0.5
	// 2) move spins in a cone using uniform pdf polar and azimuthal angles
	// 3) pick spin pairs exactly once in a random order per CMC step

	//make sure number of calls to rand doesn't match the prng period or we're asking for trouble with MC.
	prng.check_periodicity();

	for (int idx = N - 1; idx >= 0; idx -= 2) {

		if (idx == 0) break;

		int rand_idx1 = floor(prng.rand() * (idx + 1));
		//idx < N check to be sure - it is theoretically possible idx = N can be generated but the probability is small.
		if (rand_idx1 >= N) rand_idx1 = N - 1;
		int spin_idx1 = mc_indices[rand_idx1];
		mc_indices[rand_idx1] = mc_indices[idx];

		int rand_idx2 = floor(prng.rand() * idx);
		int spin_idx2 = mc_indices[rand_idx2];
		mc_indices[rand_idx2] = mc_indices[idx - 1];

		if (M1.is_not_empty(spin_idx1) && M1.is_not_empty(spin_idx2) && !M1.is_skipcell(spin_idx1) && !M1.is_skipcell(spin_idx2)) {

			//Picked spins are M1[spin_idx1], M1[spin_idx2]
			DBL3 M_old1 = M1[spin_idx1];
			DBL3 M_old2 = M1[spin_idx2];

			//rotate to a system with x axis along cmc_n : easier to apply algorithm formulas this way
			//NOTE : don't rotate x axis towards cmc_n, but rotate cmc_n towards x axis, i.e. use invrotate_polar, not rotate_polar
			DBL3 Mrot_old1 = invrotate_polar(M_old1, cmc_n);
			DBL3 Mrot_old2 = invrotate_polar(M_old2, cmc_n);

			//obtain rotated spin in a cone around the first picked spin
			double theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
			double phi_rot = prng.rand() * 2 * PI;
			DBL3 Mrot_new1 = relrotate_polar(Mrot_old1, theta_rot, phi_rot);

			//adjust second spin to keep required total moment direction
			DBL3 Mrot_new2 = DBL3(0.0, Mrot_old2.y + Mrot_old1.y - Mrot_new1.y, Mrot_old2.z + Mrot_old1.z - Mrot_new1.z);
			double sq2 = Mrot_new2.y * Mrot_new2.y + Mrot_new2.z * Mrot_new2.z;
			double sqnorm = M_old2.norm()*M_old2.norm();

			if (sq2 < sqnorm) {

				Mrot_new2.x = get_sign(Mrot_old2.x) * sqrt(sqnorm - sq2);

				//Obtain new spins in original coordinate system, i.e. rotate back
				DBL3 M_new1 = rotate_polar(Mrot_new1, cmc_n);
				DBL3 M_new2 = rotate_polar(Mrot_new2, cmc_n);

				//find energy change : new - old
				double energy_delta = 0.0;
				for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

					energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx1, M_new1) + pMod[mod_idx]->Get_EnergyChange(spin_idx2, M_new2);
				}

				double cmc_M_new = cmc_M + Mrot_new1.x + Mrot_new2.x - Mrot_old1.x - Mrot_old2.x;

				if (cmc_M_new > 0.0) {

					//Compute acceptance probability
					double P_accept = 0.0, P = 1.0;
					if (base_temperature > 0.0) {

						P_accept = (cmc_M_new / cmc_M) * (cmc_M_new / cmc_M) * (abs(Mrot_old2.x) / abs(Mrot_new2.x)) * exp(-energy_delta / (BOLTZMANN * base_temperature));
						//uniform random number between 0 and 1
						P = prng.rand();
					}
					else if (energy_delta < 0) P_accept = 1.0;

					if (P <= P_accept) {

						//move accepted (x2 since we moved 2 spins)
						mc_acceptance_rate += 2.0 / num_moves;

						//turns out this line isn't really necessary (either cmc_M ~= cmc_M_new at low temperatures, or cmc_M varies about a steady average value so the effect of just using cmc_M computed at the start of the step averages out)
						//For the serial algorithm keep it, but for the parallel algorithm keeping this line is very problematic (there are solutions, e.g. atomic writes, but at the cost of performance).
						cmc_M = cmc_M_new;

						//renormalize spin to mu_s to avoid floating point error creep
						double mu_s_val = mu_s;
						update_parameters_mcoarse(spin_idx1, mu_s, mu_s_val);
						M_new1.renormalize(mu_s_val);

						update_parameters_mcoarse(spin_idx2, mu_s, mu_s_val);
						M_new2.renormalize(mu_s_val);

						//Set new spins
						M1[spin_idx1] = M_new1;
						M1[spin_idx2] = M_new2;
					}
				}
			}
		}
	}
}

void Atom_Mesh_Cubic::Iterate_MonteCarlo_Parallel_Classic(void)
{
	//number of moves in this step : one per each spin
	int num_moves = M1.get_nonempty_cells();

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
		for (int idx_jk = 0; idx_jk < M1.n.y * M1.n.z; idx_jk++) {

			int j = idx_jk % M1.n.y;
			int k = (idx_jk / M1.n.y) % M1.n.z;

			//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
			bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

			//For red pass (first) i starts from red_nudge. For black pass (second) i starts from !red_nudge.
			for (int i = (1 - rb) * red_nudge + rb * (!red_nudge); i < M1.n.x; i += 2) {

				int spin_idx = i + j * M1.n.x + k * M1.n.x*M1.n.y;

				//only consider non-empty and non-frozen cells
				if (M1.is_not_empty(spin_idx) && !M1.is_skipcell(spin_idx)) {

					//Picked spin is M1[spin_idx]
					DBL3 M1_old = M1[spin_idx];

					//obtain rotated spin in a cone around the picked spin
					double theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
					double phi_rot = prng.rand() * 2 * PI;
					//Move spin in cone with uniform random probability distribution. This approach only requires 2 random numbers to be generated. 
					//Also using a Gaussian distribution to move spin around the initial spin is less efficient, requiring more steps to thermalize.
					DBL3 M1_new = relrotate_polar(M1_old, theta_rot, phi_rot);

					//find energy change : new - old
					double energy_delta = 0.0;
					for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

						energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx, M1_new);
					}

					//Compute acceptance probability
					double P_accept = 0.0, P = 1.0;
					if (base_temperature > 0.0) {

						P_accept = exp(-energy_delta / (BOLTZMANN * base_temperature));
						//uniform random number between 0 and 1
						P = prng.rand();
					}
					else if (energy_delta < 0) P_accept = 1.0;

					if (P <= P_accept) {

						acceptance_rate += 1.0 / num_moves;

						//renormalize spin to mu_s to avoid floating point error creep
						double mu_s_val = mu_s;
						update_parameters_mcoarse(spin_idx, mu_s, mu_s_val);
						M1_new.renormalize(mu_s_val);

						//set new spin
						M1[spin_idx] = M1_new;
					}
				}
			}
		}

		mc_acceptance_rate += acceptance_rate;
		rb++;
	}
}

void Atom_Mesh_Cubic::Iterate_MonteCarlo_Parallel_Constrained(void)
{
	//number of moves in this step : one per each spin
	int num_moves = M1.get_nonempty_cells();

	//number of cells in this mesh
	unsigned N = n.dim();

	//number of red and black squares
	int num_reds = (n.z / 2) * (n.x * n.y / 2) + (n.z - (n.z / 2)) * (n.x * n.y / 2 + (n.x * n.y) % 2);
	int num_blacks = (n.z / 2) * (n.x * n.y / 2 + (n.x * n.y) % 2) + (n.z - (n.z / 2)) * (n.x * n.y / 2);

	//make sure indices array has correct memory allocated
	if (mc_indices_red.size() != num_reds) if (!malloc_vector(mc_indices_red, num_reds)) return;
	if (mc_indices_black.size() != num_blacks) if (!malloc_vector(mc_indices_black, num_blacks)) return;

	//recalculate it
	mc_acceptance_rate = 0.0;

	//Total moment along CMC direction
	double cmc_M = 0.0;

	//calculate cmc_M
#pragma omp parallel for reduction(+:cmc_M)
	for (int idx = 0; idx < N; idx++) {

		if (M1.is_not_empty(idx)) cmc_M += M1[idx] * cmc_n;
	}

	///////////////////////////////////////////////////////////////
	// PARALLEL CONSTRAINED MONTE-CARLO METROPOLIS

	// Following Asselin et al., PRB 82, 054415 (2010) but with some differences: 
	// 1) adaptive cone angle for target acceptance of 0.5
	// 2) move spins in a cone using uniform pdf polar and azimuthal angles
	// 3) pick spin pairs exactly once in a random order per CMC step, and pairs picked both from red or from black squares, not mixed.
	// 4) do not update cmc_M after every move, but only calculate cmc_M at the start of every step (this is problematic for parallel algorithm but turns out not necessary - see comments for the serial version where we do update after every move)

	//red-black : two passes will be taken : first generate red and black indices for these passes, together with random doubles so we can shuffle them using sort-based shuffle
	
	//make sure number of calls to rand doesn't match the prng period or we're asking for trouble with MC.
	prng.check_periodicity();

	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
				bool red_nudge = (((j % 2) == 1 && (k % 2) == 0) || (((j % 2) == 0 && (k % 2) == 1)));

				//true for black square, false for red square
				bool black_square = (i - red_nudge) % 2;

				int spin_idx = i + j * n.x + k * n.x*n.y;

				if (!black_square) {

					//Form red squares index (0, 1, 2 ...) from total index
					int even_rows = 0, odd_rows;

					if (k % 2 == 0) {

						even_rows = (j / 2) * (n.x / 2);
						odd_rows = (j - (j / 2)) * (n.x / 2 + n.x % 2);
					}
					else {

						even_rows = (j / 2) * (n.x / 2 + n.x % 2);
						odd_rows = (j - (j / 2)) * (n.x / 2);
					}

					int even_planes = (k / 2) * (n.x * n.y / 2);
					int odd_planes = (k - (k / 2)) * (n.x * n.y / 2 + (n.x * n.y) % 2);

					int idx_red = (i / 2) + (even_rows + odd_rows) + (even_planes + odd_planes);

					mc_indices_red[idx_red] = { prng.rand(), spin_idx };
				}
				else {

					//Form black squares index (0, 1, 2 ...) from total index
					int even_rows = 0, odd_rows;

					if (k % 2 == 0) {

						even_rows = (j / 2) * (n.x / 2);
						odd_rows = (j - (j / 2)) * (n.x / 2 + n.x % 2);
					}
					else {

						even_rows = (j / 2) * (n.x / 2 + n.x % 2);
						odd_rows = (j - (j / 2)) * (n.x / 2);
					}

					int even_planes = (k / 2) * (n.x * n.y / 2);
					int odd_planes = (k - (k / 2)) * (n.x * n.y / 2 + (n.x * n.y) % 2);

					int idx_black = (i / 2) + (even_rows + odd_rows) + (even_planes + odd_planes);

					mc_indices_black[idx_black] = { prng.rand(), spin_idx };
				}
			}
		}
	}

	//Now do the truffle shuffle!
	//Versions with parallel sort (requries <execution>)
	//std::sort(std::execution::par_unseq, mc_indices_red.begin(), mc_indices_red.end());
	//std::sort(std::execution::par_unseq, mc_indices_black.begin(), mc_indices_black.end());
	//Without parallel sort
	std::sort(mc_indices_red.begin(), mc_indices_red.end());
	std::sort(mc_indices_black.begin(), mc_indices_black.end());

	//Finally run the algorithm in parallel using two passes
	int rb = 0;
	while (rb < 2) {

		double acceptance_rate = 0.0;
		int num = (rb == 0 ? num_reds : num_blacks);

#pragma omp parallel for reduction(+:acceptance_rate)
		for (int idx = num - 1; idx >= 0; idx -= 2) {

			if (idx == 0) continue;

			//use pre-shuffled spins : each will be picked exactly once if even number; if odd, then one spin (random) will be left untouched.
			int spin_idx1, spin_idx2;

			if (rb == 0) {
				
				spin_idx1 = mc_indices_red[idx].second;
				spin_idx2 = mc_indices_red[idx - 1].second;
			}
			else {

				spin_idx1 = mc_indices_black[idx].second;
				spin_idx2 = mc_indices_black[idx - 1].second;
			}

			//If there are empty cells then make sure to only pair non-empty ones
			if (M1.is_not_empty(spin_idx1) && M1.is_not_empty(spin_idx2) && !M1.is_skipcell(spin_idx1) && !M1.is_skipcell(spin_idx2)) {

				//Picked spins are M1[spin_idx1], M1[spin_idx2]
				DBL3 M_old1 = M1[spin_idx1];
				DBL3 M_old2 = M1[spin_idx2];

				//rotate to a system with x axis along cmc_n : easier to apply algorithm formulas this way
				//NOTE : don't rotate x axis towards cmc_n, but rotate cmc_n towards x axis, i.e. use invrotate_polar, not rotate_polar
				DBL3 Mrot_old1 = invrotate_polar(M_old1, cmc_n);
				DBL3 Mrot_old2 = invrotate_polar(M_old2, cmc_n);

				//obtain rotated spin in a cone around the first picked spin
				double theta_rot = prng.rand() * mc_cone_angledeg * PI / 180.0;
				double phi_rot = prng.rand() * 2 * PI;
				DBL3 Mrot_new1 = relrotate_polar(Mrot_old1, theta_rot, phi_rot);

				//adjust second spin to keep required total moment direction
				DBL3 Mrot_new2 = DBL3(0.0, Mrot_old2.y + Mrot_old1.y - Mrot_new1.y, Mrot_old2.z + Mrot_old1.z - Mrot_new1.z);
				double sq2 = Mrot_new2.y * Mrot_new2.y + Mrot_new2.z * Mrot_new2.z;
				double sqnorm = M_old2.norm()*M_old2.norm();

				if (sq2 < sqnorm) {

					Mrot_new2.x = get_sign(Mrot_old2.x) * sqrt(sqnorm - sq2);

					//Obtain new spins in original coordinate system, i.e. rotate back
					DBL3 M_new1 = rotate_polar(Mrot_new1, cmc_n);
					DBL3 M_new2 = rotate_polar(Mrot_new2, cmc_n);

					//find energy change : new - old
					double energy_delta = 0.0;
					for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

						energy_delta += pMod[mod_idx]->Get_EnergyChange(spin_idx1, M_new1) + pMod[mod_idx]->Get_EnergyChange(spin_idx2, M_new2);
					}

					//use abs: since we're not updating cmc_M after every spin it can become negative above the Curie temperature. with the cmc_M_new > 0.0 check this will result in solver getting stuck
					double cmc_M_new = abs(cmc_M) + Mrot_new1.x + Mrot_new2.x - Mrot_old1.x - Mrot_old2.x;

					if (cmc_M_new > 0.0) {

						//Compute acceptance probability; make sure cmc_M is not zero otherwise we'll stop accepting anything and solver gets stuck
						double P_accept = 0.0, P = 1.0;
						if (base_temperature > 0.0) {

							if (cmc_M) P_accept = (cmc_M_new / cmc_M) * (cmc_M_new / cmc_M) * (abs(Mrot_old2.x) / abs(Mrot_new2.x)) * exp(-energy_delta / (BOLTZMANN * base_temperature));
							else P_accept = (abs(Mrot_old2.x) / abs(Mrot_new2.x)) * exp(-energy_delta / (BOLTZMANN * base_temperature));
							//uniform random number between 0 and 1
							P = prng.rand();
						}
						else if (energy_delta < 0) P_accept = 1.0;

						if (P <= P_accept) {

							//move accepted (x2 since we moved 2 spins)
							acceptance_rate += 2.0 / num_moves;

							//renormalize spins to mu_s to avoid floating point error creep
							double mu_s_val = mu_s;
							update_parameters_mcoarse(spin_idx1, mu_s, mu_s_val);
							M_new1.renormalize(mu_s_val);

							update_parameters_mcoarse(spin_idx2, mu_s, mu_s_val);
							M_new2.renormalize(mu_s_val);

							//set new spins
							M1[spin_idx1] = M_new1;
							M1[spin_idx2] = M_new2;
						}
					}
				}
			}
		}

		mc_acceptance_rate += acceptance_rate;
		rb++;
	}
}

#endif