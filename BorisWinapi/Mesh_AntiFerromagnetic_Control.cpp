#include "stdafx.h"
#include "Mesh_AntiFerromagnetic.h"

//----------------------------------- FERROMAGNETIC MESH QUANTITIES CONTROL 

void AFMesh::SetMagnetisationAngle(double polar, double azim, Rect rectangle)
{
#if COMPILECUDA == 1

	if (pMeshCUDA) {

		if (rectangle.IsNull()) {

			pMeshCUDA->M()->setnonempty(Polar_to_Cartesian(cuReal3(Ms, polar, azim)));
			pMeshCUDA->M2()->setnonempty(Polar_to_Cartesian(cuReal3(-1.0 * Ms, polar, azim)));
		}
		else {

			pMeshCUDA->M()->setrectnonempty((cuRect)rectangle, Polar_to_Cartesian(cuReal3(Ms, polar, azim)));
			pMeshCUDA->M2()->setrectnonempty((cuRect)rectangle, Polar_to_Cartesian(cuReal3(-1.0 * Ms, polar, azim)));
		}

		return;
	}

#endif

	if (M.linear_size()) {

		if (rectangle.IsNull()) {

			M.setnonempty(Polar_to_Cartesian(DBL3(Ms, polar, azim)));
			M2.setnonempty(Polar_to_Cartesian(DBL3(-1.0 * Ms, polar, azim)));
		}
		else {

			M.setrectnonempty(rectangle, Polar_to_Cartesian(DBL3(Ms, polar, azim)));
			M2.setrectnonempty(rectangle, Polar_to_Cartesian(DBL3(-1.0 * Ms, polar, azim)));
		}
	}
}

//Invert magnetisation direction in given mesh (must be ferromagnetic)
void AFMesh::SetInvertedMagnetisation(void)
{
#if COMPILECUDA == 1
	//refresh M from gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

#pragma omp parallel for
	for (int idx = 0; idx < M.linear_size(); idx++) {

		if (M.is_not_empty(idx)) {

			M[idx] *= -1;
			M2[idx] *= -1;
		}
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_from_cpuvec(M);
		pMeshCUDA->M2()->copy_from_cpuvec(M2);
	}
#endif
}

//Set random magentisation distribution in given mesh (must be ferromagnetic)
void AFMesh::SetRandomMagnetisation(void)
{
#if COMPILECUDA == 1
	//refresh M from gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

	BorisRand prng(GetTickCount());

#pragma omp parallel for
	for (int idx = 0; idx < M.linear_size(); idx++) {

		if (M.is_not_empty(idx)) {

			double theta = prng.rand() * TWO_PI;
			double phi = prng.rand() * TWO_PI;

			M[idx] = M[idx].norm() * DBL3(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
			M2[idx] = -1.0 * M[idx];
		}
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_from_cpuvec(M);
		pMeshCUDA->M2()->copy_from_cpuvec(M2);
	}
#endif
}

void AFMesh::SetMagnetisationDomainWall(int longitudinal, int transverse, double width, double position)
{
#if COMPILECUDA == 1
	//refresh M from gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

	if (!M.linear_size()) return;

	//1: x, 2: y, 3: z, 1: -x, 2: -y, 3: -z
	if (longitudinal == transverse) return;
	if (width <= 0) return;

	auto set_dbl3_component = [](int selector, double component) -> DBL3 {

		DBL3 value = DBL3();

		if (abs(selector) == 1) value = DBL3(component, 0, 0);
		else if (abs(selector) == 2) value = DBL3(0, component, 0);
		else if (abs(selector) == 3) value = DBL3(0, 0, component);

		return value;
	};

	int width_cells = (int)floor_epsilon(width / h.x);
	int start_cell = int(position / h.x);

	for (int k = 0; k < n.z; k++) {
#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = start_cell; i < start_cell + width_cells && i < n.x; i++) {

				if (i < 0) continue;

				int cell_index = i + j * n.x + k * n.x*n.y;

				if (M.is_not_empty(cell_index)) {

					double magnitude = GetMagnitude(M[cell_index]);

					double longitudinal_value = -get_sign(longitudinal) * magnitude * tanh(12 * ((double)(i - start_cell) / ((double)width_cells - 1) - 0.5));
					double transverse_value = get_sign(transverse) * magnitude / cosh(12 * ((double)(i - start_cell) / ((double)width_cells - 1) - 0.5));

					M[cell_index] = set_dbl3_component(longitudinal, longitudinal_value) + set_dbl3_component(transverse, transverse_value);
					M2[cell_index] = -1.0 * M[cell_index];
				}
			}
		}
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_from_cpuvec(M);
		pMeshCUDA->M2()->copy_from_cpuvec(M2);
	}
#endif
}

//set Neel skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
void AFMesh::SetSkyrmion(int orientation, int chirality, Rect skyrmion_rect)
{
#if COMPILECUDA == 1
	//refresh M from gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

	if (!M.linear_size()) return;

	INT3 start = round(skyrmion_rect.s / h);
	INT3 end = round(skyrmion_rect.e / h);

	if (chirality != 1) chirality = -1;
	if (orientation != 1) orientation = -1;

	//skyrmion center xy coordinates
	DBL2 center_position = DBL2(skyrmion_rect.s.x + skyrmion_rect.size().x / 2, skyrmion_rect.s.y + skyrmion_rect.size().y / 2);

	//skyrmion radius : skyrmion_rect must be a square in the xy plane
	double radius = (skyrmion_rect.e.x - skyrmion_rect.s.x) / 2;

	for (int i = (start.i >= 0 ? start.i : 0); i < end.i && i < n.x; i++) {
		for (int j = (start.j >= 0 ? start.j : 0); j < end.j && j < n.y; j++) {
			for (int k = (start.k >= 0 ? start.k : 0); k < end.k && k < n.z; k++) {

				int idx = i + j * n.x + k * n.x*n.y;

				double magnitude = GetMagnitude(M[idx]);

				//distance from origin
				DBL2 position = DBL2(((double)i + 0.5) * h.x, ((double)j + 0.5) * h.y);
				double distance = (center_position - position).norm();

				if (M.is_not_empty(idx) && distance <= radius) {

					if (distance >= radius * (1 - SKYRMION_RING_WIDTH)) {

						//skyrmion ring
						double z_val = -tanh((distance - radius * (1 - SKYRMION_RING_WIDTH / 2)) * SKYRMION_TANH_CONST / (radius * SKYRMION_RING_WIDTH / 2));

						//build unit vector towards center
						DBL2 radial = DBL2();
						if (center_position != position) radial = chirality * orientation * (center_position - position).normalized();

						//now set the length of the vector towards centre so the 3D vector with z_val z component is a unit vector
						radial *= sqrt(1 - z_val * z_val);

						M[idx] = orientation * magnitude * DBL3(radial.x, radial.y, z_val);
						M2[idx] = -1.0 * M[idx];
					}
					else {

						//inside the skyrmion
						M[idx] = orientation * magnitude * DBL3(0, 0, 1);
						M2[idx] = -1.0 * M[idx];
					}
				}
			}
		}
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_from_cpuvec(M);
		pMeshCUDA->M2()->copy_from_cpuvec(M2);
	}
#endif
}

//set Bloch skyrmion with given chirality (outside is up: 1, outside is down: -1) in given rectangle (relative to mesh), calculated in the x-y plane
void AFMesh::SetSkyrmionBloch(int orientation, int chirality, Rect skyrmion_rect)
{
#if COMPILECUDA == 1
	//refresh M from gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

	if (!M.linear_size()) return;

	INT3 start = round(skyrmion_rect.s / h);
	INT3 end = round(skyrmion_rect.e / h);

	if (chirality != 1) chirality = -1;
	if (orientation != 1) orientation = -1;

	//skyrmion center xy coordinates
	DBL2 center_position = DBL2(skyrmion_rect.s.x + skyrmion_rect.size().x / 2, skyrmion_rect.s.y + skyrmion_rect.size().y / 2);

	//skyrmion radius : skyrmion_rect must be a square in the xy plane
	double radius = (skyrmion_rect.e.x - skyrmion_rect.s.x) / 2;

	for (int i = (start.i >= 0 ? start.i : 0); i < end.i && i < n.x; i++) {
		for (int j = (start.j >= 0 ? start.j : 0); j < end.j && j < n.y; j++) {
			for (int k = (start.k >= 0 ? start.k : 0); k < end.k && k < n.z; k++) {

				int idx = i + j * n.x + k * n.x*n.y;

				double magnitude = GetMagnitude(M[idx]);

				//distance from origin
				DBL2 position = DBL2(((double)i + 0.5) * h.x, ((double)j + 0.5) * h.y);
				double distance = (center_position - position).norm();

				if (M.is_not_empty(idx) && distance <= radius) {

					if (distance >= radius * (1 - SKYRMION_RING_WIDTH)) {

						//skyrmion ring
						double z_val = -tanh((distance - radius * (1 - SKYRMION_RING_WIDTH / 2)) * SKYRMION_TANH_CONST / (radius * SKYRMION_RING_WIDTH / 2));

						//build unit vector towards center
						DBL2 radial = DBL2();
						if (center_position != position) radial = chirality * (center_position - position).normalized();

						//build unit vector perpendicular to both xy_vec and z direction:
						DBL3 radial_perp = DBL3(0, 0, 1) ^ DBL3(radial.x, radial.y, 0);

						//now set the length of the radial_perp vector so the 3D vector with z_val z component is a unit vector
						radial_perp *= sqrt(1 - z_val * z_val);

						M[idx] = orientation * magnitude * DBL3(radial_perp.x, radial_perp.y, z_val);
						M2[idx] = -1.0 * M[idx];
					}
					else {

						//inside the skyrmion
						M[idx] = orientation * magnitude * DBL3(0, 0, 1);
						M2[idx] = -1.0 * M[idx];
					}
				}
			}
		}
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_from_cpuvec(M);
		pMeshCUDA->M2()->copy_from_cpuvec(M2);
	}
#endif
}

//set M from given data VEC (0 values mean empty points) -> stretch data to M dimensions if needed.
void AFMesh::SetMagnetisationFromData(VEC<DBL3>& data)
{
#if COMPILECUDA == 1
	//refresh M from gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_to_cpuvec(M);
		pMeshCUDA->M2()->copy_to_cpuvec(M2);
	}
#endif

	if (!M.linear_size()) return;

	//make sure dimensions match
	if (!data.resize(M.n)) return;

	//copy values in data, as well as shape, inverting it for M2
	M.copy_values(data);
	M2.copy_values(data);

#pragma omp parallel for
	for (int idx = 0; idx < M2.linear_size(); idx++) {

		if (M2.is_not_empty(idx)) {

			M2[idx] *= -1;
		}
	}

#if COMPILECUDA == 1
	//refresh gpu memory
	if (pMeshCUDA) {

		pMeshCUDA->M()->copy_from_cpuvec(M);
		pMeshCUDA->M2()->copy_from_cpuvec(M2);
	}
#endif
}

//set periodic boundary conditions for magnetization
BError AFMesh::Set_PBC_X(int pbc_x)
{
	BError error(__FUNCTION__);
	
	//set pbc conditions for demag module if available
	if (IsModuleSet(MOD_DEMAG)) {

		INT3 pbc_images = reinterpret_cast<Demag*>(pMod(MOD_DEMAG))->Get_PBC();

		pbc_images.x = pbc_x;

		error = reinterpret_cast<Demag*>(pMod(MOD_DEMAG))->Set_PBC(pbc_images);
	}
	else return error(BERROR_GPUERROR_CRIT);
	
	return error;
}

//set periodic boundary conditions for magnetization
BError AFMesh::Set_PBC_Y(int pbc_y)
{
	BError error(__FUNCTION__);
	
	//set pbc conditions for demag module if available
	if (IsModuleSet(MOD_DEMAG)) {

		INT3 pbc_images = reinterpret_cast<Demag*>(pMod(MOD_DEMAG))->Get_PBC();

		pbc_images.y = pbc_y;

		error = reinterpret_cast<Demag*>(pMod(MOD_DEMAG))->Set_PBC(pbc_images);
	}
	else return error(BERROR_GPUERROR_CRIT);
	
	return error;
}

//set periodic boundary conditions for magnetization
BError AFMesh::Set_PBC_Z(int pbc_z)
{
	BError error(__FUNCTION__);
	
	//set pbc conditions for demag module if available
	if (IsModuleSet(MOD_DEMAG)) {

		INT3 pbc_images = reinterpret_cast<Demag*>(pMod(MOD_DEMAG))->Get_PBC();

		pbc_images.z = pbc_z;

		error = reinterpret_cast<Demag*>(pMod(MOD_DEMAG))->Set_PBC(pbc_images);
	}
	else return error(BERROR_GPUERROR_CRIT);
	
	return error;
}