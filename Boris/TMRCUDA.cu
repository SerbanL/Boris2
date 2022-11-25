#include "TMRCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TMR

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "ManagedAtom_MeshCUDA.h"

//--------------------------------------------------------------- TMR computation Launcher

__global__ void CalculateElectricalConductivity_TMR_COS_Kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal>* pRAtmr_p_equation,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal>* pRAtmr_ap_equation,
	ManagedMeshCUDA* pMeshFM_Top, size_t pMeshFM_Top_size, ManagedMeshCUDA* pMeshFM_Bot, size_t pMeshFM_Bot_size,
	ManagedAtom_MeshCUDA* pMeshAtom_Top, size_t pMeshAtom_Top_size, ManagedAtom_MeshCUDA* pMeshAtom_Bot, size_t pMeshAtom_Bot_size)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuBReal>& V = *cuMesh.pV;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n_e = elC.n;
	cuReal3 h_e = elC.h;
	cuRect meshRect = elC.rect;

	if (idx < n_e.x * n_e.y) {

		//expecting elC to be 2D
		int i = idx % n_e.x;
		int j = idx / n_e.x;

		int idx = i + j * n_e.x;

		//skip empty cells
		if (elC.is_not_empty(idx)) {

			//m direction values for top and bottom, used to calculate TMR in this cell
			cuReal3 m_top = cuReal3(), m_bot = cuReal3();

			//TOP

			//Look at FM meshes first, if any
			for (int mesh_idx = 0; mesh_idx < pMeshFM_Top_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Top = *(pMeshFM_Top[mesh_idx].pM);
				cuRect tmeshRect = M_Top.rect;

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h_e.x + meshRect.s.x - tmeshRect.s.x,
					(j + 0.5) * h_e.y + meshRect.s.y - tmeshRect.s.y,
					M_Top.h.z / 2);

				if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

				m_top = cu_normalize(M_Top.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Top.h.z)));
				break;
			}

			if (m_top == cuReal3()) {

				//Look at Atomistic meshes, if any
				for (int mesh_idx = 0; mesh_idx < pMeshAtom_Top_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M_Top = *(pMeshAtom_Top[mesh_idx].pM1);
					cuRect tmeshRect = M_Top.rect;

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h_e.x + meshRect.s.x - tmeshRect.s.x,
						(j + 0.5) * h_e.y + meshRect.s.y - tmeshRect.s.y,
						M_Top.h.z / 2);

					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

					m_top = cu_normalize(M_Top.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Top.h.z)));
					break;
				}
			}

			//BOTTOM

			//Look at FM meshes first, if any
			for (int mesh_idx = 0; mesh_idx < pMeshFM_Bot_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMeshFM_Bot[mesh_idx].pM);
				cuRect bmeshRect = M_Bot.rect;

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h_e.x + meshRect.s.x - bmeshRect.s.x,
					(j + 0.5) * h_e.y + meshRect.s.y - bmeshRect.s.y,
					bmeshRect.e.z - bmeshRect.s.z - (M_Bot.h.z / 2));

				//can't couple to an empty cell
				if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

				m_bot = cu_normalize(M_Bot.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Bot.h.z)));
				break;
			}

			if (m_bot == cuReal3()) {

				//Look at Atomistic meshes, if any
				for (int mesh_idx = 0; mesh_idx < pMeshAtom_Bot_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M_Bot = *(pMeshAtom_Bot[mesh_idx].pM1);
					cuRect bmeshRect = M_Bot.rect;

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h_e.x + meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * h_e.y + meshRect.s.y - bmeshRect.s.y,
						bmeshRect.e.z - bmeshRect.s.z - (M_Bot.h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

					m_bot = cu_normalize(M_Bot.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Bot.h.z)));
					break;
				}
			}

			//now apply TMR formula to store conductivity value

			//cos dependence of RA product (where dRA = RAap - RAp):
			//RA = (RAp + dRA * (1 - m1.m2)/2)
			//so resistivity is (thickness t):
			//ro = (RAp + dRA * (1 - m1.m2)/2) / t
			//so set conductivity as: 1 / ro

			cuBReal RAtmr_p = *cuMesh.pRAtmr_p;
			cuBReal RAtmr_ap = *cuMesh.pRAtmr_ap;
			cuBReal elecCond = *cuMesh.pelecCond;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pRAtmr_p, RAtmr_p, *cuMesh.pRAtmr_ap, RAtmr_ap, *cuMesh.pelecCond, elecCond);

			//RA bias dependence if set
			if (pRAtmr_p_equation || pRAtmr_ap_equation) {

				cuBReal bias = 0.0;

				if (V.linear_size()) {

					cuBReal Vt_1 = V[idx + (n_e.z - 1) * n_e.x * n_e.y];
					cuBReal Vt_2 = V[idx + (n_e.z - 2) * n_e.x * n_e.y];

					cuBReal Vb_1 = V[idx];
					cuBReal Vb_2 = V[idx + n_e.x * n_e.y];

					cuBReal Vt = 1.5 * Vt_1 - 0.5 * Vt_2;
					cuBReal Vb = 1.5 * Vb_1 - 0.5 * Vb_2;

					bias = Vt - Vb;
				}

				cuBReal RAtmr_p0 = RAtmr_p;
				cuBReal RAtmr_ap0 = RAtmr_ap;

				if (pRAtmr_p_equation) RAtmr_p = pRAtmr_p_equation->evaluate(RAtmr_p0, RAtmr_ap0, bias);
				if (pRAtmr_ap_equation) RAtmr_ap = pRAtmr_ap_equation->evaluate(RAtmr_p0, RAtmr_ap0, bias);
			}

			for (int k = 0; k < n_e.k; k++) {

				if (elecCond > 0.0) {

					//Metallic pinholes
					elC[idx + k * n_e.x * n_e.y] = elecCond;
				}
				else {

					//TMR
					elC[idx + k * n_e.x * n_e.y] = meshRect.height() / (RAtmr_p + (RAtmr_ap - RAtmr_p) * (1 - m_top * m_bot) / 2);
				}
			}
		}
	}
}

__global__ void CalculateElectricalConductivity_TMR_SLONCZEWSKI_Kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal>* pRAtmr_p_equation,
	ManagedFunctionCUDA<cuBReal, cuBReal, cuBReal>* pRAtmr_ap_equation,
	ManagedMeshCUDA* pMeshFM_Top, size_t pMeshFM_Top_size, ManagedMeshCUDA* pMeshFM_Bot, size_t pMeshFM_Bot_size,
	ManagedAtom_MeshCUDA* pMeshAtom_Top, size_t pMeshAtom_Top_size, ManagedAtom_MeshCUDA* pMeshAtom_Bot, size_t pMeshAtom_Bot_size)
{
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;
	cuVEC_VC<cuBReal>& V = *cuMesh.pV;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n_e = elC.n;
	cuReal3 h_e = elC.h;
	cuRect meshRect = elC.rect;

	if (idx < n_e.x * n_e.y) {

		//expecting elC to be 2D
		int i = idx % n_e.x;
		int j = idx / n_e.x;

		int idx = i + j * n_e.x;

		//skip empty cells
		if (elC.is_not_empty(idx)) {

			//m direction values for top and bottom, used to calculate TMR in this cell
			cuReal3 m_top = cuReal3(), m_bot = cuReal3();

			//TOP

			//Look at FM meshes first, if any
			for (int mesh_idx = 0; mesh_idx < pMeshFM_Top_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Top = *(pMeshFM_Top[mesh_idx].pM);
				cuRect tmeshRect = M_Top.rect;

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h_e.x + meshRect.s.x - tmeshRect.s.x,
					(j + 0.5) * h_e.y + meshRect.s.y - tmeshRect.s.y,
					M_Top.h.z / 2);

				if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

				m_top = cu_normalize(M_Top.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Top.h.z)));
				break;
			}

			if (m_top == cuReal3()) {

				//Look at Atomistic meshes, if any
				for (int mesh_idx = 0; mesh_idx < pMeshAtom_Top_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M_Top = *(pMeshAtom_Top[mesh_idx].pM1);
					cuRect tmeshRect = M_Top.rect;

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h_e.x + meshRect.s.x - tmeshRect.s.x,
						(j + 0.5) * h_e.y + meshRect.s.y - tmeshRect.s.y,
						M_Top.h.z / 2);

					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s)) continue;

					m_top = cu_normalize(M_Top.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Top.h.z)));
					break;
				}
			}

			//BOTTOM

			//Look at FM meshes first, if any
			for (int mesh_idx = 0; mesh_idx < pMeshFM_Bot_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMeshFM_Bot[mesh_idx].pM);
				cuRect bmeshRect = M_Bot.rect;

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h_e.x + meshRect.s.x - bmeshRect.s.x,
					(j + 0.5) * h_e.y + meshRect.s.y - bmeshRect.s.y,
					bmeshRect.e.z - bmeshRect.s.z - (M_Bot.h.z / 2));

				//can't couple to an empty cell
				if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

				m_bot = cu_normalize(M_Bot.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Bot.h.z)));
				break;
			}

			if (m_bot == cuReal3()) {

				//Look at Atomistic meshes, if any
				for (int mesh_idx = 0; mesh_idx < pMeshAtom_Bot_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M_Bot = *(pMeshAtom_Bot[mesh_idx].pM1);
					cuRect bmeshRect = M_Bot.rect;

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h_e.x + meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * h_e.y + meshRect.s.y - bmeshRect.s.y,
						bmeshRect.e.z - bmeshRect.s.z - (M_Bot.h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s)) continue;

					m_bot = cu_normalize(M_Bot.weighted_average(cell_rel_pos, cuReal3(h_e.x, h_e.y, M_Bot.h.z)));
					break;
				}
			}

			//now apply TMR formula to store conductivity value

			//Slonczewski form : cos dependence of conductivity
			//RA = 2*RAp / [ (1 + RAp/RAap) + (1 - RAp/RAap)cos(theta) ]

			cuBReal RAtmr_p = *cuMesh.pRAtmr_p;
			cuBReal RAtmr_ap = *cuMesh.pRAtmr_ap;
			cuBReal elecCond = *cuMesh.pelecCond;
			cuMesh.update_parameters_ecoarse(idx, *cuMesh.pRAtmr_p, RAtmr_p, *cuMesh.pRAtmr_ap, RAtmr_ap, *cuMesh.pelecCond, elecCond);

			//RA bias dependence if set
			if (pRAtmr_p_equation || pRAtmr_ap_equation) {

				cuBReal bias = 0.0;

				if (V.linear_size()) {

					cuBReal Vt_1 = V[idx + (n_e.z - 1) * n_e.x * n_e.y];
					cuBReal Vt_2 = V[idx + (n_e.z - 2) * n_e.x * n_e.y];

					cuBReal Vb_1 = V[idx];
					cuBReal Vb_2 = V[idx + n_e.x * n_e.y];

					cuBReal Vt = 1.5 * Vt_1 - 0.5 * Vt_2;
					cuBReal Vb = 1.5 * Vb_1 - 0.5 * Vb_2;

					bias = Vt - Vb;
				}

				cuBReal RAtmr_p0 = RAtmr_p;
				cuBReal RAtmr_ap0 = RAtmr_ap;

				if (pRAtmr_p_equation) RAtmr_p = pRAtmr_p_equation->evaluate(RAtmr_p0, RAtmr_ap0, bias);
				if (pRAtmr_ap_equation) RAtmr_ap = pRAtmr_ap_equation->evaluate(RAtmr_p0, RAtmr_ap0, bias);
			}

			for (int k = 0; k < n_e.k; k++) {

				if (elecCond > 0.0) {

					//Metallic pinholes
					elC[idx + k * n_e.x * n_e.y] = elecCond;
				}
				else {

					//TMR
					elC[idx + k * n_e.x * n_e.y] = meshRect.height() * ((1 + RAtmr_p / RAtmr_ap) + (1 - RAtmr_p / RAtmr_ap) * m_top * m_bot) / (2 * RAtmr_p);
				}
			}
		}
	}
}

//calculate electrical conductivity based on TMR formula
void TMRCUDA::CalculateElectricalConductivity_TMR(TMR_ TMR_type)
{
	switch (TMR_type) {
		
	case TMR_COS:
		CalculateElectricalConductivity_TMR_COS_Kernel <<< (pMeshCUDA->n_e.x * pMeshCUDA->n_e.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh,
			(RAtmr_p_equation.is_set() ? &RAtmr_p_equation.get_x() : nullptr), (RAtmr_ap_equation.is_set() ? &RAtmr_ap_equation.get_x() : nullptr),
			pMeshFM_Top, pMeshFM_Top.size(), pMeshFM_Bot, pMeshFM_Bot.size(),
			pMeshAtom_Top, pMeshAtom_Top.size(), pMeshAtom_Bot, pMeshAtom_Bot.size());
		break;
		
	case TMR_SLONCZEWSKI:
		CalculateElectricalConductivity_TMR_SLONCZEWSKI_Kernel <<< (pMeshCUDA->n_e.x * pMeshCUDA->n_e.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
			(pMeshCUDA->cuMesh,
			(RAtmr_p_equation.is_set() ? &RAtmr_p_equation.get_x() : nullptr), (RAtmr_ap_equation.is_set() ? &RAtmr_ap_equation.get_x() : nullptr),
			pMeshFM_Top, pMeshFM_Top.size(), pMeshFM_Bot, pMeshFM_Bot.size(),
			pMeshAtom_Top, pMeshAtom_Top.size(), pMeshAtom_Bot, pMeshAtom_Bot.size());
		break;
	}
}

//-------------------Calculation Methods : Electric Field

__global__ void CalculateElectricField_TMR_Kernel(cuVEC<cuReal3>& E, cuVEC_VC<cuBReal>& V)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < V.linear_size()) {

		//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
		if (V.is_not_empty(idx)) {

			E[idx] = -1.0 * V.grad_diri(idx);
		}
		else E[idx] = cuReal3(0.0);
	}
}

//calculate electric field as the negative gradient of V
void TMRCUDA::CalculateElectricField(bool open_potential)
{
	CalculateElectricField_TMR_Kernel <<< (pMeshCUDA->n_e.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->E, pMeshCUDA->V);
}

//-------------------Others

BError TMRCUDA::SetBiasEquation_Parallel(const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	BError error(CLASS_STR(TMRCUDA));
	
	if (!fspec.size()) RAtmr_p_equation.clear();
	else if (!RAtmr_p_equation.make_scalar(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

BError TMRCUDA::SetBiasEquation_AntiParallel(const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	BError error(CLASS_STR(TMRCUDA));
	
	if (!fspec.size()) RAtmr_ap_equation.clear();
	else if (!RAtmr_ap_equation.make_scalar(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif

#endif