#include "stdafx.h"
#include "Exchange.h"

#ifdef MODULE_COMPILATION_EXCHANGE

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "ExchangeCUDA.h"
#endif

/////////////////////////////////////////////////////////////////
//Exch_Const_6ngbr
//
//Hex_i = 2*A/(MU0*Ms^2) * delsq M_i       (delsq is the Laplace operator)
//
//delsq approximated using the 6 neighbor method (with Neumann boundary condition for the boundary cells): 
//
//delsq M_i = sum(j belongs to Ni) (M_j - M_i)/h^2, where Ni is the set of available neighbors for M_i : cells j that share a face with cell i
//
//Exchange energy density contribution from cell i is given by :
//
//Eex,i = -(MU0/2)*M_i.Hex_i = A/Ms^2 * sum(j belongs to Ni) M_i.(M_i-M_j)/h^2
//
// Using Neumann boundary conditions : at boundary, derivative normal to the boundary is zero. Thus to evaluate Laplacian we use:
//
// On inner cells (i, j, k) : M(i-1,j,k) + M(i+1,j,k) + M(i,j-1,k) + M(i,j+1,k) + M(i,j,k-1) + M(i,j,k+1) - 6 * M(i,j,k). Number of neighbors is 6 in this case.
// On boundary cells (including corners) : replace fictitious cell outside of mesh with value of opposite neighbor along same axis and use same equation as for inner cells.
// If both neighbors are missing along an axis (e.g. a 2D problem) then reduce number of neighbors used by 2 (i.e. decrease dimensionality of inner cell Laplacian equation)

Exch_6ngbr_Neu::Exch_6ngbr_Neu(Mesh *pMesh_) : 
	Modules(),
	ExchangeBase(pMesh_),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Exch_6ngbr_Neu::~Exch_6ngbr_Neu() 
{
}

BError Exch_6ngbr_Neu::Initialize(void) 
{
	BError error(CLASS_STR(Exch_6ngbr_Neu));

	error = ExchangeBase::Initialize();
	
	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_EXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_EXCH), 
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_EXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError Exch_6ngbr_Neu::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Exch_6ngbr_Neu));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Exch_6ngbr_Neu::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Exch_6ngbr_Neu));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Exch_6ngbr_NeuCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Exch_6ngbr_Neu::UpdateField(void) 
{
	double energy = 0;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double A = pMesh->A;
				pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->Ms, Ms);

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_neu evaluates to zero in the CMBND coupling direction.
				DBL3 Hexch;

				if (pMesh->base_temperature > 0.0 && pMesh->T_Curie > 0.0) {

					//for finite temperature simulations the magnetization length may have a spatial variation
					//this will not affect the transverse torque (mxH), but will affect the longitudinal term in the sLLB equation (m.H) and cannot be neglected when close to Tc.

					DBL33 Mg = pMesh->M.grad_neu(idx);
					DBL3 dMdx = Mg.x, dMdy = Mg.y, dMdz = Mg.z;

					double delsq_Msq = 2 * pMesh->M[idx] * (pMesh->M.dxx_neu(idx) + pMesh->M.dyy_neu(idx) + pMesh->M.dzz_neu(idx)) + 2 * (dMdx * dMdx + dMdy * dMdy + dMdz * dMdz);
					double Mnorm = pMesh->M[idx].norm();
					Hexch = (2 * A / (MU0*Ms*Ms)) * (pMesh->M.delsq_neu(idx) - pMesh->M[idx] * delsq_Msq / (2 * Mnorm*Mnorm));
				}
				else {

					//zero temperature simulations : magnetization length could still vary but will only affect mxH term, so not needed for 0K simulations.
					Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(idx);
				}

				pMesh->Heff[idx] += Hexch;

				energy += pMesh->M[idx] * Hexch;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Hexch;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * (pMesh->M[idx] * Hexch) / 2;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				DBL2 A_AFM = pMesh->A_AFM;
				DBL2 Ah = pMesh->Ah;
				DBL2 Anh = pMesh->Anh;
				pMesh->update_parameters_mcoarse(idx, pMesh->A_AFM, A_AFM, pMesh->Ah, Ah, pMesh->Anh, Anh, pMesh->Ms_AFM, Ms_AFM);

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_neu evaluates to zero in the CMBND coupling direction.
				
				DBL3 delsq_M_A = pMesh->M.delsq_neu(idx);
				DBL3 delsq_M_B = pMesh->M2.delsq_neu(idx);

				DBL2 M = DBL2(pMesh->M[idx].norm(), pMesh->M2[idx].norm());

				DBL3 Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (-4 * Ah.i * (pMesh->M[idx] ^ (pMesh->M[idx] ^ pMesh->M2[idx])) / (M.i*M.i) + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
				DBL3 Hexch2 = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (-4 * Ah.j * (pMesh->M2[idx] ^ (pMesh->M2[idx] ^ pMesh->M[idx])) / (M.j*M.j) + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

				pMesh->Heff[idx] += Hexch;
				pMesh->Heff2[idx] += Hexch2;

				energy += (pMesh->M[idx] * Hexch + pMesh->M2[idx] * Hexch2) / 2;

				if (Module_Heff.linear_size()) Module_Heff[idx] = Hexch;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = Hexch2;
				if (Module_energy.linear_size()) Module_energy[idx] = -MU0 * (pMesh->M[idx] * Hexch) / 2;
				if (Module_energy2.linear_size()) Module_energy2[idx] = -MU0 * (pMesh->M2[idx] * Hexch2) / 2;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// COUPLING ACROSS MULTIPLE MESHES ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//if exchange coupling across multiple meshes, this is calculation method to use
	std::function<double(int, int, DBL3, DBL3, DBL3, Mesh&, Mesh&)> calculate_coupling = [&](int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 stencil, DBL3 hshift_primary, Mesh& Mesh_pri, Mesh& Mesh_sec) -> double {

		double energy_ = 0.0;

		double hR = hshift_primary.norm();
		double hRsq = hR * hR;

		//AFM to AFM
		if (Mesh_pri.GetMeshType() == MESH_ANTIFERROMAGNETIC && Mesh_sec.GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			//Both meshes antiferromagnetic : both sub-lattices couple, but do not include anti AF coupling here as it has already been included in the main loop above
			//here we only compute differential operators across a boundary.

			DBL2 Ms_AFM = Mesh_pri.Ms_AFM;
			DBL2 A_AFM = Mesh_pri.A_AFM;
			DBL2 Anh = pMesh->Anh;
			Mesh_pri.update_parameters_mcoarse(cell1_idx, Mesh_pri.A_AFM, A_AFM, Mesh_pri.Ms_AFM, Ms_AFM, pMesh->Anh, Anh);

			DBL3 Hexch, Hexch_B;

			//values at cells -1, 1
			DBL3 M_1 = Mesh_pri.M[cell1_idx];
			DBL3 M_m1 = Mesh_sec.M.weighted_average(relpos_m1, stencil);

			DBL3 M_1_B = Mesh_pri.M2[cell1_idx];
			DBL3 M_m1_B = Mesh_sec.M2.weighted_average(relpos_m1, stencil);

			DBL3 delsq_M_A, delsq_M_B;

			if (cell2_idx < Mesh_pri.n.dim() && Mesh_pri.M.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				DBL3 M_2 = Mesh_pri.M[cell2_idx];
				DBL3 M_2_B = Mesh_pri.M2[cell2_idx];

				delsq_M_A = (M_2 + M_m1 - 2 * M_1) / hRsq;
				delsq_M_B = (M_2_B + M_m1_B - 2 * M_1_B) / hRsq;
			}
			else {

				delsq_M_A = (M_m1 - M_1) / hRsq;
				delsq_M_B = (M_m1_B - M_1_B) / hRsq;
			}

			//set effective field value contribution at cell 1 : direct exchange coupling
			Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (Anh.i / (MU0*Ms_AFM.i*Ms_AFM.j)) * delsq_M_B;
			Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (Anh.j / (MU0*Ms_AFM.i*Ms_AFM.j)) * delsq_M_A;

			Mesh_pri.Heff[cell1_idx] += Hexch;
			Mesh_pri.Heff2[cell1_idx] += Hexch_B;

			energy_ = (M_1 * Hexch + M_1_B * Hexch_B) / 2;
		}

		//FM to FM
		else if (Mesh_pri.GetMeshType() == MESH_FERROMAGNETIC && Mesh_sec.GetMeshType() == MESH_FERROMAGNETIC) {

			//both meshes are ferromagnetic
			//here we only compute differential operators across a boundary.

			double Ms, A;
			Ms = Mesh_pri.Ms;
			A = Mesh_pri.A;
			Mesh_pri.update_parameters_mcoarse(cell1_idx, Mesh_pri.A, A, Mesh_pri.Ms, Ms);

			DBL3 Hexch;

			//values at cells -1, 1
			DBL3 M_1 = Mesh_pri.M[cell1_idx];
			DBL3 M_m1 = Mesh_sec.M.weighted_average(relpos_m1, stencil);

			if (cell2_idx < Mesh_pri.n.dim() && Mesh_pri.M.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				DBL3 M_2 = Mesh_pri.M[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling 
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_2 + M_m1 - 2 * M_1) / hRsq;
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_m1 - M_1) / hRsq;
			}

			Mesh_pri.Heff[cell1_idx] += Hexch;

			energy_ = M_1 * Hexch;
		}

		return energy_;
	};

	//if exchange coupled to other meshes calculate the exchange field at marked cmbnd cells and accumulate energy density contribution
	if (pMesh->GetMeshExchangeCoupling()) CalculateExchangeCoupling(energy, calculate_coupling);

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// FINAL ENERGY DENSITY VALUE //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//average energy density, this is correct, see notes - e_ex ~= -(mu0/2) M.H as can be derived from the full expression for e_ex. It is not -mu0 M.H, which is the energy density for a dipole in a field !
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * (pMesh->M.get_nonempty_cells()));
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

//-------------------Energy methods

//FM mesh
double Exch_6ngbr_Neu::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double Ms = pMesh->Ms;
		double A = pMesh->A;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->A, A, pMesh->Ms, Ms);

		DBL3 Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(spin_index);
		double energy_ = pMesh->M[spin_index] * Hexch;

		if (Mnew != DBL3()) {

			//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
			Mnew.renormalize(pMesh->M[spin_index].norm());

			DBL3 Mold = pMesh->M[spin_index];
			pMesh->M[spin_index] = Mnew;
			Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(spin_index);
			double energynew_ = pMesh->M[spin_index] * Hexch;
			pMesh->M[spin_index] = Mold;

			//do not divide by 2 as we are not double-counting here
			return -MU0 * pMesh->h.dim() * (energynew_ - energy_);
		}
		//If Mnew is null then this method is used to obtain current energy only, not energy change
		else return -MU0 * pMesh->h.dim() * energy_;
	}
	else return 0.0;
}

//AFM mesh
DBL2 Exch_6ngbr_Neu::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL2 A_AFM = pMesh->A_AFM;
		DBL2 Ah = pMesh->Ah;
		DBL2 Anh = pMesh->Anh;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->A_AFM, A_AFM, pMesh->Ah, Ah, pMesh->Anh, Anh, pMesh->Ms_AFM, Ms_AFM);

		auto Get_Energy = [&](void) -> DBL2
		{
			DBL3 delsq_M_A = pMesh->M.delsq_neu(spin_index);
			DBL3 delsq_M_B = pMesh->M2.delsq_neu(spin_index);

			DBL3 Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (4 * Ah.i * pMesh->M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
			DBL3 Hexch2 = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (4 * Ah.j * pMesh->M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

			return DBL2(pMesh->M[spin_index] * Hexch, pMesh->M2[spin_index] * Hexch2);
		};

		DBL2 energy_ = Get_Energy();

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
			Mnew_A.renormalize(pMesh->M[spin_index].norm());
			Mnew_B.renormalize(pMesh->M2[spin_index].norm());

			DBL3 Mold_A = pMesh->M[spin_index];
			DBL3 Mold_B = pMesh->M2[spin_index];

			pMesh->M[spin_index] = Mnew_A;
			pMesh->M2[spin_index] = Mnew_B;
			
			DBL2 energynew_ = Get_Energy();

			pMesh->M[spin_index] = Mold_A;
			pMesh->M2[spin_index] = Mold_B;

			//do not divide by 2 as we are not double-counting here
			return -MU0 * pMesh->h.dim() * (energynew_ - energy_);
		}
		//If Mnew is null then this method is used to obtain current energy only, not energy change
		else return -MU0 * pMesh->h.dim() * energy_;
	}
	else return DBL2();
}

//-------------------Torque methods

DBL3 Exch_6ngbr_Neu::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif
