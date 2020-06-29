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
				DBL3 Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(idx);

				pMesh->Heff[idx] += Hexch;

				energy += pMesh->M[idx] * Hexch;
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

//-------------------Energy density methods

double Exch_6ngbr_Neu::GetEnergyDensity(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return dynamic_cast<Exch_6ngbr_NeuCUDA*>(pModuleCUDA)->GetEnergyDensity(avRect);
#endif

	double energy = 0;
	int num_points = 0;

	INT3 n = pMesh->n;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy, num_points)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			//only average over values in given rectangle
			if (!avRect.contains(pMesh->M.cellidx_to_position(idx))) continue;

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double A = pMesh->A;
				pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->Ms, Ms);

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_neu evaluates to zero in the CMBND coupling direction.
				DBL3 Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(idx);

				energy += pMesh->M[idx] * Hexch;
				num_points++;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy, num_points) 
		for (int idx = 0; idx < n.dim(); idx++) {

			//only average over values in given rectangle
			if (!avRect.contains(pMesh->M.cellidx_to_position(idx))) continue;

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

				energy += (pMesh->M[idx] * Hexch + pMesh->M2[idx] * Hexch2) / 2;
				num_points++;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// FINAL ENERGY DENSITY VALUE //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//average energy density, this is correct, see notes - e_ex ~= -(mu0/2) M.H as can be derived from the full expression for e_ex. It is not -mu0 M.H, which is the energy density for a dipole in a field !
	if (num_points) energy *= -MU0 / (2 * num_points);
	else energy = 0;

	return energy;
}

double Exch_6ngbr_Neu::GetEnergy_Max(Rect& rectangle)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return dynamic_cast<Exch_6ngbr_NeuCUDA*>(pModuleCUDA)->GetEnergy_Max(rectangle);
#endif

	INT3 n = pMesh->n;

	OmpReduction<double> emax;
	emax.new_minmax_reduction();

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for
		for (int idx = 0; idx < n.dim(); idx++) {

			//only obtain max in given rectangle
			if (!rectangle.contains(pMesh->M.cellidx_to_position(idx))) continue;

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double A = pMesh->A;
				pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->Ms, Ms);

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_neu evaluates to zero in the CMBND coupling direction.
				DBL3 Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(idx);

				emax.reduce_max(fabs(pMesh->M[idx] * Hexch));
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) 
	{

#pragma omp parallel for
		for (int idx = 0; idx < n.dim(); idx++) {

			//only obtain max in given rectangle
			if (!rectangle.contains(pMesh->M.cellidx_to_position(idx))) continue;

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

				emax.reduce_max(fabs((pMesh->M[idx] * Hexch + pMesh->M2[idx] * Hexch2) / 2));
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// FINAL ENERGY DENSITY VALUE //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	return emax.maximum() * MU0 / 2;
}

//Compute exchange energy density and store it in displayVEC
void Exch_6ngbr_Neu::Compute_Exchange(VEC<double>& displayVEC)
{
	displayVEC.resize(pMesh->h, pMesh->meshRect);

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		dynamic_cast<Exch_6ngbr_NeuCUDA*>(pModuleCUDA)->Compute_Exchange(displayVEC);
		return;
	}
#endif

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double A = pMesh->A;
				pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->Ms, Ms);

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_neu evaluates to zero in the CMBND coupling direction.
				DBL3 Hexch = (2 * A / (MU0*Ms*Ms)) * pMesh->M.delsq_neu(idx);

				displayVEC[idx] = -(MU0 / 2) * (pMesh->M[idx] * Hexch);
			}
			else displayVEC[idx] = 0.0;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for
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

				displayVEC[idx] = -(MU0 / 2) * (pMesh->M[idx] * Hexch + pMesh->M2[idx] * Hexch2) / 2;
			}
			else displayVEC[idx] = 0.0;
		}
	}
}

#endif
