#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include <vector>
#include <tuple>

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

#include "TransportCUDA.h"
#include "STransportCUDA_GInterf_V.h"
#include "STransportCUDA_GInterf_S.h"

class SuperMesh;
class STransport;

using namespace std;

class STransportCUDA :
	public ModulesCUDA
{

	friend STransport;

private:

	//pointer to cpu version of this super-mesh module
	STransport * pSTrans;

	SuperMesh* pSMesh;

	//used to monitor Poisson equation convergence : max_error is the maximum change in solved quantity from one iteration to the next, max_value is the maximum value in that iteration
	//Find global max_value to divide and obtain normalized maximum error
	cu_obj<cuReal> max_error, max_value;

	//---------------------- CMBND data

	//CMBND contacts for all contacting transport meshes - these are ordered by first vector index; for each transport mesh there could be multiple contacting meshes and these are ordered by second vector index
	//CMBNDInfo describes the contact between 2 meshes, allowing calculation of values at cmbnd cells based on continuity of a potential and flux
	//we need a cu_obj managed copy of CMBNDcontacts from STransport so we can pass it to cuda kernels efficiently
	vector< vector< cu_obj<CMBNDInfoCUDA> > > CMBNDcontactsCUDA;
	//...and we also need a cpu-memory version, even though we can access it using pSTrans - the problem is, we need the cpu data in .cu files where we cannot define STransport (as nvcc will then attempt to compile BorisLib)
	vector< vector<CMBNDInfoCUDA> > CMBNDcontacts;

	//list of all transport modules in transport meshes (same ordering as first vector in CMBNDcontacts)
	vector<TransportCUDA*> pTransport;

	//vector of pointers to all V - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	vector<cu_obj<cuVEC_VC<cuReal>>*> pV;

	//vector of pointers to all S - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	vector<cu_obj<cuVEC_VC<cuReal3>>*> pS;

	//Object with functions for calculating cmbnd values when using interface conductances
	cu_obj<STransportCUDA_GInterf_V_Funcs> gInterf_V;
	cu_obj<STransportCUDA_GInterf_S_NF_Funcs> gInterf_S_NF;
	cu_obj<STransportCUDA_GInterf_S_FN_Funcs> gInterf_S_FN;

	//fixed SOR damping to use for V and S (second value) Poisson equations - copied from SOR_damping in STransport
	cu_obj<cuReal> SOR_damping_V;
	cu_obj<cuReal> SOR_damping_S;

private:

	void Zero_Errors(void);

	//scale all potential values in all V cuVECs by given scaling value
	void scale_potential_values(cuReal scaling);

	//set potential values using a slope between the potential values of ground and another electrode (if set)
	void initialize_potential_values(void);

	//-----Charge Transport only

	//solve for V and Jc in all meshes using SOR
	void solve_charge_transport_sor(void);

	//calculate and set values at composite media boundaries after all other cells have been computed and set
	void set_cmbnd_charge_transport(void);

	//-----Spin and Charge Transport

	//solve for V, Jc and S in all meshes using SOR for Poisson equation and FTCS for S equation
	void solve_spin_transport_sor(void);

	//calculate and set values at composite media boundaries for V (when using spin transport solver)
	void set_cmbnd_spin_transport_V(void);

	//calculate and set values at composite media boundaries for S
	void set_cmbnd_spin_transport_S(void);

	//Calculate interface spin accumulation torque (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(void);

public:

	STransportCUDA(SuperMesh* pSMesh_, STransport* pSTrans_);
	~STransportCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------Getters

	//return interfacial spin torque in given mesh with matching transport module
	cu_obj<cuVEC<cuReal3>>& GetInterfacialSpinTorque(TransportCUDA* pMeshTrans);
};

#else

class STransportCUDA
{
};

#endif

#endif
