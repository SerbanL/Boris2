#pragma once

//TO DO:
//pbc
//cmbnd
//mesh transfer

//////////////// mcuObject POLICY CLASS for managed cuVEC

#include "mcuObject.h"

#include "cuTypes.h"
#include "alloc_cpy.h"
#include "mGPU_transfer.h"

#include "mcuValue.h"
#include "mcuArray.h"
#include "cuVEC_VC.h"

//VType : VEC value type, e.g. float, cuFLT3, etc.
//MType : managed type, e.g. cuVEC, or cuVEC_VC
//Some methods need to be handled differently for cuVEC and cuVEC_VC, so use the std::enable_if_t mechanism
template <typename VType, typename MType>
class mcuVEC
{
private:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// SPECIAL OBJECTS FROM CONSTRUCTOR (all Policy classes)

	//reference to mcu_obj manager for which this is a policy class
	mcu_obj<MType, mcuVEC<VType, MType>>& mng;

	//multi-GPU configuration (list of physical devices configured with memory transfer type configuration set)
	mGPUConfig& mGPU;

	//////////////////////////////////////////////////////////////////////////////////
	//
	// POLICY CLASS DATA (specific)

	//----------------- (cuVEC)

	//dimensions for each device
	cuSZ3* pn_d = nullptr;

	//rectangle for each device (in absolute coordinates)
	cuRect* prect_d = nullptr;

	//box for each device (box relative to entire cuVEC, i.e. cell start and end box coordinates make sense in cuBox(n), where n dimensions of entire cuVEC)
	cuBox* pbox_d = nullptr;

	//----------------- (cuVEC) line profiles

	//used to extract profiles
	mcu_arr<VType> profile_aux;
	mcu_arr<cuReal2> profile_component_aux;

	//----------------- (cuVEC) Histograms

	//gpu memory allocated on base gpu for histogram data
	//size of vector is number devices, each holding an auxiliary cu_arr with size histogram_size (see cuVEC) allocated on base gpu only. used to add histograms on base gpu.
	std::vector<cu_arr<cuBReal>> histogram_base_aux;
	//cu_arr on base gpu storing collection of cu_arrs above (so we can pass them to a cuda kernel)
	cu_arr<cuBReal*> histogram_base_aux_col;

	//store handles to histogram data from multiple devices
	//base handles store histogram (see cuVEC). additional handles defined for base gpu only, storing handles to cu_arr held in histogram_base_aux
	mGPU_Transfer<cuBReal> histogram_transf;

	//----------------- (cuVEC_VC) PBCs
	
	//Periodic boundary conditions for evaluating differential operators.
	//pbcs are recorded here and implemented through appropriate halo exchanges, if pbc is set along the halo stacking direction.
	//for other pbc directions just set flags in cuVEC_VC. Note that halo flags take precedence over pbc flags, so setting pbc flags in all devices is fine.
	int pbc_x = 0;
	int pbc_y = 0;
	int pbc_z = 0;

	//----------------- (cuVEC_VC) Halo data

	//Halo exchanges : array of cu_arr, set to number of devices

	//phalo_..._n : used for a given mesh for transfers on n side.
	//phalo_..._p : used for a given mesh for transfers on p side.
	//Thus e.g. : if a given mesh (index idx) needs transfer on its n side, then adjacent mesh needs to place values in the linear memory phalo_..._n[idx - 1], which will then be transfered to phalo_..._n[idx], etc.

	//spaces for halo exchanges : ngbrFlags, negative and positive sides. 
	std::vector<cu_arr<int>*> phalo_ngbr_n;
	std::vector<cu_arr<int>*> phalo_ngbr_p;

	//halo depth - allowed values 0 or 1, but will be extended later to > 1.
	int halo_depth = 1;

	//in which direction are the devices stacked? this determines where the halos are placed.
	//e.g. halo_flag = NF2_HALOX, NF2_HALOY or NF2_HALOZ, and halo_flag_n, halo_flag_p are respective values for na dn p sides.
	int halo_flag_n = NF2_HALONX, halo_flag_p = NF2_HALOPX, halo_flag = NF2_HALOX;

	//----------------- (cuVEC_VC) Halo transfers

	//objects used to handle memory transfers between devices far halos
	mGPU_Transfer<int> halo_ngbr_n_transf;
	mGPU_Transfer<int> halo_ngbr_p_transf;

	mGPU_Transfer<VType> halo_quant_n_transf;
	mGPU_Transfer<VType> halo_quant_p_transf;

	//-----------------

public:

	//----------------- (cuVEC)

	//overall dimensions along x, y and z of the quantity
	cuSZ3 n = cuSZ3(0);

	//cellsize of structured mesh (same for all devices)
	cuReal3 h = cuReal3(0);

	//overall rectangle, same units as h. VEC has n number of cells, so n * h gives the rect dimensions
	cuRect rect = cuRect();

private:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// AUXILIARY (specific)

	//--------------------------------------------AUXILIARY : mcuVEC_halo.h

	//allocate memory for halos
	
	//not applicable to cuVEC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC<VType_>>::value>* = nullptr>
	void allocate_halos(void) {}

	//for cuVEC_VC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void allocate_halos(void);

	//coordinate ngbr flag exchanges to set halo conditions in managed devices
	
	//not applicable to cuVEC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC<VType_>>::value>* = nullptr>
	void set_halo_conditions(void) {}

	//for cuVEC_VC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void set_halo_conditions(void);
	
	//--------------------------------------------AUXILIARY : mcuVEC_histo.cuh, mcuVEC_hist.h

	void setup_histogram_transfers(void);

	//add arrays in histogram_base_aux_col together, setting result in first one
	void add_histograms(void);

	//--------------------------------------------AUXILIARY : mcuVEC_aux.h

	//for given h and rect values, find n
	cuSZ3 get_n_from_h_and_rect(const cuReal3& h_, const cuRect& rect_);

	//for configured number of devices and given n value, calculate the number of cells (as cuSZ3) each device should be assigned
	//this is returned as first: n value for all devices but last. second: n value for last device (this can differ if devices cannot be assigned all same n value)
	std::pair<cuSZ3, cuSZ3> get_devices_n_values(cuSZ3 new_n);

	//--------------------------------------------AUXILIARY : mcuVEC_flags.h

	//set pbc from VEC_VC
	template <typename cpuVEC_VC, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void set_pbc_from(cpuVEC_VC& vec_vc);

	//pbcs not defined for VEC
	template <typename cpuVEC, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC<VType_>>::value>* = nullptr>
	void set_pbc_from(cpuVEC& vec) {}

	//set pbc in VEC_VC
	template <typename cpuVEC_VC, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void set_pbc_to(cpuVEC_VC& vec_vc);

	//pbcs not defined for VEC
	template <typename cpuVEC, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC<VType_>>::value>* = nullptr>
	void set_pbc_to(cpuVEC& vec) {}

	//////////////////////////////////////////////////////////////////////////////////
	//
	// AUXILIARY (all Policy classes)

	//clear all allocated memory
	void clear_memory_aux(void);

public:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// CONSTRUCTORS (all Policy classes)

	//--------------------------------------------CONSTRUCTORS : mcuVEC_mng.h

	//constructor for this policy class : void
	mcuVEC<VType, MType>(mcu_obj<MType, mcuVEC<VType, MType>>& mng_, mGPUConfig& mGPU_);

	void construct_policy(void) {}

	//construct to given number of cells. called from mng after void constructor finished.
	void construct_policy(const cuSZ3& n_) { resize(n_); }

	//construct to given dimensions. called from mng after void constructor finished.
	void construct_policy(const cuReal3& h_, const cuRect& rect_) { resize(h_, rect_); }

	//construct to given dimensions and initialize to given value. called from mng after void constructor finished.
	void construct_policy(const cuReal3& h_, const cuRect& rect_, VType value) { assign(h_, rect_, value); }

	//assignment operator
	mcuVEC<VType, MType>& operator=(const mcuVEC<VType, MType>& copyThis);

	//destructor
	virtual ~mcuVEC<VType, MType>() { clear_memory_aux(); }

	//////////////////////////////////////////////////////////////////////////////////
	//
	// POLICY CLASS METHODS (specific)

	//--------------------------------------------Halos : mcuVEC_halo.h

	//exchange values in all halos
	
	//not applicable to cuVEC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC<VType_>>::value>* = nullptr>
	void exchange_halos(void) {}

	//for cuVEC_VC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void exchange_halos(void);

	//--------------------------------------------Others

	//number of points allocated for a given GPU (linear index, not physical device index)
	size_t device_size(int idx) { return pn_d[idx].dim(); }
	
	//smallest allocated number of points from all devices
	size_t min_device_size(void)
	{
		size_t min_size = n.dim();
		for (int idx = 0; idx < mGPU.get_num_devices(); idx++) if (device_size(idx) < min_size) min_size = device_size(idx);
		return min_size;
	}

	//////////////////////////////////////////////////////////////////////////////////
	//"Overload" cuVEC, cuVEC_VC public __host__ methods. Should have the exact order and names and parameters as methods appear in cuVEC / cuVEC_VC.

	//--------------------------------------------COPY TO / FROM VEC : mcuVEC_mng.h

	//copy everything from a VEC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC>
	bool set_from_cpuvec(cpuVEC& vec);

	//copy everything to a VEC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC>
	bool set_cpuvec(cpuVEC& vec);

	//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC
	template <typename cpuVEC>
	bool copy_from_cpuvec(cpuVEC& vec);

	//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC
	template <typename cpuVEC>
	bool copy_to_cpuvec(cpuVEC& vec);

	//load_cuarr, store_cuarr : not available, consider doing versions with mcu_arr

	//copy values from a std::vector (cpu memory)
	template <typename SType>
	bool copy_from_vector(std::vector<SType>& vec);

	//copy values to a std::vector (cpu memory)
	template <typename SType>
	bool copy_to_vector(std::vector<SType>& vec);

	//copy flags only from vec_vc, where sizes must match. Applicable to cuVEC_VC only
	template <typename cpuVEC_VC, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool copyflags_from_cpuvec(cpuVEC_VC& vec_vc);

	//--------------------------------------------COPY TO ANOTHER cuVEC : mcuVEC_mng.h

	//extract values from this and place them in cuvec : both must have same rectangle, but can differ in cuVEC<VType>::h - cuvec.h <= this->cuVEC<VType>::h needed (and hence cuVEC<VType>::n); e.g. this method allows extraction of a coarser cuvec.
	void extract_cuvec(mcu_obj<MType, mcuVEC<VType, MType>>& mcuvec);

	//--------------------------------------------EXTRACT A LINE PROFILE : mcuVEC_extract.h

	//no wrap-around is currently used if multiple devices are available (TO DO). for a single device, wrap-around is used as defined in cuVEC methods.
	//cordinates are relative

	//1. Profile values only, without stencil operation

	//extract profile to a mcu_arr : extract size points starting at start in the direction step for the given number of points (size); use weighted average to extract profile with h stencil only
	//profile_gpu resized as needed for each device
	void extract_profilevalues(size_t size, mcu_arr<VType>& profile_gpu, cuReal3 start, cuReal3 step);

	//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
	void extract_profilevalues_component_x(size_t size, mcu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
	void extract_profilevalues_component_y(size_t size, mcu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
	void extract_profilevalues_component_z(size_t size, mcu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

	//2. Profile values only, with stencil operation around profile point

	//extract profile components: extract starting at start in the direction end - step, with given step; use weighted average to extract profile with given stencil
	//all coordinates are relative positions. Return profile values in profile_gpu (which is resized as needed for each device)
	bool extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<VType>& profile_gpu);

	//as above, but only store profile in internal memory (line_profile) so we can read it out later as needed
	bool extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil);

	//extract profile components: extract starting at start in the direction end - step, with given step; use weighted average to extract profile with given stencil
	//all coordinates are relative positions. Return profile values in profile_cpu.
	template <typename SType>
	bool extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<SType>& profile_cpu);

	//3. Profile values for individual components together with profile position returned in cuReal2/DBL2, with stencil operation around profile point (capped to mesh size)

	//as above but only component x and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu (which is resized as needed for each device)
	bool extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu);
	//as above but only component x and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
	bool extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

	//as above but only component y and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu (which is resized as needed for each device)
	bool extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu);
	//as above but only component y and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
	bool extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

	//as above but only component z and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu (which is resized as needed for each device)
	bool extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu);
	//as above but only component z and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
	bool extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

	//as above but only component which has largest value for the first point and pack in profile position too (after stencil averaging) (for VAL3 floating types only). Return data in profile_gpu (which is resized as needed for each device)
	bool extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu);
	//as above but only component which has largest value for the first point and pack in profile position too (after stencil averaging) (for VAL3 floating types only). Return data in profile_cpu. profile_cpu resized as needed.
	bool extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

	//--------------------------------------------HISTOGRAMS : mcuVEC_histo.h

	//compute magnitude histogram data
	//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first.
	//if num_bins not given then use default value of 100
	//if macrocell_dims greater than 1 in any dimension then first average mesh data in macrocells of given size
	//without macrocell then pass in num_nonempty_cells (number of nonempty cells); if this is not passed it's counted first (costs another kernel launch)
	//output transferred to cpu in histogram_cpu
	bool get_mag_histogram(
		std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu,
		int num_bins, double& min, double& max, size_t num_nonempty_cells = 0, cuINT3 macrocell_dims = cuINT3(1));

	//get angular deviation histogram. deviation from ndir direction is calculated, or the average direction if ndir not specified (IsNull)
	bool get_ang_histogram(
		std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu,
		int num_bins, double& min, double& max, size_t num_nonempty_cells = 0, cuINT3 macrocell_dims = cuINT3(1), VType ndir = VType());

	//--------------------------------------------ITERATORS

	//Not needed at Policy level

	//--------------------------------------------SIZING : mcuVEC_mng.h

	//resize number of cells, breaking up space amongst devices
	bool resize(cuSZ3 new_n);
	
	//resize and set shape using linked vec. Applicable to cuVEC_VC only
	template <typename LVType, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool resize(cuSZ3 new_n, mcu_obj<cuVEC_VC<LVType>, mcuVEC<LVType, cuVEC_VC<LVType>>>& linked_vec);

	//set dimensions and assign value, breaking up space amongst devices
	bool resize(cuReal3 new_h, cuRect new_rect);

	//resize and set shape using linked vec. Applicable to cuVEC_VC only
	template <typename LVType, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool resize(cuReal3 new_h, cuRect new_rect, mcu_obj<cuVEC_VC<LVType>, mcuVEC<LVType, cuVEC_VC<LVType>>>& linked_vec);

	//resize number of cells and assign value, breaking up space amongst devices
	bool assign(cuSZ3 new_n, VType value);
	
	//set value and shape from linked vec
	template <typename LVType, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool assign(cuSZ3 new_n, VType value, mcu_obj<cuVEC_VC<LVType>, mcuVEC<LVType, cuVEC_VC<LVType>>>& linked_vec);

	//set dimensions and assign value, breaking up space amongst devices
	bool assign(cuReal3 new_h, cuRect new_rect, VType value);

	//set value and shape from linked vec
	template <typename LVType, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool assign(cuReal3 new_h, cuRect new_rect, VType value, mcu_obj<cuVEC_VC<LVType>, mcuVEC<LVType, cuVEC_VC<LVType>>>& linked_vec);

	//set everything to zero but h
	void clear(void);

	//set rect start (i.e. shift the entire rectangle to align with given absolute starting coordinates
	void set_rect_start(const cuReal3& rect_start);
	void shift_rect_start(const cuReal3& shift);

	//--------------------------------------------FLAG CHECKING : mcuVEC_flags.h

	int get_nonempty_cells_cpu(void);

	//--------------------------------------------SET CELL FLAGS - EXTERNAL USE : mcuVEC_flags.h

	//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
	//return false on memory allocation failure only, otherwise return true even if surface_rect was not valid. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool set_dirichlet_conditions(cuRect surface_rect, VType value);

	//clear all dirichlet flags and vectors. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void clear_dirichlet_flags(void);

	//set pbc conditions : setting any to false clears flags. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void set_pbc(int pbc_x_, int pbc_y_, int pbc_z_);

	//not applicable to cuVEC
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC<VType_>>::value>* = nullptr>
	void set_pbc(int pbc_x_, int pbc_y_, int pbc_z_) {}
	
	//clear all pbc flags : can also be achieved setting all flags to false in set_pbc but this one is more readable. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void clear_pbc(void);

	//clear all composite media boundary flags. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void clear_cmbnd_flags(void);

	//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations (if status true, else clear the skip cells flags in this rectangle). Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void set_skipcells(cuRect rectangle, bool status = true);

	//clear all skip cell flags. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void clear_skipcells(void);

	//Robin conditions. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void set_robin_conditions(cuReal2 robin_v_, cuReal2 robin_px_, cuReal2 robin_nx_, cuReal2 robin_py_, cuReal2 robin_ny_, cuReal2 robin_pz_, cuReal2 robin_nz_);

	//clear all Robin boundary conditions and values. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void clear_robin_conditions(void);

	//--------------------------------------------MULTIPLE ENTRIES SETTERS : mcuVEC_oper.h

	//set value in box
	void setbox(cuBox box, VType value = VType());

	//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this cuVEC's rectangle.
	void setrect(const cuRect& rectangle, VType value = VType());

	//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void delrect(cuRect rectangle);

	//set value in all cells
	void set(VType value);

	//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	void renormalize(PType new_norm);

	//mask values in cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0). Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	bool apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth = 0);

	//exactly the same as assign value - do not use assign as it is slow (sets flags). Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void setnonempty(VType value = VType());

	//set value in non-empty cells only in given rectangle (relative coordinates). Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void setrectnonempty(const cuRect& rectangle, VType value = VType());

	//shift all the values in this cuVEC by the given delta (units same as cuVEC<VType>::h). Shift values in given shift_rect (absolute coordinates). Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void shift_x(cuBReal delta, cuRect shift_rect);

	//--------------------------------------------VEC GENERATORS : mcuVEC_generate.h

	//linear : use interpolation to set values in this VEC based on projected distance between position1 and position2 and given fixed end values.
	bool generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, VType value1, cuRect contact2, VType value2);

	//similar to generate_linear except new dimensions not set
	//also allow 'degeneracy' : multiple linear generators may be superimposed with use of degeneracy : degeneracy.first is index, degeneracy.second is number of geenerators
	//if using degeneracy make sure these are called in order, or at least index 0 goes first
	void set_linear(cuRect contact1, VType value1, cuRect contact2, VType value2, cuReal2 degeneracy = cuReal2());

	//--------------------------------------------GETTERS : mcuVEC_aux.h

	size_t linear_size_cpu(void)  const { return n.dim(); }

	cuSZ3 size_cpu(void) { return n; }
	cuReal3 cellsize_cpu(void) { return h; }

	//from cell index return cell center coordinates (relative to start of rectangle)
	cuReal3 cellidx_to_position_cpu(int idx);

	//from cell index return cell center coordinates (relative to start of rectangle)
	cuReal3 cellidx_to_position_cpu(cuINT3 ijk);

	//return cell index from relative position : the inverse of cellidx_to_position
	int position_to_cellidx_cpu(const cuReal3& position);

	//get index of cell which contains position (absolute value, not relative to start of rectangle), capped to mesh size
	cuINT3 cellidx_from_position_cpu(cuReal3 absolute_position);

	//get cell rectangle (absolute values, not relative to start of mesh rectangle) for cell with index ijk
	cuRect get_cellrect_cpu(cuINT3 ijk);

	//get_cellrect using single index.
	cuRect get_cellrect_cpu(int idx);

	//extract box of cells intersecting with the given rectangle (rectangle is in absolute coordinates). Cells in box : from and including start, up to but not including end; Limited to cuVEC sizes.
	cuBox box_from_rect_max_cpu(cuRect rectangle);

	//extract box of cells completely included in the given rectangle (rectangle is in absolute coordinates).
	cuBox box_from_rect_min_cpu(cuRect rectangle);

	//count cells which don't have a null value set : i.e. non-empty; set result in aux_integer
	void count_nonempty_cells(void);

	//--------------------------------------------ARITHMETIC OPERATIONS ON ENTIRE VEC : mcuVEC_arith.h

	//add to this vec the values in add_this : must have same size : size
	void add_values(mcu_obj<MType, mcuVEC<VType, MType>>& add_this);

	//subtract from this vec the values in sub_this : must have same size : size
	void sub_values(mcu_obj<MType, mcuVEC<VType, MType>>& sub_this);

	//scale all stored values by the given constant. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void scale_values(cuBReal constant);

	//--------------------------------------------AVERAGING OPERATIONS : mcuVEC_avg.h

	//average in a box (which should be contained in the cuVEC dimensions)
	VType average(cuBox box);

	//average over given rectangle (relative to this cuVEC's rect)
	VType average(cuRect rectangle = cuRect());

	//as above but exclude empty points from averaging
	VType average_nonempty(cuBox box);
	VType average_nonempty(cuRect rectangle = cuRect());

	//just sum over non-empty cells instead of averaging. Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	VType sum_nonempty(cuBox box);

	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	VType sum_nonempty(cuRect rectangle = cuRect());

	//--------------------------------------------NUMERICAL PROPERTIES : mcuVEC_nprops.h

	//Find min and max values. rectangles are relative to this VEC.
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax(cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax(cuRect rectangle = cuRect());

	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax_component_x(cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax_component_x(cuRect rectangle = cuRect());

	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax_component_y(cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax_component_y(cuRect rectangle = cuRect());

	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax_component_z(cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	cuVAL2<PType> get_minmax_component_z(cuRect rectangle = cuRect());

	//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : mcuVEC_cmbnd.h

	//TO DO
	
	//--------------------------------------------LAPLACE / POISSON EQUATION : mcuVEC_solve.h, mcuVEC_solve.cuh

	//LAPLACE

	//Take one SOR iteration for Laplace equation on this VEC. Return error (maximum change in quantity from one iteration to the next) by reference.
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this).
	//Applicable to cuVEC_VC only
	template <typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void IterateLaplace_SOR(mcu_val<cuBReal>& damping, mcu_val<cuBReal>& max_error, mcu_val<cuBReal>& max_val);

	//POISSON with homogeneous Neumann boundaries

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = Poisson_RHS
	//Poisson_RHS must be a member const method of Class_Poisson_RHS taking an index value (the index ranges over this VEC) and returning a cuBReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//obj of type Class_Poisson_RHS must be mcu_obj managed so it is entirely in gpu memory (mcu_obj for multiple devices if available)
	//Return error(maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a mcu_val.
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	//Applicable to cuVEC_VC only
	template <typename Class_Poisson_RHS, typename Policy_Class_Poisson_RHS, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void IteratePoisson_SOR(mcu_obj<Class_Poisson_RHS, Policy_Class_Poisson_RHS>& obj, mcu_val<cuBReal>& damping, mcu_val<cuBReal>& max_error, mcu_val<cuBReal>& max_val);

	//POISSON with non-homogeneous Neumann boundaries

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = Poisson_RHS
	//Poisson_RHS must be a member const method of Class_Poisson_NNeu taking an index value (the index ranges over this VEC) and returning a cuBReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//Class_Poisson_NNeu must also define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index) - this is the non-homogeneous Neumann boundary condition at that cell
	//obj of type Class_Poisson_NNeu must be mcu_obj managed so it is entirely in gpu memory (mcu_obj for multiple devices if available).
	//Return error(maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a mcu_val.
	//Dirichlet boundary conditions used, defaulting to non-homogeneous Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	//Applicable to cuVEC_VC only
	template <typename Class_Poisson_NNeu, typename Policy_Class_Poisson_NNeu, typename VType_ = VType, typename MType_ = MType, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>* = nullptr>
	void IteratePoisson_NNeu_SOR(mcu_obj<Class_Poisson_NNeu, Policy_Class_Poisson_NNeu>& obj, mcu_val<cuBReal>& damping, mcu_val<cuBReal>& max_error, mcu_val<cuBReal>& max_val);

	//--------------------------------------------MESH TRANSFER : mcuVEC_MeshTransfer.h

	//TO DO
};

//Macros to simplify declarations, e.g. mcu_VEC(cuReal3) defines a cuVEC<cuReal3>, mcu_obj managed across multiple GPUs etc.
#define mcu_VEC(type) mcu_obj<cuVEC<type>, mcuVEC<type, cuVEC<type>>>
#define mcu_VEC_VC(type) mcu_obj<cuVEC_VC<type>, mcuVEC<type, cuVEC_VC<type>>>
