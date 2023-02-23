#pragma once

#include <vector>

#include "cuVEC.h"

////////////////////////////////////////////////////////////////////////////////////////////////// VEC_VC<VType>
//
// extends VEC with vector calculus operations

/////////////////////////////////////////////////////////////////////

//CMBNDInfo describes a composite media boundary contact between 2 meshes of same type, used to calculate values at CMBND cells using boundary conditions
struct CMBNDInfoCUDA;

template <typename VType>
class cuVEC_VC :
	public cuVEC<VType>
{

//the following are used as masks for ngbrFlags. 32 bits in total (4 bytes for an int)

//neighbor existence masks (+x, -x, +y, -y, +z, -z). Bits 0, 1, 2, 3, 4, 5
#define NF_NPX	1
#define NF_NNX	2
#define NF_NPY	4
#define NF_NNY	8
#define NF_NPZ	16
#define NF_NNZ	32

//Existance of at least one neighbor along a given axis : use as masks (test using &).
#define NF_NGBRX	(NF_NPX + NF_NNX)	//test bits 0 and 1
#define NF_NGBRY	(NF_NPY + NF_NNY)	//test bits 2 and 3
#define NF_NGBRZ	(NF_NPZ + NF_NNZ)   //test bits 4 and 5

//existence of both neighbors along axes x, y, z
//the use this check mask with value and see if the result is the same as the mask, e.g. if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) { both neighbors are present }
#define NF_BOTHX	(NF_NPX + NF_NNX)
#define NF_BOTHY	(NF_NPY + NF_NNY)
#define NF_BOTHZ	(NF_NPZ + NF_NNZ)

//periodic boundary condition along x. Set at x sides only if there is a neighbor present at the other side - this is how we know which side to use, +x or -x : bit 6
#define NF_PBCX	64

//periodic boundary condition along y. Set at y sides only if there is a neighbor present at the other side - this is how we know which side to use, +y or -y : bit 7
#define NF_PBCY	128

//periodic boundary condition along z. Set at z sides only if there is a neighbor present at the other side - this is how we know which side to use, +z or -z : bit 8
#define NF_PBCZ	256

//mask for all pbc flags
#define NF_PBC (NF_PBCX + NF_PBCY + NF_PBCZ)

//mask to check for cell with zero value set : bit 9
#define NF_NOTEMPTY	512

//this is not necessarily an empty cell, but mark them to be skipped during computations for some algorithms (e.g. moving mesh algorithm where the ends of the magnetic mesh must not be updated by the ODE solver) : bit 10
#define NF_SKIPCELL	1024

//composite media boundary cells (used to flag cells where boundary conditions must be applied). These flags are not set using set_ngbrFlags, but must be externally set
//cell on positive x side of boundary : bit 11
#define NF_CMBNDPX	2048
//cell on negative x side of boundary : bit 12
#define NF_CMBNDNX	4096
//cell on positive y side of boundary : bit 13
#define NF_CMBNDPY	8192
//cell on negative y side of boundary : bit 14
#define NF_CMBNDNY	16384
//cell on positive z side of boundary : bit 15
#define NF_CMBNDPZ	32768
//cell on negative z side of boundary : bit 16
#define NF_CMBNDNZ	65536

//mask for all cmbnd flags
#define NF_CMBND	(NF_CMBNDPX + NF_CMBNDNX + NF_CMBNDPY + NF_CMBNDNY + NF_CMBNDPZ + NF_CMBNDNZ)
//masks for cmbnd flags along the x, y, or z axes
#define NF_CMBNDX	(NF_CMBNDPX + NF_CMBNDNX)
#define NF_CMBNDY	(NF_CMBNDPY + NF_CMBNDNY)
#define NF_CMBNDZ	(NF_CMBNDPZ + NF_CMBNDNZ)

//off-axis neighbor at +x, +y, 0z (xy differentials) : bit 17
#define NF_XY_PXPY	131072
//off-axis neighbor at +x, -y, 0z (xy differentials) : bit 18
#define NF_XY_PXNY	262144
//off-axis neighbor at -x, +y, 0z (xy differentials) : bit 19
#define NF_XY_NXPY	524288
//off-axis neighbor at -x, -y, 0z (xy differentials) : bit 20
#define NF_XY_NXNY	1048576
//off-axis neighbor at +x, 0y, +z (xz differentials) : bit 21
#define NF_XZ_PXPZ	2097152
//off-axis neighbor at +x, 0y, -z (xz differentials) : bit 22
#define NF_XZ_PXNZ	4194304
//off-axis neighbor at -x, 0y, +z (xz differentials) : bit 23
#define NF_XZ_NXPZ	8388608
//off-axis neighbor at -x, 0y, -z (xz differentials) : bit 24
#define NF_XZ_NXNZ	16777216
//off-axis neighbor at 0x, +y, +z (yz differentials) : bit 25
#define NF_YZ_PYPZ	33554432
//off-axis neighbor at 0x, +y, -z (yz differentials) : bit 26
#define NF_YZ_PYNZ	67108864
//off-axis neighbor at 0x, -y, +z (yz differentials) : bit 27
#define NF_YZ_NYPZ	134217728
//off-axis neighbor at 0x, -y, -z (yz differentials) : bit 28
#define NF_YZ_NYNZ	268435456

//for off-axis neighbors stencil, indicate if mixed second order differential stencil is available, i.e. at least 2 columns must have at least 2 non-empty cells.
//Better to use 3 additional bits to speed up these checks rather than build it every time from the above bits.

//off-axis stencil available in XY plane : bit 29
#define NF_XY_OASTENCIL	536870912
//off-axis stencil available in XZ plane : bit 30
#define NF_XZ_OASTENCIL	1073741824
//off-axis stencil available in YZ plane : bit 31
#define NF_YZ_OASTENCIL	2147483648

//check for full off-axis stencil in a plane as (ngbrFlags[idx] & NF_XY_FULL) == NF_XY_FULL
#define NF_XY_FULL	(NF_XY_PXPY + NF_XY_PXNY + NF_XY_NXPY + NF_XY_NXNY)
#define NF_XZ_FULL	(NF_XZ_PXPZ + NF_XZ_PXNZ + NF_XZ_NXPZ + NF_XZ_NXNZ)
#define NF_YZ_FULL	(NF_YZ_PYPZ + NF_YZ_PYNZ + NF_YZ_NYPZ + NF_YZ_NYNZ)

//Extended flags

//Robin boundary conditions flags
//cell on positive x side of boundary : bit 0 -> use robin_nx values
#define NF2_ROBINPX	1
//cell on negative x side of boundary : bit 1 -> use robin_px values, etc.
#define NF2_ROBINNX	2
//cell on positive y side of boundary : bit 2
#define NF2_ROBINPY	4
//cell on negative y side of boundary : bit 3
#define NF2_ROBINNY	8
//cell on positive z side of boundary : bit 4
#define NF2_ROBINPZ	16
//cell on negative z side of boundary : bit 5
#define NF2_ROBINNZ	32
//flag Robin boundary with a void cell (use robin_v values) : bit 6
#define NF2_ROBINV	64

//mask for all Robin flags
#define NF2_ROBIN	(NF2_ROBINPX + NF2_ROBINNX + NF2_ROBINPY + NF2_ROBINNY + NF2_ROBINPZ + NF2_ROBINNZ + NF2_ROBINV)
//masks for Robin flags along the x, y, or z axes
#define NF2_ROBINX	(NF2_ROBINPX + NF2_ROBINNX)
#define NF2_ROBINY	(NF2_ROBINPY + NF2_ROBINNY)
#define NF2_ROBINZ	(NF2_ROBINPZ + NF2_ROBINNZ)

//these are used in conjunction with dirichlet vectors to indicate dirichlet boundary conditions should be used
//cell on +x side of boundary : bit 7
#define NF2_DIRICHLETPX	128
//cell on -x side of boundary : bit 8
#define NF2_DIRICHLETNX	256
//cell on +y side of boundary : bit 9
#define NF2_DIRICHLETPY	512
//cell on -y side of boundary : bit 10
#define NF2_DIRICHLETNY	1024
//cell on +z side of boundary : bit 11
#define NF2_DIRICHLETPZ	2048
//cell on -z side of boundary : bit 12
#define NF2_DIRICHLETNZ	4096

//masks for x, y, z directions for Dirichlet cells
#define NF2_DIRICHLETX (NF2_DIRICHLETPX + NF2_DIRICHLETNX)
#define NF2_DIRICHLETY (NF2_DIRICHLETPY + NF2_DIRICHLETNY)
#define NF2_DIRICHLETZ (NF2_DIRICHLETPZ + NF2_DIRICHLETNZ)
#define NF2_DIRICHLET (NF2_DIRICHLETX + NF2_DIRICHLETY + NF2_DIRICHLETZ)

//halo on +x side of boundary, and halo cell not empty : bit 13
#define NF2_HALOPX	8192
//halo on -x side of boundary, and halo cell not empty : bit 14
#define NF2_HALONX	16384
//halo on +y side of boundary, and halo cell not empty : bit 15
#define NF2_HALOPY	32768
//halo on -y side of boundary, and halo cell not empty : bit 16
#define NF2_HALONY	65536
//halo on +z side of boundary, and halo cell not empty : bit 17
#define NF2_HALOPZ	131072
//halo on -z side of boundary, and halo cell not empty : bit 18
#define NF2_HALONZ	262144

//masks for x, y, z directions for halo cells
#define NF2_HALOX (NF2_HALOPX + NF2_HALONX)
#define NF2_HALOY (NF2_HALOPY + NF2_HALONY)
#define NF2_HALOZ (NF2_HALOPZ + NF2_HALONZ)
#define NF2_HALO (NF2_HALOX + NF2_HALOY + NF2_HALOZ)

private:

	//mark cells with various flags to indicate properties of neighboring cells
	int* ngbrFlags;

	//ngbrFlags2 defines additional flags. Only allocate memory if these additional flags are enabled - this is more memory efficient + I need to do it this way to keep older save files backward compatible.
	int* ngbrFlags2;

	//allocated size for ngbrFlags (value stored on gpu)
	size_t ngbrFlags_size;

	//indicates if ngbrFlags2 is in use, and is meant to speed up checks during computations only. This is strictly linked to the size of ngbrFlags2: if nullptr this is false, else it is true.
	bool using_extended_flags;

	//number of non-empty cells (i.e. cells  not marked with NF_NOTEMPTY)
	int nonempty_cells;

	//store dirichlet boundary conditions at mesh sides - only allocate memory as required
	//these vectors are of sizes equal to 1 cell deep at each respective side. dirichlet_nx are the dirichlet values at the -x side of the mesh, etc.
	VType* dirichlet_px;
	VType* dirichlet_nx;
	VType* dirichlet_py;
	VType* dirichlet_ny;
	VType* dirichlet_pz;
	VType* dirichlet_nz;

	//allocated sizes for dirichlet values (size values stored on gpu)
	size_t dirichlet_px_size;
	size_t dirichlet_nx_size;
	size_t dirichlet_py_size;
	size_t dirichlet_ny_size;
	size_t dirichlet_pz_size;
	size_t dirichlet_nz_size;

	//Robin boundary conditions values : diff_norm(u) = alpha * (u - cuVEC<VType>::h), where diff_norm means differential along surface normal (e.g. positive sign at +x boundary, negative sign at -x boundary).
	//alpha is a positive constant : robins_nx.i. Note, if this is zero then homogeneous Neumann boundary condition results.
	//cuVEC<VType>::h is a value (e.g. ambient temperature for heat equation): robins_nx.j
	//nx, px, etc... for mesh boundaries - use these when flagged with NF2_ROBINNX etc. and not flagged with NF2_ROBINV
	cuReal2 robin_px, robin_nx;
	cuReal2 robin_py, robin_ny;
	cuReal2 robin_pz, robin_nz;
	//robin_v applies for boundary conditions at void cells - more precisely for cells next to a void cell. Use this when flagged with NF2_ROBINNX etc. and also flagged with NF2_ROBINV
	cuReal2 robin_v;

	//when used with moving mesh algorithms calls to shift... functions may be used. If the shift requested is smaller than the cellsize then we cannot perform the shift. 
	//Add it to shift_debt and on next shift call we might be able to shift the mesh values.
	cuReal3 shift_debt;

	//Periodic boundary conditions for evaluating differential operators. If these are set then neighbor flags are calculated accordingly, and applied when evaluating operators.
	int pbc_x;
	int pbc_y;
	int pbc_z;

	//halo arrays applicable to boundaries of mesh at each respective side. either empty or at least 1 cell deep.
	//halos for each axis are allocated in pairs, i.e. both n and p sides are allocated
	VType* halo_px;
	VType* halo_nx;
	VType* halo_py;
	VType* halo_ny;
	VType* halo_pz;
	VType* halo_nz;

	//allocated sizes for halo arrays for each axis (size values stored on gpu) - n and p halos have same dimensions
	size_t halo_x_size;
	size_t halo_y_size;
	size_t halo_z_size;

private:

	//--------------------------------------------MEMORY MANAGEMENT HELPER METHODS : in cuVEC_VC_mng.h

	//memory allocation for objects and initialize to default - only call at start of managed constructors
	__host__ void alloc_initialize_data(void);

	//set size of ngbrFlags, allocating memory
	__host__ cudaError_t set_ngbrFlags_size(size_t size);

	//get ngbrFlags_size value in cpu memory
	__host__ size_t get_ngbrFlags_size(void) const;

	//get ngbrFlags2_size value in cpu memory
	__host__ size_t get_ngbrFlags2_size(void) const;

	//set robin values using a robin id
	__host__ void set_robin(cuReal2 robin_value, int robin_id);

	//set dirichlet array size, allocating memory as required
	__host__ cudaError_t set_dirichlet_size(size_t size, int dirichlet_id);

	//get dirichlet array size
	__host__ size_t get_dirichlet_size(int dirichlet_id) const;

	//set halo array size, allocating memory as required
	//halo_id is one of NF2_HALOX, NF2_HALOY, NF2_HALOZ, since both n and p halos are allocated for a given axis.
	//also accepts individual halo flags for halo_id
	__host__ cudaError_t set_halo_size(size_t size, int halo_id);

	//get halo array size
	//halo_id is one of NF2_HALOX, NF2_HALOY, NF2_HALOZ, since both n and p halos are allocated for a given axis.
	//also accepts individual halo flags for halo_id
	__host__ size_t get_halo_size(int halo_id) const;

	//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS : in cuVEC_VC_flags.cuh and cuVEC_VC_flags.h
	
	//set size of ngbrFlags to new_n also mapping shape from current size to new size (if current zero size set solid shape). Memory must be reserved in ngbrFlags to guarantee success. Also cuVEC<VType>::n should still have the old value : call this before changing it.
	__host__ cudaError_t resize_ngbrFlags(cuSZ3 new_n);

	//initialization method for neighbor flags : set flags at size cuVEC<VType>::n, counting neighbors etc.
	//Set empty cell values using information in linked_vec (keep same shape) - this must have same rectangle
	__host__ void set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);

	//initialization method for neighbor flags : set flags at size cuVEC<VType>::n, counting neighbors etc. Use current shape in ngbrFlags
	__host__ void set_ngbrFlags(void);
	
	//from NF2_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
	__device__ VType get_dirichlet_value(int dirichlet_flag, int idx) const;
	__device__ VType get_dirichlet_value(int dirichlet_flag, const cuINT3& ijk) const;
	
	//set robin flags from robin values and shape. Doesn't affect any other flags. Call from set_ngbrFlags after counting neighbors, and after setting robin values
	__host__ void set_robin_flags(void);

	//set pbc flags depending on set conditions and currently calculated flags - ngbrFlags must already be calculated before using this
	__host__ void set_pbc_flags(void);

	//mark cell as not empty / empty : internal use only; routines that use these must finish with recalculating ngbrflags as neighbours will have changed
	__device__ void mark_not_empty(int index) { ngbrFlags[index] |= NF_NOTEMPTY; }
	__device__ void mark_empty(int index) { ngbrFlags[index] &= ~NF_NOTEMPTY; cuVEC<VType>::quantity[index] = VType(); }

	//check if we need to use ngbrFlags2 (allocate memory etc.)
	__host__ bool use_extended_flags(void);

public:

	//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only. Real constructors are never called since you should never make a real instance of a cuVEC_VC. : cuVEC_VC_mng.h

	//void constructor
	__host__ void construct_cu_obj(void);

	//construct to given number of cells : n_ is in cpu memory
	__host__ void construct_cu_obj(const cuSZ3& n_);

	//construct to given dimensions : h_ and rect_ are in cpu memory
	__host__ void construct_cu_obj(const cuReal3& h_, const cuRect& rect_);

	//construct to given dimensions and initialize to given value : h_ and rect_ are in cpu memory
	__host__ void construct_cu_obj(const cuReal3& h_, const cuRect& rect_, VType value);

	//copy constructor
	__host__ void construct_cu_obj(const cuVEC_VC& copyThis);

	//assignment operator
	__host__ void assign_cu_obj(const cuVEC_VC& copyThis);

	//destructor
	__host__ void destruct_cu_obj(void);

	//--------------------------------------------SPECIAL ACCESS

	__host__ int*& ngbrFlags_ref(void) { return ngbrFlags; }

	__host__ VType*& halo_px_ref(void) { return halo_px; }
	__host__ VType*& halo_nx_ref(void) { return halo_nx; }
	__host__ VType*& halo_py_ref(void) { return halo_py; }
	__host__ VType*& halo_ny_ref(void) { return halo_ny; }
	__host__ VType*& halo_pz_ref(void) { return halo_pz; }
	__host__ VType*& halo_nz_ref(void) { return halo_nz; }

	//--------------------------------------------COPY TO / FROM VEC_VC

	//copy everything from a VEC_VC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC_VC>
	__host__ bool set_from_cpuvec(cpuVEC_VC& vec_vc);

	//copy everything to a VEC_VC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC_VC>
	__host__ bool set_cpuvec(cpuVEC_VC& vec_vc);

	//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC_VC
	template <typename cpuVEC_VC>
	__host__ bool copy_from_cpuvec(cpuVEC_VC& vec_vc);

	//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC_VC
	template <typename cpuVEC_VC>
	__host__ bool copy_to_cpuvec(cpuVEC_VC& vec_vc);

	//copy flags only from vec_vc, where sizes must match
	template <typename cpuVEC_VC>
	__host__ bool copyflags_from_cpuvec(cpuVEC_VC& vec_vc);

	//--------------------------------------------COPY TO ANOTHER cuVEC :  cuVEC_VC_mng.cuh

	//extract values from this and place them in cuvec : both must have same rectangle, but can differ in cuVEC<VType>::h: cuvec.h <= this->cuVEC<VType>::h needed (and hence cuVEC<VType>::n, where cuvec.cuVEC<VType>::n.dim() = size); e.g. this method allows extraction of a coarser cuvec.
	__host__ void extract_cuvec(size_t size, cuVEC_VC<VType>& cuvec);

	//--------------------------------------------EXTRACT REGION (useful for halo management):  cuVEC_VC_halo.cuh

	//extract ngbrFlags values within the given cells_box (end point not included, so e.g. cuBox(n) is the entire mesh), and place them inside the provided linear_storage space
	__host__ bool extract_ngbrFlags(cuBox cells_box, cu_arr<int>& linear_storage);

	//extract values from quantity into halo arrays as specified by halo_id
	//halo_id can be an individual flag, and identifies the required halos to extract values to, e.g. NF2_HALOPX, NF2_HALONX
	//values from p side of quantity are extracted to n halo and values from n side are extracted to p halo
	//e.g. if called with NF2_HALOPX, then values from n side are extracted to halo_px; if called with NF2_HALONX then values from p side are extracted to halo_nx
	__host__ void extract_halo(int halo_id, size_t size_transfer);

	//--------------------------------------------SIZING : cuVEC_VC_shape.h and cuVEC_VC_shape.cuh

	//sizing methods return true or false (failed to resize) - if failed then no changes made

	//resize and set shape using linked vec
	template <typename LVType>
	__host__ bool resize(cuSZ3 new_n, cuVEC_VC<LVType>& linked_vec);
	//resize but keep shape
	__host__ bool resize(cuSZ3 new_n);
	
	//resize and set shape using linked vec
	template <typename LVType>
	__host__ bool resize(cuReal3 new_h, cuRect new_rect, cuVEC_VC<LVType>& linked_vec);
	//resize but keep shape
	__host__ bool resize(cuReal3 new_h, cuRect new_rect);

	//set value and shape from linked vec
	template <typename LVType>
	__host__ bool assign(cuSZ3 new_n, VType value, cuVEC_VC<LVType>& linked_vec);
	//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
	__host__ bool assign(cuSZ3 new_n, VType value);

	//set value and shape from linked vec
	template <typename LVType>
	__host__ bool assign(cuReal3 new_h, cuRect new_rect, VType value, cuVEC_VC<LVType>& linked_vec);
	//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
	__host__ bool assign(cuReal3 new_h, cuRect new_rect, VType value);

	void clear(void);

	//--------------------------------------------FLAG CHECKING : cuVEC_VC_flags.h

	__device__ int get_nonempty_cells(void) const { return nonempty_cells; }
	__host__ int get_nonempty_cells_cpu(void);
	
	__device__ bool is_not_empty(int index) const { return (ngbrFlags[index] & NF_NOTEMPTY); }
	__device__ bool is_not_empty(const cuINT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*cuVEC<VType>::n.x + ijk.k*cuVEC<VType>::n.x*cuVEC<VType>::n.y] & NF_NOTEMPTY); }
	__device__ bool is_not_empty(const cuReal3& rel_pos) const { return (ngbrFlags[int(rel_pos.x / cuVEC<VType>::h.x) + int(rel_pos.y / cuVEC<VType>::h.y) * cuVEC<VType>::n.x + int(rel_pos.z / cuVEC<VType>::h.z) * cuVEC<VType>::n.x * cuVEC<VType>::n.y] & NF_NOTEMPTY); }

	__device__ bool is_empty(int index) const { return !(ngbrFlags[index] & NF_NOTEMPTY); }
	__device__ bool is_empty(const cuINT3& ijk) const { return !(ngbrFlags[ijk.i + ijk.j*cuVEC<VType>::n.x + ijk.k*cuVEC<VType>::n.x*cuVEC<VType>::n.y] & NF_NOTEMPTY); }
	__device__ bool is_empty(const cuReal3& rel_pos) const { return !(ngbrFlags[int(rel_pos.x / cuVEC<VType>::h.x) + int(rel_pos.y / cuVEC<VType>::h.y) * cuVEC<VType>::n.x + int(rel_pos.z / cuVEC<VType>::h.z) * cuVEC<VType>::n.x * cuVEC<VType>::n.y] & NF_NOTEMPTY); }

	__device__ bool is_empty(const cuRect& rectangle) const;

	__device__ bool is_not_cmbnd(int index) const { return !(ngbrFlags[index] & NF_CMBND); }
	__device__ bool is_not_cmbnd(const cuINT3& ijk) const { return !(ngbrFlags[ijk.i + ijk.j*cuVEC<VType>::n.x + ijk.k*cuVEC<VType>::n.x*cuVEC<VType>::n.y] & NF_CMBND); }
	__device__ bool is_not_cmbnd(const cuReal3& rel_pos) const { return !(ngbrFlags[int(rel_pos.x / cuVEC<VType>::h.x) + int(rel_pos.y / cuVEC<VType>::h.y) * cuVEC<VType>::n.x + int(rel_pos.z / cuVEC<VType>::h.z) * cuVEC<VType>::n.x * cuVEC<VType>::n.y] & NF_CMBND); }

	__device__ bool is_cmbnd(int index) const { return (ngbrFlags[index] & NF_CMBND); }
	__device__ bool is_cmbnd(const cuINT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*cuVEC<VType>::n.x + ijk.k*cuVEC<VType>::n.x*cuVEC<VType>::n.y] & NF_CMBND); }
	__device__ bool is_cmbnd(const cuReal3& rel_pos) const { return (ngbrFlags[int(rel_pos.x / cuVEC<VType>::h.x) + int(rel_pos.y / cuVEC<VType>::h.y) * cuVEC<VType>::n.x + int(rel_pos.z / cuVEC<VType>::h.z) * cuVEC<VType>::n.x * cuVEC<VType>::n.y] & NF_CMBND); }

	__device__ bool is_cmbnd_px(int index) const { return (ngbrFlags[index] & NF_CMBNDPX); }
	__device__ bool is_cmbnd_nx(int index) const { return (ngbrFlags[index] & NF_CMBNDNX); }
	__device__ bool is_cmbnd_py(int index) const { return (ngbrFlags[index] & NF_CMBNDPY); }
	__device__ bool is_cmbnd_ny(int index) const { return (ngbrFlags[index] & NF_CMBNDNY); }
	__device__ bool is_cmbnd_pz(int index) const { return (ngbrFlags[index] & NF_CMBNDPZ); }
	__device__ bool is_cmbnd_nz(int index) const { return (ngbrFlags[index] & NF_CMBNDNZ); }

	__device__ bool is_cmbnd_x(int index) const { return (ngbrFlags[index] & NF_CMBNDX); }
	__device__ bool is_cmbnd_y(int index) const { return (ngbrFlags[index] & NF_CMBNDY); }
	__device__ bool is_cmbnd_z(int index) const { return (ngbrFlags[index] & NF_CMBNDZ); }

	__device__ bool is_skipcell(int index) const { return (ngbrFlags[index] & NF_SKIPCELL); }
	__device__ bool is_skipcell(const cuINT3& ijk) const { return (ngbrFlags[ijk.i + ijk.j*cuVEC<VType>::n.x + ijk.k*cuVEC<VType>::n.x*cuVEC<VType>::n.y] & NF_SKIPCELL); }
	__device__ bool is_skipcell(const cuReal3& rel_pos) const { return (ngbrFlags[int(rel_pos.x / cuVEC<VType>::h.x) + int(rel_pos.y / cuVEC<VType>::h.y) * cuVEC<VType>::n.x + int(rel_pos.z / cuVEC<VType>::h.z) * cuVEC<VType>::n.x * cuVEC<VType>::n.y] & NF_SKIPCELL); }

	__device__ bool is_dirichlet(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLET); }

	__device__ bool is_dirichlet_px(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETPX); }
	__device__ bool is_dirichlet_nx(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETNX); }
	__device__ bool is_dirichlet_py(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETPY); }
	__device__ bool is_dirichlet_ny(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETNY); }
	__device__ bool is_dirichlet_pz(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETPZ); }
	__device__ bool is_dirichlet_nz(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETNZ); }

	__device__ bool is_dirichlet_x(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETX); }
	__device__ bool is_dirichlet_y(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETY); }
	__device__ bool is_dirichlet_z(int index) const { return using_extended_flags && (ngbrFlags2[index] & NF2_DIRICHLETZ); }

	//are all neighbors available? (for 2D don't check the z neighbors)
	__device__ bool is_interior(int index) const { return (((ngbrFlags[index] & NF_BOTHX) == NF_BOTHX) && ((ngbrFlags[index] & NF_BOTHY) == NF_BOTHY) && (cuVEC<VType>::n.z == 1 || ((ngbrFlags[index] & NF_BOTHZ) == NF_BOTHZ))); }
	
	//are all neighbors in the xy plane available?
	__device__ bool is_plane_interior(int index) const { return (((ngbrFlags[index] & NF_BOTHX) == NF_BOTHX) && ((ngbrFlags[index] & NF_BOTHY) == NF_BOTHY)); }

	//return number of neighbors present (pbcs not taken into consideration)
	__device__ int ngbr_count(int index) const { return ((ngbrFlags[index] & NF_NPX) == NF_NPX) + ((ngbrFlags[index] & NF_NNX) == NF_NNX) + ((ngbrFlags[index] & NF_NPY) == NF_NPY) + ((ngbrFlags[index] & NF_NNY) == NF_NNY) + ((ngbrFlags[index] & NF_NPZ) == NF_NPZ) + ((ngbrFlags[index] & NF_NNZ) == NF_NNZ); }

	//--------------------------------------------SET CELL FLAGS - EXTERNAL USE : cuVEC_VC_flags.h and cuVEC_VC_flags.cuh
	
	//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
	//return false on memory allocation failure only, otherwise return true even if surface_rect was not valid
	__host__ bool set_dirichlet_conditions(cuRect surface_rect, VType value);

	//clear all dirichlet flags and vectors
	__host__ void clear_dirichlet_flags(void);

	//set halo boundary conditions depending on halo_flag value (one of the 6 NF2_HALO.. flags); halo_depth is number of cells deep
	//also allocate memory and set size for respective hallo array (return false if memory allocation fails)
	//set halo flags using the halo_ngbrFlags values (this is an array with ngbrFlags values matching the halo region, so must be same size as halo array allocated here)
	//a halo flag is set if cell in halo is not empty (and current cell also not empty)
	//halo management is done externally (in policy class for mcuObject - mcuVEC - used to manage a cuVEC_VC across multiple GPUs).
	//This means halos are always cleared whenever set_ngbrFlags is triggered, so must be reset (e.g. after dimensions change or shape change, etc.) as managed by policy class
	//NOTE: this method does not set values in halo array, which is initialized to zero - halo values refresh must be done externally by computation routines as needed, through the policy class
	__host__ bool set_halo_conditions(cu_arr<int>& halo_ngbrFlags, int halo_flag, int halo_depth);

	//clear all halo flags and vectors
	__host__ void clear_halo_flags(void);

	//set pbc conditions : setting any to false clears flags
	__host__ void set_pbc(int pbc_x_, int pbc_y_, int pbc_z_);

	//clear all pbc flags : can also be achieved setting all flags to false in set_pbc but this one is more readable
	__host__ void clear_pbc(void);

	//clear all composite media boundary flags
	__host__ void clear_cmbnd_flags(void);

	//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations (if status true, else clear the skip cells flags in this rectangle)
	__host__ void set_skipcells(cuRect rectangle, bool status = true);

	//clear all skip cell flags
	__host__ void clear_skipcells(void);

	__host__ void set_robin_conditions(cuReal2 robin_v_, cuReal2 robin_px_, cuReal2 robin_nx_, cuReal2 robin_py_, cuReal2 robin_ny_, cuReal2 robin_pz_, cuReal2 robin_nz_);

	//clear all Robin boundary conditions and values
	__host__ void clear_robin_conditions(void);
	
	//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS : cuVEC_VC_shape.h and cuVEC_VC_shape.cuh

	//set value in box (i.e. in cells entirely included in box) - all cells become non-empty cells irrespective of value set
	__host__ void setbox(cuBox box, VType value = VType());

	//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle - all cells become non-empty cells irrespective of value set
	__host__ void setrect(cuRect rectangle, VType value = VType());

	//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells
	__host__ void delrect(cuRect rectangle);

	//mask values in cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
	__host__ bool apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth = 0);
	
	//--------------------------------------------MULTIPLE ENTRIES SETTERS - OTHERS : cuVEC_oper.h and cuVEC_oper.cuh

	//exactly the same as assign value - do not use assign as it is slow (sets flags)
	__host__ void setnonempty(VType value = VType());

	//set value in non-empty cells only in given rectangle (relative coordinates)
	__host__ void setrectnonempty(const cuRect& rectangle, VType value = VType());

	//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
	//Launch it with arr_size = cuVEC<VType>::n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(cuVEC<VType>::n).dim()
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ void renormalize(size_t arr_size, PType new_norm);

	//shift all the values in this cuVEC by the given delta (units same as cuVEC<VType>::h). Shift values in given shift_rect (absolute coordinates).
	__host__ void shift_x(size_t size, cuBReal delta, cuRect shift_rect);

	//--------------------------------------------ARITHMETIC OPERATIONS ON ENTIRE VEC : cuVEC_VC_arith.cuh

	//scale all stored values by the given constant
	__host__ void scale_values(size_t size, cuBReal constant);
	__host__ void scale_values(cuBReal constant) { scale_values(get_gpu_value(cuVEC<VType>::n).dim(), constant); }

	//--------------------------------------------AVERAGING / SUMMING OPERATIONS : cuVEC_VC_avg.cuh

	//overload VEC method : use NF_NOTEMPTY flags instead here
	//Launch it with arr_size = cuVEC<VType>::n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(cuVEC<VType>::n).dim()
	__host__ VType average_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu = true);
	//average over non-empty cells over given rectangle (relative to this VEC's rect)
	__host__ VType average_nonempty(size_t arr_size, cuRect rectangle = cuRect(), bool transfer_to_cpu = true);

	//Launch it with arr_size = cuVEC<VType>::n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(cuVEC<VType>::n).dim()
	__host__ VType sum_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu = true);
	//sum over non-empty cells over given rectangle (relative to this VEC's rect)
	__host__ VType sum_nonempty(size_t arr_size, cuRect rectangle = cuRect(), bool transfer_to_cpu = true);

	//--------------------------------------------NUMERICAL PROPERTIES : cuVEC_VC_nprops.cuh
	
	//Find min and max values from non-empty points. rectangles are relative to this VEC.
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax(size_t arr_size, cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax(size_t arr_size, cuRect rectangle = cuRect());

	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax_component_x(size_t arr_size, cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax_component_x(size_t arr_size, cuRect rectangle = cuRect());

	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax_component_y(size_t arr_size, cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax_component_y(size_t arr_size, cuRect rectangle = cuRect());

	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax_component_z(size_t arr_size, cuBox box);
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	__host__ cuVAL2<PType> get_minmax_component_z(size_t arr_size, cuRect rectangle = cuRect());

	//--------------------------------------------SPECIAL NUMERICAL PROPERTIES : cuVEC_VC_nprops.cuh

	//Robin value is of the form alpha * (Tb - Ta). alpha and Ta are known from values set robin_nx, robin_px, ...
	//Tb is quantity value at boundary. Here we will return Robin value for x, y, z axes for which any shift component is nonzero (otherwise zero for that component)
	//e.g. if shift.x is non-zero then Tb value is obtained at rel_pos + (shift.x, 0, 0) using extrapolation from values at rel_pos and rel_pos - (shift.x, 0, 0) -> both these values should be inside the mesh, else return zero.
	__device__ cuReal3 get_robin_value(const cuReal3& rel_pos, const cuReal3& shift);

	//for given cell index, find if any neighboring cells are empty and get distance (shift) valuer to them along each axis
	//if any shift is zero this means both cells are present either side, or both are missing
	//NOTE : this is intended to be used with get_robin_value method to determine the shift value, and rel_pos will be position corresponding to idx
	__device__ cuReal3 get_shift_to_emptycell(int idx);

	//--------------------------------------------CALCULATE COMPOSITE MEDIA BOUNDARY VALUES : cuVEC_VC_cmbnd.h

	//calculate and set values at composite media boundary cells in this mesh for a given contacting mesh (the "secondary" V_sec) and given contact description (previously calculated using set_cmbnd_flags)
	//The values are calculated based on the continuity of a potential and flux normal to the interface. The potential is this VEC_VC, call it V, and the flux is the function f(V) = a_func + b_func * V', where the V' differential direction is perpendicular to the interface.
	//The second order differential of V perpendicular to the interface, V'', is also used and specified using the methods diff2 in Class_CMBND
	//Boundary with labelled cells either side, on the left is the secondary mesh, on the right is the primary mesh (i.e. the mesh for which we are setting boundary values using this method now): -2 -1 | 1 2 
	//Functions which must be defined in Class_CMBND:
	//a_func_sec is for the secondary side and takes a position for cell -1, a position shift to add to position to reach cell -2 and finally a stencil to use when obtaining values at cells -1 and -2 (use weighted_average)
	//b_func_sec similar
	//diff2_sec takes a position and stencil (since only need cell -1). It also takes a position shift vector perpendicular to the interface and pointing from primary to secondary.
	//a_func_pri takes indexes for cells 1 and 2. It also takes a position shift vector perpendicular to the interface and pointing from primary to secondary.
	//b_func_pri takes indexes for cells 1 and 2
	//diff2_pri takes index for cell 1. It also takes a position shift vector perpendicular to the interface and pointing from primary to secondary.
	//also need instances for the secondary and primary objects whose classes contain the above methods
	template <typename cuVEC_sec, typename Class_CMBNDs, typename Class_CMBNDp>
	__host__ void set_cmbnd_continuous(size_t size, cu_obj<cuVEC_sec>& V_sec, Class_CMBNDs& cmbndFuncs_sec, Class_CMBNDp& cmbndFuncs_pri, CMBNDInfoCUDA& contact);

	//calculate cmbnd values based on continuity of flux only. The potential is allowed to drop across the interface as:
	//f_sec(V) = f_pri(V) = A + B * delV, where f_sec and f_pri are the fluxes on the secondary and primary sides of the interface, and delV = V_pri - V_sec, the drop in potential across the interface.
	//Thus in addition to the functions in set_cmbnd_continuous we need two extra functions A_func, B_func, defined in Class_CMBND_S
	//Note, this reduces to the fully continuous case for B tending to infinity and A = VType(0)
	template <typename Class_CMBNDs, typename Class_CMBNDp, typename Class_CMBND_S>
	__host__ void set_cmbnd_continuousflux(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBNDs& cmbndFuncs_sec, Class_CMBNDp& cmbndFuncs_pri, Class_CMBND_S& cmbndFuncs_s, CMBNDInfoCUDA& contact);

	//most general case of composite media boundary conditions
	//calculate cmbnd values based on set boundary flux values; both the flux and potential is allowed to be discontinuous across the interface.
	//Fluxes at the interface are specified as: f_sec(V) = A_sec + B_sec * delVs, f_pri(V) = A_pri + B_pri * delVs, with directions from secondary to primary
	//B functions may return a double, or a DBL33 (3x3 matrix) if VType is DBL3 (Cartesian vector).
	//delVs = c_pri * V_pri - c_sec * V_sec, c are double values specified by given functions
	template <typename Class_CMBNDs, typename Class_CMBNDp, typename Class_CMBND_S>
	__host__ void set_cmbnd_discontinuous(size_t size, cuVEC_VC<VType>& V_sec, Class_CMBNDs& cmbndFuncs_sec, Class_CMBNDp& cmbndFuncs_pri, Class_CMBND_S& cmbndFuncs_s, CMBNDInfoCUDA& contact);

	//--------------------------------------------OPERATORS and ALGORITHMS

	//----LAPLACE OPERATOR : cuVEC_VC_del.h

	//calculate Laplace operator at cell with given index. Use Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	__device__ VType delsq_neu(int idx) const;

	//calculate Laplace operator at cell with given index. Use non-homogeneous Neumann boundary conditions with the specified boundary differential.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Returns zero at composite media boundary cells.
	template <typename Class_BDiff>
	__device__ VType delsq_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType delsq_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	__device__ VType delsq_diri(int idx) const;

	//calculate Laplace operator at cell with given index. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Returns zero at composite media boundary cells.
	template <typename Class_BDiff>
	__device__ VType delsq_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType delsq_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//calculate Laplace operator at cell with given index. Use Robin boundary conditions (defaulting to Neumann if not set).
	//Returns zero at composite media boundary cells.
	//The K constant is used in Robin boundary condition calculations, where -K*diff_norm(T) = alpha*(Tboundary - Tambient) is the flux normal to the boundary - K is the thermal conductivity in the heat equation
	__device__ VType delsq_robin(int idx, cuBReal K) const;

	//----GRADIENT OPERATOR : cuVEC_VC_grad.h

	//gradient operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ cuVAL3<VType> grad_neu(int idx) const;

	//gradient operator. Use non-homogeneous Neumann boundary conditions.
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	template <typename Class_BDiff>
	__device__ cuVAL3<VType> grad_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuVAL3<VType> grad_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//gradient operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ cuVAL3<VType> grad_diri(int idx) const;

	//gradient operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	template <typename Class_BDiff>
	__device__ cuVAL3<VType> grad_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuVAL3<VType> grad_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//gradient operator. Use sided differentials at boundaries (including at composite media boundaries)
	__device__ cuVAL3<VType> grad_sided(int idx) const;

	//----DIVERGENCE OPERATOR : cuVEC_VC_div.h

	//divergence operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	__device__ cuBReal div_neu(int idx) const;

	//divergence operator. Use non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	template <typename Class_BDiff>
	__device__ cuBReal div_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuBReal div_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//divergence operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	__device__ cuBReal div_diri(int idx) const;

	//divergence operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	template <typename Class_BDiff>
	__device__ cuBReal div_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ cuBReal div_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//divergence operator. Use sided differentials (also at composite media boundaries)
	//div operator can be applied if VType is a VAL3<Type>, returning Type
	__device__ cuBReal div_sided(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ VType diveps3_neu(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	__device__ VType diveps3_diri(int idx) const;

	//divergence operator of epsilon3(quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
	//Use sided differentials (also at composite media boundaries)
	__device__ VType diveps3_sided(int idx) const;

	//----CURL OPERATOR : cuVEC_VC_curl.h

	//curl operator. Use Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	__device__ VType curl_neu(int idx) const;

	//curl operator. Use non-homogeneous Neumann boundary conditions.
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//can only be applied if VType is a VAL3
	template <typename Class_BDiff>
	__device__ VType curl_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType curl_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//curl operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	__device__ VType curl_diri(int idx) const;

	//curl operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//can only be applied if VType is a VAL3
	template <typename Class_BDiff>
	__device__ VType curl_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;

	//Same as above but boundary conditions specified using a constant
	__device__ VType curl_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//curl operator. Use sided differentials at boundaries (including at composite media boundaries)
	//can only be applied if VType is a VAL3
	__device__ VType curl_sided(int idx) const;

	//---- SECOND ORDER DIFFERENTIALS : cuVEC_VC_diff2.h

	//homogeneous second order.
	//Use Neumann boundary conditions.
	//Returns zero at composite media boundary cells
	__device__ VType dxx_neu(int idx) const;
	__device__ VType dyy_neu(int idx) const;
	__device__ VType dzz_neu(int idx) const;

	//Use non-homogeneous Neumann boundary conditions.
	//Returns zero at composite media boundary cells
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	template <typename Class_BDiff>
	__device__ VType dxx_nneu(int idx, const Class_BDiff& bdiff_class) const;

	template <typename Class_BDiff>
	__device__ VType dyy_nneu(int idx, const Class_BDiff& bdiff_class) const;

	template <typename Class_BDiff>
	__device__ VType dzz_nneu(int idx, const Class_BDiff& bdiff_class) const;

	__device__ VType dxx_nneu(int idx, const cuVAL3<VType>& bdiff) const;
	__device__ VType dyy_nneu(int idx, const cuVAL3<VType>& bdiff) const;
	__device__ VType dzz_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//Use Dirichlet boundary conditions, else Neumann boundary conditions (homogeneous).
	//Returns zero at composite media boundary cells.
	__device__ VType dxx_diri(int idx) const;
	__device__ VType dyy_diri(int idx) const;
	__device__ VType dzz_diri(int idx) const;

	//Use Dirichlet boundary conditions, else non-homogeneous Neumann boundary conditions.
	//Returns zero at composite media boundary cells.
	//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
	template <typename Class_BDiff>
	__device__ VType dxx_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;
	
	template <typename Class_BDiff>
	__device__ VType dyy_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;
	
	template <typename Class_BDiff>
	__device__ VType dzz_diri_nneu(int idx, const Class_BDiff& bdiff_class) const;

	__device__ VType dxx_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;
	__device__ VType dyy_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;
	__device__ VType dzz_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const;

	//Use Robin boundary conditions (defaulting to Neumann if not set).
	//Returns zero at composite media boundary cells.
	//The K constant is used in Robin boundary condition calculations, where -K*diff_norm(T) = alpha*(Tboundary - Tambient) is the flux normal to the boundary - K is the thermal conductivity in the heat equation
	__device__ VType dxx_robin(int idx, cuBReal K) const;
	__device__ VType dyy_robin(int idx, cuBReal K) const;
	__device__ VType dzz_robin(int idx, cuBReal K) const;

	//mixed second order

	//Use Neumann boundary conditions(homogeneous).
	//Can be used at composite media boundaries where sided differentials will be used instead.
	//dxy same as dyx
	__device__ VType dxy_neu(int idx) const;

	//dxz same as dzx
	__device__ VType dxz_neu(int idx) const;

	//dyz same as dzy
	__device__ VType dyz_neu(int idx) const;

	//----NEIGHBOR SUM : cuVEC_VC_ngbrsum.h

	//calculate 6-point neighbor sum at given index
	//missing neighbors not added, including at boundaries, but taking into account pbc
	__device__ VType ngbr_sum(int idx) const;

	//same as ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
	__device__ VType ngbr_dirsum(int idx) const;

	//calculate 6-point anisotropic neighbor sum at given index as rij x Vj over j points neighboring the point i at this index.
	//missing neighbors not added, including at boundaries, but taking into account pbc
	//only used if VType is a cuVAL3
	__device__ VType anisotropic_ngbr_sum(int idx) const;

	//same as anisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
	__device__ VType anisotropic_ngbr_dirsum(int idx) const;

	//calculate 6-point anisotropic neighbor sum at given index as (rij x z) x Vj over j points neighboring the point i at this index.
	//missing neighbors not added, including at boundaries, but taking into account pbc
	//only used if VType is a cuVAL3
	__device__ VType zanisotropic_ngbr_sum(int idx) const;
	//calculate 6-point anisotropic neighbor sum at given index as (rij x y) x Vj over j points neighboring the point i at this index.
	__device__ VType yanisotropic_ngbr_sum(int idx) const;
	//calculate 6-point anisotropic neighbor sum at given index as (rij x x) x Vj over j points neighboring the point i at this index.
	__device__ VType xanisotropic_ngbr_sum(int idx) const;

	//same as zanisotropic_ngbr_sum but sum normalised values only; for scalar values this is a trivial operation, but for vectors it's not.
	__device__ VType zanisotropic_ngbr_dirsum(int idx) const;
	__device__ VType yanisotropic_ngbr_dirsum(int idx) const;
	__device__ VType xanisotropic_ngbr_dirsum(int idx) const;

	//--------------------------------------------LAPLACE / POISSON EQUATION : cuVEC_VC_solve.cuh
	
	//LAPLACE

	//Take one SOR iteration for Laplace equation on this VEC. Return error (maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a cu_obj managed cuBReal.
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	//Launch it with arr_size = cuVEC<VType>::n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(cuVEC<VType>::n).dim()
	__host__ void IterateLaplace_SOR(size_t arr_size, cuBReal& damping, cuBReal& max_error, cuBReal& max_val);

	//as above, but allow separate control of red and black passes launch
	__host__ void IterateLaplace_SOR_red(size_t arr_size, cuBReal& damping);
	__host__ void IterateLaplace_SOR_black(size_t arr_size, cuBReal& damping, cuBReal& max_error, cuBReal& max_val);

	//red and black SOR passes launched from kernels : the flow is : 
	//IterateLaplace_SOR -> 1. global red pass kernel with *this as parameter -> call IterateLaplace_SOR_red_ondevice, 
	//2. global black pass kernel with *this as parameter -> call IterateLaplace_SOR_black_ondevice
	__device__ void IterateLaplace_SOR_red_ondevice(cuBReal damping);
	__device__ void IterateLaplace_SOR_black_ondevice(cuBReal damping, cuBReal& max_error, cuBReal& max_val);

	//POISSON with homogeneous Neumann boundaries

	//NOTE : Poisson variants can only be called in cu files (where it is possible to include cuVEC_VC_solve.cuh), otherwise explicit template parameters would have to be declared, which is too restrictive.

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = Poisson_RHS
	//Poisson_RHS must be a member const method of Class_Poisson_RHS taking an index value (the index ranges over this VEC) and returning a cuBReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//obj of type Class_Poisson_RHS must be cu_obj managed so it is entirely in gpu memory.
	//Return error(maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a cu_obj managed cuBReal.
	//Dirichlet boundary conditions used, defaulting to Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	template <typename Class_Poisson_RHS>
	__host__ void IteratePoisson_SOR(size_t arr_size, Class_Poisson_RHS& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val);

	//as above, but allow separate control of red and black passes launch
	template <typename Class_Poisson_RHS>
	__host__ void IteratePoisson_SOR_red(size_t arr_size, Class_Poisson_RHS& obj, cuBReal& damping);
	template <typename Class_Poisson_RHS>
	__host__ void IteratePoisson_SOR_black(size_t arr_size, Class_Poisson_RHS& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val);

	//red and black SOR passes launched from kernels : the flow is : 
	//IteratePoisson_SOR -> 1. global red pass kernel with *this as parameter -> call IteratePoisson_SOR_red_ondevice, 
	//2. global black pass kernel with *this as parameter -> call IteratePoisson_SOR_black_ondevice
	template <typename Class_Poisson_RHS>
	__device__ void IteratePoisson_SOR_red_ondevice(Class_Poisson_RHS& obj, cuBReal damping);
	template <typename Class_Poisson_RHS>
	__device__ void IteratePoisson_SOR_black_ondevice(Class_Poisson_RHS& obj, cuBReal damping, cuBReal& max_error, cuBReal& max_val);

	//POISSON with non-homogeneous Neumann boundaries

	//For Poisson equation we need a function to specify the RHS of the equation delsq V = Poisson_RHS
	//Poisson_RHS must be a member const method of Class_Poisson_NNeu taking an index value (the index ranges over this VEC) and returning a cuBReal value : Poisson_RHS(index) evaluated at the index-th cell.
	//Class_Poisson_NNeu must also define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index) - this is the non-homogeneous Neumann boundary condition at that cell
	//obj of type Class_Poisson_NNeu must be cu_obj managed so it is entirely in gpu memory.
	//Return error(maximum change in quantity from one iteration to the next) by reference, where max_error is already allocated on the gpu - pass in a cu_obj managed cuBReal.
	//Dirichlet boundary conditions used, defaulting to non-homogeneous Neumann boundary conditions where not set, and composite media boundary cells skipped (use boundary conditions to set cmbnd cells after calling this)
	template <typename Class_Poisson_NNeu>
	__host__ void IteratePoisson_NNeu_SOR(size_t arr_size, Class_Poisson_NNeu& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val);

	//as above, but allow separate control of red and black passes launch
	template <typename Class_Poisson_NNeu>
	__host__ void IteratePoisson_NNeu_SOR_red(size_t arr_size, Class_Poisson_NNeu& obj, cuBReal& damping);
	template <typename Class_Poisson_NNeu>
	__host__ void IteratePoisson_NNeu_SOR_black(size_t arr_size, Class_Poisson_NNeu& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val);

	//red and black SOR passes launched from kernels : the flow is : 
	//IteratePoisson_NNeu_SOR -> 1. global red pass kernel with *this as parameter -> call IteratePoisson_NNeu_SOR_red_ondevice, 
	//2. global black pass kernel with *this as parameter -> call IteratePoisson_NNeu_SOR_black_ondevice
	template <typename Class_Poisson_NNeu>
	__device__ void IteratePoisson_NNeu_SOR_red_ondevice(Class_Poisson_NNeu& obj, cuBReal damping);
	template <typename Class_Poisson_NNeu>
	__device__ void IteratePoisson_NNeu_SOR_black_ondevice(Class_Poisson_NNeu& obj, cuBReal damping, cuBReal& max_error, cuBReal& max_val);
};

