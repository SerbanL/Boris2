#pragma once

#include "cuTypes.h"
#include "cuFuncs_Aux.h"
#include "alloc_cpy.h"
#include "cuArray.h"
#include "cuObject.h"

#include "Types_VAL.h"

#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////// cuVEC<VType>
//
// n-component quantity with 3 dimensions for CUDA gpu computations.
// Important : This is meant to handle only gpu-addressable memory and you should manage it using a cu_obj. 
// Cannot have any data stored in cpu memory held here since cu_obj managed objects exist only in gpu memory.
// Also cu_obj managed objects are never actually created or destroyed with convential ctors and dtor : use construct_cu_obj and destruct_cu_obj instead.

template <typename VType> class cuTransfer;

template <typename VType> 
class cuVEC
{
	friend cuTransfer<VType>;

private:

	//Line profile extraction data

	//all profile components internal storage
	VType* line_profile;

	//extract a single component, together with profile position, for cuVAL3 types : auxiliary storage
	cuReal2* line_profile_component_x;
	cuReal2* line_profile_component_y;
	cuReal2* line_profile_component_z;

	//number of points included in line profile average (need to keep track in case there are empty points) for each point
	size_t* line_profile_avpoints;

	//line_profile_component memory size, stored in GPU memory : transfer to CPU before checking if correct size
	size_t line_profile_component_size;

	//Histogram extraction data

	//used for extracting histogram data from mesh
	cuBReal* histogram;
	size_t histogram_size;

	//used if we need to pre-average mesh values with a macrocell before extracting histogram
	VType* histogram_preaverage;
	size_t histogram_preaverage_size;

protected:

	//the actual mesh quantity: addresses gpu memory
	VType* quantity;

	//separate object to manage transfer info
	cuTransfer<VType> transfer;

	//pre-allocated memory used for calculations
	VType aux_value, aux_value2, aux_value3;
	cuBReal aux_real, aux_real2;
	size_t aux_integer;

	//pre-allocate array used for reductions : this array has size (n.dim() + CUDATHREADS) / CUDATHREADS), i.e. one entry per block used
	VType* aux_block_values;

public:

	//NOTE : n, h, rect must not be modifed in device code. Only modify them in host code using the set_ methods.

	//dimensions along x, y and z of the quantity : held in gpu memory
	cuSZ3 n;

	//cellsize of structured mesh : held in gpu memory
	cuReal3 h;

	//rectangle, same units as h. cuVEC has n number of cells, so n * h gives the rect dimensions. All that is really needed is the rect start coordinates : held in gpu memory
	cuRect rect;

private:
	
	//--------------------------------------------MEMORY MANAGEMENT HELPER METHODS : in cuVEC_mng.h

	//memory allocation for objects and initialize to default - only call at start of managed constructors
	__host__ void alloc_initialize_data(void);

	//--------------------------------------------HELPER METHODS : cuVEC_mng.h (and cuVEC_mng.cuh)

	__host__ void SetMeshRect(void);

	//set new size and map mesh values to new dimension, keeping magnitude (so don't use an average). Return outcome; if failed then no changes made. This launches a kernel mapmesh_newdims_kernel in cuVEC.cuh.
	__host__ bool mapmesh_newdims(const cuSZ3& new_n);

	//from current rectangle and h value set n. h may also need to be adjusted since n must be an integer. Resize quantity to new n value : return success or fail. If failed then nothing changed.
	__host__ bool set_n_adjust_h(void);

	//memory management for line_profile_component array : attempt to resize to new size if not already exactly given size
	__host__ bool allocate_profile_component_memory(size_t size);

	//memory management for histogram array : attempt to resize to new size if not already exactly given size
	__host__ bool allocate_histogram_memory(size_t histogram_size_cpu, size_t histogram_preaverage_size_cpu);

protected:

	//--------------------------------------------GET/SET FROM/TO GPU MEMORY : cuVEC_mng.h

	__host__ cudaError allocate_quantity(cuSZ3 new_n);

	//set n in gpu memory from cpu memory value - can take both l-values and r-values
	__host__ void set_n(const cuSZ3& n_);

	//set h in gpu memory from cpu memory value - can take both l-values and r-values
	__host__ void set_h(const cuReal3& h_);

	//set rect gpu memory from cpu memory value - can take both l-values and r-values
	__host__ void set_rect(const cuRect& rect_);

	//from h_ and rect_ (in cpu memory) calculate what n value results (in cpu memory) - but do not make any changes
	__host__ cuSZ3 get_n_from_h_and_rect(const cuReal3& h_, const cuRect& rect_);

public:

	//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only. Real constructors are never called since you should never make a real instance of a cuVEC. : cuVEC_mng.h

	//void constructor
	__host__ void construct_cu_obj(void);

	//construct to given number of cells : n_ is in cpu memory
	__host__ void construct_cu_obj(const cuSZ3& n_);
	
	//construct to given dimensions : h_ and rect_ are in cpu memory
	__host__ void construct_cu_obj(const cuReal3& h_, const cuRect& rect_);

	//construct to given dimensions and initialize to given value : h_ and rect_ are in cpu memory
	__host__ void construct_cu_obj(const cuReal3& h_, const cuRect& rect_, VType value);

	//copy constructor
	__host__ void construct_cu_obj(const cuVEC& copyThis);

	//assignment operator
	__host__ void assign_cu_obj(const cuVEC& copyThis);

	__host__ void destruct_cu_obj(void);

	//--------------------------------------------COPY TO / FROM VEC : cuVEC_mng.h

	//copy everything from a VEC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC>
	__host__ bool set_from_cpuvec(cpuVEC& vec);

	//copy everything to a VEC - type must be convertible. Return false if failed (memory could not be allocated)
	template <typename cpuVEC>
	__host__ bool set_cpuvec(cpuVEC& vec);

	//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC : only quantity is copied.
	template <typename cpuVEC>
	__host__ bool copy_from_cpuvec(cpuVEC& vec);

	//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC : only quantity is copied.
	template <typename cpuVEC>
	__host__ bool copy_to_cpuvec(cpuVEC& vec);

	//copy values from a cu_arr of same type -> sizes must match
	__host__ void load_cuarr(size_t size, cu_arr<VType>& input);

	//copy values to a cu_arr of same type -> sizes must match
	__host__ void store_cuarr(size_t size, cu_arr<VType>& output);

	template <typename SType>
	__host__ bool copy_from_vector(std::vector<SType>& vec);

	//as above but specifically for std::vector
	template <typename SType>
	__host__ bool copy_to_vector(std::vector<SType>& vec);

	//--------------------------------------------COPY TO ANOTHER cuVEC :  cuVEC_mng.cuh

	//extract values from this and place them in cuvec : both must have same rectangle, but can differ in h - cuvec.h <= this->h needed (and hence n, where cuvec.n.dim() = size); e.g. this method allows extraction of a coarser cuvec.
	__host__ void extract_cuvec(size_t size, cuVEC<VType>& cuvec);

	//--------------------------------------------EXTRACT A LINE PROFILE : cuVEC_extract.cuh

	//auxiliary access
	__device__ VType*& get_line_profile(void) { return line_profile; }
	__device__ cuReal2*& get_line_profile_component_x(void) { return line_profile_component_x; }
	__device__ cuReal2*& get_line_profile_component_y(void) { return line_profile_component_y; }
	__device__ cuReal2*& get_line_profile_component_z(void) { return line_profile_component_z; }

	__device__ size_t get_line_profile_size(void) { return line_profile_component_size; }

	//1. Profile values only, without stencil operation, and with mesh wrap-around

	//extract profile to a cu_arr : extract size points starting at start in the direction step for the given number of points (size); use weighted average to extract profile with h stencil only
	__host__ void extract_profilevalues(size_t size, cu_arr<VType>& profile_gpu, cuReal3 start, cuReal3 step);

	//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
	__host__ void extract_profilevalues_component_x(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
	__host__ void extract_profilevalues_component_y(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);
	__host__ void extract_profilevalues_component_z(size_t size, cu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step);

	//2. Profile values only, with stencil operation around profile point (capped to mesh size), and without mesh wrap-around

	//extract profile components: extract starting at start in the direction end - step, with given step; use weighted average to extract profile with given stencil
	//all coordinates are relative positions. Return profile values in profile_gpu.
	__host__ bool extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<VType>& profile_gpu);

	//as above, but only store profile in internal memory (line_profile) so we can read it out later as needed
	__host__ bool extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil);

	//extract profile components: extract starting at start in the direction end - step, with given step; use weighted average to extract profile with given stencil
	//all coordinates are relative positions. Return profile values in profile_cpu.
	template <typename SType>
	__host__ bool extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<SType>& profile_cpu);

	//3. Profile values for individual components together with profile position returned in cuReal2/DBL2, with stencil operation around profile point (capped to mesh size), and without mesh wrap-around

	//as above but only component x and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu. profile_gpu must have correct size: (end - start).norm() / step + 1.
	__host__ bool extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
	//as above but only component x and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
	__host__ bool extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);
	
	//as above but only component y and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu. profile_gpu must have correct size: (end - start).norm() / step + 1.
	__host__ bool extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
	//as above but only component y and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
	__host__ bool extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

	//as above but only component z and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu. profile_gpu must have correct size: (end - start).norm() / step + 1.
	__host__ bool extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
	//as above but only component z and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
	__host__ bool extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);
	
	//as above but only component which has largest value for the first point and pack in profile position too (after stencil averaging) (for VAL3 floating types only). Return data in profile_gpu. profile_gpu must have correct size: (end - start) / step.norm() + 1.
	__host__ bool extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, cu_arr<cuReal2>& profile_gpu);
	//as above but only component which has largest value for the first point and pack in profile position too (after stencil averaging) (for VAL3 floating types only). Return data in profile_cpu. profile_cpu resized as needed.
	__host__ bool extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu);

	//--------------------------------------------HISTOGRAMS : cuVEC_histo.cuh

	//compute magnitude histogram data
	//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first.
	//if num_bins not given then use default value of 100
	//if macrocell_dims greater than 1 in any dimension then first average mesh data in macrocells of given size
	//without macrocell then pass in num_nonempty_cells (number of nonempty cells); if this is not passed it it's counted first (costs another kernel launch)
	//output transferred to cpu in histogram_cpu
	__host__ bool get_mag_histogram(
		std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, 
		int num_bins, double& min, double& max, size_t num_nonempty_cells = 0, cuINT3 macrocell_dims = cuINT3(1));

	//get angular deviation histogram. deviation from ndir direction is calculated, or the average direction if ndir not specified (IsNull)
	__host__ bool get_ang_histogram(
		std::vector<double>& histogram_x_cpu, std::vector<double>& histogram_p_cpu, 
		int num_bins, double& min, double& max, size_t num_nonempty_cells = 0, cuINT3 macrocell_dims = cuINT3(1), VType ndir = VType());

	//--------------------------------------------INDEXING

	//Index using a single combined index (use e.g. when more convenient to use a single for loop to iterate over the quantity's elements)
	__device__ VType& operator[](int idx) { return quantity[idx]; }

	//index using a cuVAL3, integral type (e.g. use with nested loops)
	__device__ VType& operator[](const cuINT3& idx) { return quantity[idx.i + idx.j*n.x + idx.k*n.x*n.y]; }

	//index by position relative to cuVEC rect
	__device__ VType& operator[](const cuReal3& rel_pos) { return quantity[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y]; }

	//--------------------------------------------PROPERTY CHECKING

	__device__ bool is_not_empty(int index) { return (quantity[index] != VType()) ; }
	__device__ bool is_not_empty(const cuINT3& ijk) { return (quantity[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] != VType()); }
	__device__ bool is_not_empty(const cuReal3& rel_pos) { return (quantity[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y] != VType()); }

	__device__ bool is_empty(int index) { return (quantity[index] == VType()); }
	__device__ bool is_empty(const cuINT3& ijk) { return (quantity[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] == VType()); }
	__device__ bool is_empty(const cuReal3& rel_pos) { return (quantity[int(rel_pos.x / h.x) + int(rel_pos.y / h.y) * n.x + int(rel_pos.z / h.z) * n.x * n.y] == VType()); }

	//--------------------------------------------ITERATORS
	
	__host__ __device__ VType* begin(void) { return &quantity[0]; }
	__host__ __device__ VType* end(void) { return &quantity[linear_size()]; }
	__host__ __device__ VType* data(void) { return quantity; }
	
	__host__ VType*& quantity_ref(void) { return quantity; }

	//--------------------------------------------SIZING : cuVEC_mng.h

	//all sizing methods (apart from clear) return true (success) or false (could not resize). If failed then previous settings are maintained.

	//change to new number of cells : keep h and rect.s the same but adjust rect.e. Also map values to new size.
	__host__ bool resize(cuSZ3 new_n);
	//set rect and h; n is obtained from them and h also may be adjusted. Also map values to new size.
	__host__ bool resize(cuReal3 new_h, cuRect new_rect);

	//works like resize but sets given value also
	__host__ bool assign(cuSZ3 new_n, VType value);
	//works like resize but sets given value also
	__host__ bool assign(cuReal3 new_h, cuRect new_rect, VType value);

	//set everything to zero but h
	__host__ void clear(void);

	//--------------------------------------------MULTIPLE ENTRIES SETTERS : cuVEC_oper.h (and cuVEC_oper.cuh)

	//set value in box
	__host__ void setbox(cuBox box, VType value = VType());

	//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this cuVEC's rectangle.
	__host__ void setrect(const cuRect& rectangle, VType value = VType());

	//set value in all cells
	__host__ void set(size_t size, VType value);
	__host__ void set(VType value) { set(get_gpu_value(n).dim(), value); }

	//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
	template <typename PType = decltype(cu_GetMagnitude(std::declval<VType>()))>
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	__host__ void renormalize(size_t arr_size, PType new_norm);
	
	//--------------------------------------------VEC GENERATORS : cuVEC_generate.cuh

	//linear : use interpolation to set values in this VEC based on projected distance between position1 and position2 and given fixed end values.
	__host__ bool generate_linear(cuReal3 new_h, cuRect new_rect, cuReal3 position1, VType value1, cuReal3 position2, VType value2);

	//similar to generate_linear except new dimensions not set
	__host__ void set_linear(cuReal3 position1, VType value1, cuReal3 position2, VType value2);

	//--------------------------------------------GETTERS : cuVEC_aux.h

	__device__ cuSZ3 size(void)  const { return n; }
	__device__ size_t linear_size(void)  const { return n.dim(); }
	__host__ size_t linear_size_cpu(void)  const { return get_gpu_value(n).dim(); }

	__host__ cuSZ3 size_cpu(void) { return get_gpu_value(n); }
	__host__ cuReal3 cellsize_cpu(void) { return get_gpu_value(h); }

	//from cell index return cell center coordinates (relative to start of rectangle)
	__device__ cuReal3 cellidx_to_position(int idx)  const;
	__host__ cuReal3 cellidx_to_position_cpu(int idx);
	
	//from cell index return cell center coordinates (relative to start of rectangle)
	__device__ cuReal3 cellidx_to_position(const cuINT3& ijk)  const;
	__host__ cuReal3 cellidx_to_position_cpu(cuINT3 ijk);

	//return cell index from relative position : the inverse of cellidx_to_position
	__device__ int position_to_cellidx(const cuReal3& position) const { return cu_floor_epsilon(position.x / h.x) + cu_floor_epsilon(position.y / h.y) * n.x + cu_floor_epsilon(position.z / h.z) * n.x*n.y; }

	//get index of cell which contains position (absolute value, not relative to start of rectangle), capped to mesh size
	__device__ cuINT3 cellidx_from_position(const cuReal3& absolute_position)  const;
	__host__ cuINT3 cellidx_from_position_cpu(cuReal3 absolute_position);

	//get cell rectangle (absolute values, not relative to start of mesh rectangle) for cell with index ijk
	__device__ cuRect get_cellrect(const cuINT3& ijk)  const;
	__host__ cuRect get_cellrect_cpu(cuINT3 ijk);
	
	//get_cellrect using single index.
	__device__ cuRect get_cellrect(int idx)  const;
	__host__ cuRect get_cellrect_cpu(int idx);

	//extract box of cells intersecting with the given rectangle (rectangle is in absolute coordinates). Cells in box : from and including start, up to but not including end; Limited to cuVEC sizes.
	__device__ cuBox box_from_rect_max(const cuRect& rectangle) const;
	__host__ cuBox box_from_rect_max_cpu(cuRect rectangle);
	
	//extract box of cells completely included in the given rectangle (rectangle is in absolute coordinates).
	__device__ cuBox box_from_rect_min(const cuRect& rectangle)  const;
	__host__ cuBox box_from_rect_min_cpu(cuRect rectangle);
	
	//count cells which don't have a null value set : i.e. non-empty; set result in aux_integer
	__host__ void count_nonempty_cells(size_t arr_size);
	//after using the above method, call this to get the number of non-empty points : e.g. count them just before launching a kernel, then inside the kernel it is available in aux_integer
	__device__ size_t get_aux_integer(void) const { return aux_integer; }

	//--------------------------------------------ARITHMETIC OPERATIONS ON ENTIRE VEC : cuVEC_arith.cuh

	//add to this vec the values in add_this : must have same size : size
	__host__ void add_values(size_t size, cu_obj<cuVEC<VType>>& add_this);

	//subtract from this vec the values in sub_this : must have same size : size
	__host__ void sub_values(size_t size, cu_obj<cuVEC<VType>>& sub_this);

	//--------------------------------------------AVERAGING OPERATIONS : cuVEC_avg.h (and cuVEC_avg.cuh)

	//average in a box (which should be contained in the cuVEC dimensions)
	//Launch it with arr_size = n.dim() : quicker to pass in this value rather than get it internally using get_gpu_value(n).dim()
	__host__ VType average(size_t arr_size, cuBox box);
	//average over given rectangle (relative to this cuVEC's rect)
	__host__ VType average(size_t arr_size, cuRect rectangle = cuRect());

	//even though cuVEC doesn't hold a shape we might want to obtain averages by excluding zero-value cells
	__host__ VType average_nonempty(size_t arr_size, cuBox box);
	__host__ VType average_nonempty(size_t arr_size, cuRect rectangle = cuRect());

	//smoother : obtain a weighted average value at coord, over a stencil of given size. All dimension units are same as h and rect. Include values from all cells which intersect the stencil.
	///the coord is taken as the centre value and is relative to the mesh rectangle start coordinate which might not be 0,0,0 : i.e. not an absolute value.
	//the weights vary linearly with distance from coord
	__device__ VType weighted_average(const cuReal3& coord, const cuReal3& stencil);
	
	//ijk is the cell index in a mesh with cellsize cs and same rect as this cuVEC; if cs is same as h then just read the value at ijk - much faster! If not then get the usual weighted average.
	__device__ VType weighted_average(const cuINT3& ijk, const cuReal3& cs);

	//full average in given rectangle (relative coordinates).
	__device__ VType average(const cuRect& rectangle);

	//average in given rectangle (relative coordinates), excluding zero points (assumed empty).
	__device__ VType average_nonempty(const cuRect& rectangle);

	//--------------------------------------------NUMERICAL PROPERTIES : cuVEC_nprops.cuhj

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

	//--------------------------------------------MESH TRANSFER : cuVEC_MeshTransfer.h

	//SINGLE INPUT, SINGLE OUTPUT

	//copy pre-calculated transfer info from cpu memory. return false if not enough memory to copy
	template <typename cpuVEC>
	__host__ bool copy_transfer_info(cu_arr<cuVEC<VType>>& mesh_in_arr, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuVEC& vec) 
	{ 
		return transfer.copy_transfer_info(mesh_in_arr, mesh_out_arr, vec); 
	}

	//MULTIPLE INPUTS, SINGLE OUTPUT

	//copy pre-calculated transfer info from cpu memory. return false if not enough memory to copy
	//mesh_in and mesh_in2 vectors must have same sizes
	//All VECs in mesh_in should be non-empty
	//Some VECs in mesh_in2 allowed to be non-empty (in this case single input is used), but otherwise should have exactly same dimensions as the corresponding VECs in mesh_in
	template <typename cpuVEC>
	__host__ bool copy_transfer_info_averagedinputs(cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<VType>>& mesh_in_arr2, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuVEC& vec) 
	{ 
		return transfer.copy_transfer_info_averagedinputs(mesh_in_arr1, mesh_in_arr2, mesh_out_arr, vec);
	}

	template <typename cpuVEC>
	__host__ bool copy_transfer_info_multipliedinputs(cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<cuBReal>>& mesh_in_arr2_real, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuVEC& vec)
	{
		return transfer.copy_transfer_info_multipliedinputs(mesh_in_arr1, mesh_in_arr2_real, mesh_out_arr, vec);
	}

	//MULTIPLE INPUTS, MULTIPLE OUTPUT

	//copy pre-calculated transfer info from cpu memory. return false if not enough memory to copy
	//mesh_in and mesh_in2 vectors must have same sizes; same as mesh_out, mesh_out2
	//All VECs in mesh_in and mesh_out should be non-empty
	//Some VECs in mesh_in2 and mesh_out2 allowed to be non-empty (in this single input/output is used), but otherwise should have exactly same dimensions as the corresponding VECs in mesh_in, mesh_out
	//Also if a VEC in mesh_in2 is non-empty the corresponding VEC in mesh_out2 should also be non-empty.
	template <typename cpuVEC>
	__host__ bool copy_transfer_info_averagedinputs_duplicatedoutputs(cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<VType>>& mesh_in_arr2, cu_arr<cuVEC<VType>>& mesh_out_arr1, cu_arr<cuVEC<VType>>& mesh_out_arr2, cpuVEC& vec)
	{
		return transfer.copy_transfer_info_averagedinputs_duplicatedoutputs(mesh_in_arr1, mesh_in_arr2, mesh_out_arr1, mesh_out_arr2, vec);
	}

	//SINGLE INPUT, SINGLE OUTPUT

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in(size_t size, size_t size_transfer) 
	{ 
		//first zero smesh quantity as we'll be adding in values from meshes in mesh_in
		set(size, VType());
		transfer.transfer_in(size_transfer, quantity); 
	}

	//transfer to output meshes. Pass in size_transfer (transfer_info_size) and number of output meshes if you want to zero the output meshes first (leave this default zero not to clear output meshes first)
	void transfer_out(size_t size_transfer, int mesh_out_num = 0) 
	{ 
		transfer.transfer_out(size_transfer, quantity, mesh_out_num); 
	}

	//AVERAGED INPUT

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in_averaged(size_t size, size_t size_transfer)
	{
		//first zero smesh quantity as we'll be adding in values from meshes in mesh_in
		set(size, VType());
		transfer.transfer_in_averaged(size_transfer, quantity);
	}

	//MULTIPLIED INPUTS

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in_multiplied(size_t size, size_t size_transfer)
	{
		//first zero smesh quantity as we'll be adding in values from meshes in mesh_in
		set(size, VType());
		transfer.transfer_in_multiplied(size_transfer, quantity);
	}

	//DUPLICATED OUTPUT

	//transfer to output meshes. Pass in size_transfer (transfer_info_size) and number of output meshes if you want to zero the output meshes first (leave this default zero not to clear output meshes first)
	void transfer_out_duplicated(size_t size_transfer, int mesh_out_num = 0)
	{
		transfer.transfer_out_duplicated(size_transfer, quantity, mesh_out_num);
	}
};
