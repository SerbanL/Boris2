#include "stdafx.h"
#include "DataProcessing.h"

DPArrays::DPArrays(int _max_arrays) :
	max_arrays(_max_arrays)
{
	dpA.resize(max_arrays);
}

//------------------------------------------------------------------------------------------ load_arrays

BError DPArrays::load_arrays(std::string fileName, std::vector<int> all_indexes, int* prows_read)
{
	BError error(__FUNCTION__);

	*prows_read = 0;

	//vector of indexes : first half are file data column indexes, second half are dp arrays indexes

	//must have an even number of indexes
	if (all_indexes.size() % 2) return error(BERROR_INCORRECTARRAYS);

	//check all indexes are valid
	if(!GoodArrays(all_indexes)) return error(BERROR_INCORRECTARRAYS);

	//load columns from data file
	std::vector<std::vector<double>> data_cols;
	ReadDataColumns(fileName, "\t", data_cols, subvec(all_indexes, 0, (int)all_indexes.size() / 2));

	//save columns in respective dp arrays
	for (int idx = 0; idx < (int)data_cols.size(); idx++) {

		dpA[all_indexes[all_indexes.size() / 2 + idx]] = data_cols[idx];
	}

	if (!data_cols.size()) return error(BERROR_COULDNOTLOADFILE);

	*prows_read = (int)data_cols[0].size();

	return error;
}

//------------------------------------------------------------------------------------------ set

//set 3 consecutive arrays by extracting from DBL3 in vector
BError DPArrays::set_arrays(int arr_idx, std::vector<DBL3>& data)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(arr_idx, arr_idx + 1, arr_idx + 2)) return error(BERROR_INCORRECTARRAYS);

	if (!resize(arr_idx, data.size())) return error(BERROR_OUTOFMEMORY_NCRIT);
	if (!resize(arr_idx + 1, data.size())) return error(BERROR_OUTOFMEMORY_NCRIT);
	if (!resize(arr_idx + 2, data.size())) return error(BERROR_OUTOFMEMORY_NCRIT);

#pragma omp parallel for
	for (int idx = 0; idx < data.size(); idx++) {

		dpA[arr_idx][idx] = data[idx].x;
		dpA[arr_idx + 1][idx] = data[idx].y;
		dpA[arr_idx + 2][idx] = data[idx].z;
	}

	return error;
}

//------------------------------------------------------------------------------------------ save_arrays

BError DPArrays::save_arrays(std::string fileName, std::vector<int> all_indexes, bool append)
{
	BError error(__FUNCTION__);

	//check all indexes are valid
	if (!GoodArrays(all_indexes)) return error(BERROR_INCORRECTARRAYS);

	if (!SaveDataColumns(fileName, '\t', dpA, all_indexes, append)) return error(BERROR_COULDNOTSAVEFILE);

	return error;
}

//save dp array as a single tab-spaced row
BError DPArrays::save_array_transposed(std::string fileName, int dp_index, bool append)
{
	BError error(__FUNCTION__);

	//check all indexes are valid
	if (!GoodArrays(dp_index)) return error(BERROR_INCORRECTARRAYS);

	if (!SaveDataRow(fileName, '\t', dpA[dp_index], append)) return error(BERROR_COULDNOTSAVEFILE);

	return error;
}

//--------------------- testing

BError DPArrays::dump_tdep(SuperMesh *pSMesh, std::string meshName, std::string paramName, double max_temperature, int dp_arr)
{
	BError error(__FUNCTION__);

	if(!pSMesh->contains(meshName) || !(*pSMesh)[meshName]->contains_param(paramName)) return error(BERROR_INCORRECTNAME);

	if (!GoodArrays(dp_arr + 2)) return error(BERROR_INCORRECTARRAYS);

	PARAM_ paramID = (PARAM_)(*pSMesh)[meshName]->get_meshparam_id(paramName);

	if (!(*pSMesh)[meshName]->get_meshparam_tempscaling(paramID, max_temperature, dpA[dp_arr], dpA[dp_arr + 1], dpA[dp_arr + 2])) return error(BERROR_OPERATIONFAILED);

	return error;
}

//------------------------------------------------------------------------------------------ get_profile

BError DPArrays::get_profile(DBL3 start, DBL3 end, SuperMesh *pSMesh, int arr_idx)
{
	BError error(__FUNCTION__);

	//allow space quantities (4 consecutive dp arrays used, first for distance along line, next for physical quantities)
	if (!GoodArrays(arr_idx, arr_idx + 1, arr_idx + 2, arr_idx + 3)) return error(BERROR_INCORRECTARRAYS);

	//physical quantity displayed in the current mesh (as we travel along the line this is updated as required)
	PhysQ *pcurrPhysQ = nullptr;

	//start and end points must be contained in the supermesh rectangle
	//Rect sMeshRect = pSMesh->GetSMeshRect();
	//if (!sMeshRect.contains(start) || !sMeshRect.contains(end)) return error(BERROR_INCORRECTVALUE);

	//point on line to get value at
	DBL3 point = start;

	//get line direction unit vector (making sure line has non-zero length)
	DBL3 direction = (end - start);
	double distance = get_distance(start, end);

	if (IsZ(distance)) return error(BERROR_INCORRECTVALUE);
	else direction /= distance;

	//length of cell along the line direction for the current mesh (calculated from a corner). 
	//Start from a reasonable non-zero value, just in case the starting point is not an any mesh, so cannot calculate a proper value to start with.
	double length = distance / 100;

	dpA[arr_idx].resize(0); 
	dpA[arr_idx + 1].resize(0);  
	dpA[arr_idx + 2].resize(0); 
	dpA[arr_idx + 3].resize(0);

	//travel along the line to extract profile
	while (point <= end) {

		//is the current point still in current mesh?
		if (pcurrPhysQ == nullptr || !pcurrPhysQ->rectangle().contains(point)) {

			//No, so need to update pcurrPhysQ
			if (pcurrPhysQ) delete pcurrPhysQ;
			pcurrPhysQ = nullptr;

			//get mesh containing this point, if any, unless super-mesh display is enabled
			if (pSMesh->GetDisplayedPhysicalQuantity() == MESHDISPLAY_NONE) {

				MeshBase* pcurrMesh = (*pSMesh)[point];
				if (pcurrMesh) {

					//found a mesh. Update pcurrPhysQ.
					pcurrPhysQ = new PhysQ(pcurrMesh->FetchOnScreenPhysicalQuantity());

					//update length for current mesh
					DBL3 h = pcurrPhysQ->cellsize();
					DBL3 hx = DBL3(h.x, 0, 0); DBL3 hy = DBL3(0, h.y, 0); DBL3 hz = DBL3(0, 0, h.z);

					//set displacement length (remember direction is normalized)
					length = DBL3(h.x * direction.x, h.y * direction.y, h.z * direction.z).norm();
				}
			}
			else {

				//super-mesh display is enabled so get pcurrPhysQ from it
				pcurrPhysQ = new PhysQ((pSMesh->FetchOnScreenPhysicalQuantity())[0]);

				//update length for current mesh display
				DBL3 h = pcurrPhysQ->cellsize();
				DBL3 hx = DBL3(h.x, 0, 0); DBL3 hy = DBL3(0, h.y, 0); DBL3 hz = DBL3(0, 0, h.z);

				//set displacement length (remember direction is normalized)
				length = DBL3(h.x * direction.x, h.y * direction.y, h.z * direction.z).norm();
			}
		}

		if (pcurrPhysQ && pcurrPhysQ->rectangle().contains(point)) {

			push_value(arr_idx, get_distance(start, point));

			DBL3 value;

			//get value at current point and store in dp array (can be vectorial - use 3 dp arrays - or scalar quantity - use 1 dp array)
			if (pcurrPhysQ->is_vectorial()) {

				value = pcurrPhysQ->get_vec_point(point);
			}
			else {

				value = DBL3(pcurrPhysQ->get_sca_point(point), 0, 0);
			}

			push_value(arr_idx + 1, value);
		}

		//next point on line for given current length increment
		point += direction * length;
	}

	if (pcurrPhysQ) delete pcurrPhysQ;

	return error;
}

//calculate topological charge in M and given rect, using equation Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
BError DPArrays::get_topological_charge(VEC_VC<DBL3>& M, double x, double y, double radius, double* pQ)
{
	BError error(__FUNCTION__);

	double Q = 0.0;

#pragma omp parallel for reduction(+:Q)
	for (int idx = 0; idx < M.linear_size(); idx++) {

		if (M.is_not_empty(idx)) {

			DBL3 pos = M.cellidx_to_position(idx);

			if (get_distance(DBL2(pos.x, pos.y), DBL2(x, y)) < radius) {

				double M_mag = M[idx].norm();

				DBL33 M_grad = M.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / M_mag;
				DBL3 dm_dy = M_grad.y / M_mag;

				Q += (M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y;
			}
		}
	}

	//divide by number of z cells (this is intended for 2D layers)
	*pQ = Q / (4 * PI * M.n.z);

	return error;
}

//count skyrmions in M and given rect, using equation Q = Integral(|m.(dm/dx x dm/dy)| dxdy) / 4PI
BError DPArrays::count_skyrmions(VEC_VC<DBL3>& M, double x, double y, double radius, double* pQ)
{
	BError error(__FUNCTION__);

	double Q = 0.0;

#pragma omp parallel for reduction(+:Q)
	for (int idx = 0; idx < M.linear_size(); idx++) {

		if (M.is_not_empty(idx)) {

			DBL3 pos = M.cellidx_to_position(idx);

			if (get_distance(DBL2(pos.x, pos.y), DBL2(x, y)) < radius) {

				double M_mag = M[idx].norm();

				DBL33 M_grad = M.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / M_mag;
				DBL3 dm_dy = M_grad.y / M_mag;

				Q += fabs((M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y);
			}
		}
	}

	//divide by number of z cells (this is intended for 2D layers)
	*pQ = Q / (4 * PI * M.n.z);

	return error;
}

//calculate histogram for |M1| using given parameters if the corresponding value on |M2| is within specified bounds of [M2val - deltaM2val, M2val + deltaM2val]
BError DPArrays::calculate_histogram2(VEC_VC<DBL3>& M1, VEC_VC<DBL3>& M2, int dp_x, int dp_y, int num_bins, double min, double max, double M2val, double deltaM2val)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y)) return error(BERROR_INCORRECTARRAYS);
	if (M1.size() != M2.size()) return error(BERROR_INCORRECTARRAYS);

	if (num_bins < 2) {

		num_bins = 100;

		OmpReduction<double> omp_reduction;
		omp_reduction.new_minmax_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < M1.linear_size(); idx++) {

			omp_reduction.reduce_minmax(M1[idx].norm());
		}

		DBL2 minmax = omp_reduction.minmax();

		min = minmax.i;
		max = minmax.j;
	}

	double bin = (max - min) / (num_bins - 1);

	resize(dp_x, num_bins);
	resize(dp_y, num_bins);

#pragma omp parallel for
	for (int idx = 0; idx < num_bins; idx++) {

		dpA[dp_x][idx] = min + idx * bin;
		dpA[dp_y][idx] = 0;
	}

	for (int idx = 0; idx < M1.linear_size(); idx++) {

		if (M1.is_not_empty(idx) && M2.is_not_empty(idx)) {

			double M2norm = M2[idx].norm();
			if (M2norm < M2val - deltaM2val || M2norm > M2val + deltaM2val) continue;

			double value = M1[idx].norm();

			int bin_idx = floor((value - (min - bin / 2)) / bin);

			if (bin_idx >= 0 && bin_idx < num_bins) {

				dpA[dp_y][bin_idx] += 1.0 / M1.get_nonempty_cells();
			}
		}
	}

	return error;
}

//--------------------- dp array manipulation

//append data in dp_new at the end of dp_original
BError DPArrays::append_array(int dp_original, int dp_new)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_original, dp_new)) return error(BERROR_INCORRECTARRAYS);

	int original_size = dpA[dp_original].size();

	resize(dp_original, dpA[dp_original].size() + dpA[dp_new].size());

#pragma omp parallel for
	for (int idx = original_size; idx < dpA[dp_original].size(); idx++) {

		dpA[dp_original][idx] = dpA[dp_new][idx - original_size];
	}

	return error;
}

//Pick elements from dp_in using the skip value (1 by default) and set them in dp_out; e.g. with skip = 2 every 3rd data point is picked; skip = 1 picks every other point.
BError DPArrays::rarefy(int dp_in, int dp_out, int skip)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_in, dp_out) || skip < 1) return error(BERROR_INCORRECTARRAYS);

	resize(dp_out, ceil_epsilon((double)dpA[dp_in].size() / (skip + 1)));

#pragma omp parallel for
	for (int idx = 0; idx < dpA[dp_out].size(); idx++) {

		dpA[dp_out][idx] = dpA[dp_in][idx * (skip + 1)];
	}

	return error;
}

//From  dp_in extract a number of points (length) starting at start_index and place them in dp_out
BError DPArrays::extract(int dp_in, int dp_out, int start_index, int length)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_in, dp_out) || start_index < 0 || start_index >= dpA[dp_in].size()) return error(BERROR_INCORRECTARRAYS);

	//if length is negative then extract all the remaining points
	if (length < 0 || length > dpA[dp_in].size() - start_index) length = dpA[dp_in].size() - start_index;

	dpA[dp_out] = subvec(dpA[dp_in], start_index, start_index + length);

	return error;
}

//From dp_index erase a number of points (length) starting at start_index
BError DPArrays::erase(int dp_index, int start_index, int length)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_index) || start_index < 0 || start_index >= dpA[dp_index].size()) return error(BERROR_INCORRECTARRAYS);

	//if length is negative then extract all the remaining points
	if (length < 0 || length > dpA[dp_index].size() - start_index) length = dpA[dp_index].size() - start_index;

	dpA[dp_index].erase(dpA[dp_index].begin() + start_index, dpA[dp_index].begin() + start_index + length);

	return error;
}

//--------------------- data generation

//generate a sequence of points in dp_idx from start_value using increment
BError DPArrays::generate_sequence(int dp_idx, double start_value, double increment, int points)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_idx) || points < 1) return error(BERROR_INCORRECTARRAYS);

	resize(dp_idx, points);

#pragma omp parallel for
	for (int idx = 0; idx < points; idx++) {

		dpA[dp_idx][idx] = start_value + increment * idx;
	}

	return error;
}

//------------------------------------------------------------------------------------------ add, subtract, multiply, divide

//single source versions

BError DPArrays::add(int dp_source, int dp_dest, double value)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	resize(dp_dest, dpA[dp_source].size());

	#pragma omp parallel for
	for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

		dpA[dp_dest][idx] = dpA[dp_source][idx] + value;
	}

	return error;
}

BError DPArrays::subtract(int dp_source, int dp_dest, double value)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	resize(dp_dest, dpA[dp_source].size());

	#pragma omp parallel for
	for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

		dpA[dp_dest][idx] = dpA[dp_source][idx] - value;
	}

	return error;
}

BError DPArrays::divide(int dp_source, int dp_dest, double value)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	resize(dp_dest, dpA[dp_source].size());

	#pragma omp parallel for
	for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

		dpA[dp_dest][idx] = dpA[dp_source][idx] / value;
	}

	return error;
}

BError DPArrays::exponentiate(int dp_source, int dp_dest, double exponent)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	resize(dp_dest, dpA[dp_source].size());

#pragma omp parallel for
	for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

		dpA[dp_dest][idx] = pow(dpA[dp_source][idx], exponent);
	}

	return error;
}

BError DPArrays::multiply(int dp_source, int dp_dest, double value)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	resize(dp_dest, dpA[dp_source].size());

	#pragma omp parallel for
	for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

		dpA[dp_dest][idx] = dpA[dp_source][idx] * value;
	}

	return error;
}

BError DPArrays::dotproduct(int dp_x, int dp_y, int dp_z, DBL3 u, int dp_out)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y, dp_z, dp_out)) return error(BERROR_INCORRECTARRAYS);

	int vec_size = minimum(dpA[dp_x].size(), dpA[dp_y].size(), dpA[dp_z].size());
	resize(dp_out, vec_size);

#pragma omp parallel for
	for (int idx = 0; idx < vec_size; idx++) {

		dpA[dp_out][idx] = u * DBL3(dpA[dp_x][idx], dpA[dp_y][idx], dpA[dp_z][idx]);
	}

	return error;
}

//multiple sources versions

BError DPArrays::add_dparr(int dp_x1, int dp_x2, int dp_dest)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_x1, dp_x2, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	int vec_size = minimum(dpA[dp_x1].size(), dpA[dp_x2].size());

	resize(dp_dest, vec_size);

#pragma omp parallel for
	for (int idx = 0; idx < vec_size; idx++) {

		dpA[dp_dest][idx] = dpA[dp_x1][idx] + dpA[dp_x2][idx];
	}

	return error;
}

BError DPArrays::subtract_dparr(int dp_x1, int dp_x2, int dp_dest)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_x1, dp_x2, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	int vec_size = minimum(dpA[dp_x1].size(), dpA[dp_x2].size());

	resize(dp_dest, vec_size);

#pragma omp parallel for
	for (int idx = 0; idx < vec_size; idx++) {

		dpA[dp_dest][idx] = dpA[dp_x1][idx] - dpA[dp_x2][idx];
	}

	return error;
}

BError DPArrays::multiply_dparr(int dp_x1, int dp_x2, int dp_dest)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_x1, dp_x2, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	int vec_size = minimum(dpA[dp_x1].size(), dpA[dp_x2].size());

	resize(dp_dest, vec_size);

#pragma omp parallel for
	for (int idx = 0; idx < vec_size; idx++) {

		dpA[dp_dest][idx] = dpA[dp_x1][idx] * dpA[dp_x2][idx];
	}

	return error;
}

BError DPArrays::divide_dparr(int dp_x1, int dp_x2, int dp_dest)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_x1, dp_x2, dp_dest)) return error(BERROR_INCORRECTARRAYS);

	int vec_size = minimum(dpA[dp_x1].size(), dpA[dp_x2].size());

	resize(dp_dest, vec_size);

#pragma omp parallel for
	for (int idx = 0; idx < vec_size; idx++) {

		if(dpA[dp_x2][idx]) dpA[dp_dest][idx] = dpA[dp_x1][idx] / dpA[dp_x2][idx];
		else dpA[dp_dest][idx] = 0.0;
	}

	return error;
}

BError DPArrays::dotproduct_dparr(int dp_x1, int dp_x2, double *pvalue)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_x1, dp_x2)) return error(BERROR_INCORRECTARRAYS);

	int vec_size = minimum(dpA[dp_x1].size(), dpA[dp_x2].size());

	double value = 0.0;

#pragma omp parallel for reduction(+:value)
	for (int idx = 0; idx < vec_size; idx++) {

		value += dpA[dp_x1][idx] * dpA[dp_x2][idx];
	}

	*pvalue = value;

	return error;
}

//------------------------------------------------------------------------------------------ get_min_max, get_mean

BError DPArrays::get_min_max(int dp_source, DBL2* pmin_max_values, INT2* pmin_max_indexes)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source) || !dpA[dp_source].size()) return error(BERROR_INCORRECTARRAYS);

	auto min_max_pair = std::minmax_element(dpA[dp_source].begin(), dpA[dp_source].end());

	if(pmin_max_indexes)
		*pmin_max_indexes = INT2((int)distance(dpA[dp_source].begin(), min_max_pair.first), 
								 (int)distance(dpA[dp_source].begin(), min_max_pair.second));

	*pmin_max_values = DBL2(*(min_max_pair.first), *(min_max_pair.second));

	return error;
}

BError DPArrays::get_mean(int dp_source, DBL2* pmean_err, double exclusion_ratio)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source) || !dpA[dp_source].size()) return error(BERROR_INCORRECTARRAYS);

	//mean : sum of points / number of points
	double mean = sum_KahanNeumaier(dpA[dp_source]) / dpA[dp_source].size();
	double stdev_sq = 0;

	//second pass required to exclude points if included
	if (exclusion_ratio > 0 && mean != 0.0) {

		double mean2 = 0.0;
		int num_mean_points = 0;

#pragma omp parallel for reduction(+:mean2, num_mean_points)
		for (int idx = 0; idx < dpA[dp_source].size(); idx++) {

			double value = dpA[dp_source][idx];
			if (fabs(value / mean - 1.0) < exclusion_ratio) {

				mean2 += value;
				num_mean_points++;
			}
		}

		if (num_mean_points) {

			mean2 /= num_mean_points;

#pragma omp parallel for reduction(+:stdev_sq)
			for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

				double value = dpA[dp_source][idx];
				if (fabs(value / mean - 1.0) < exclusion_ratio) {

					stdev_sq += pow(value - mean2, 2);
				}
			}

			stdev_sq /= num_mean_points;
		}

		mean = mean2;
	}
	else {

		//unweighted standard deviation squared : sum of squares of distances between points and mean, divided by number of points
#pragma omp parallel for reduction(+:stdev_sq)
		for (int idx = 0; idx < (int)dpA[dp_source].size(); idx++) {

			stdev_sq += pow(dpA[dp_source][idx] - mean, 2);
		}

		stdev_sq /= dpA[dp_source].size();
	}

	*pmean_err = DBL2(mean, sqrt(stdev_sq));

	return error;
}

BError DPArrays::get_sum(int dp_source, double* psum)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source) || !dpA[dp_source].size()) return error(BERROR_INCORRECTARRAYS);

	double sum = 0.0;

#pragma omp parallel for reduction(+:sum)
	for (int idx = 0; idx < dpA[dp_source].size(); idx++) {

		sum += dpA[dp_source][idx];
	}
	
	*psum = sum;

	return error;
}

//find standard deviation on the mean in chunks: i.e. obtain it from every chunk number of values in dp_source, then average them to obtain final std
//in the limit chunk = dp_source number of points (set chunk = 0) for this special case, the std is the usual std on the mean
BError DPArrays::get_chunkedstd(int dp_source, int chunk, double* std)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source) || !dpA[dp_source].size()) return error(BERROR_INCORRECTARRAYS);

	if (chunk <= 0 || chunk > dpA[dp_source].size()) chunk = dpA[dp_source].size();

	double cstd = 0.0;
	int num_chunks = dpA[dp_source].size() / chunk;
	if (num_chunks == 0) num_chunks = 1;

#pragma omp parallel for reduction(+:cstd)
	for (int idx = 0; idx < dpA[dp_source].size(); idx += chunk) {

		double av = 0.0, avsq = 0.0;

		if (idx + chunk <= dpA[dp_source].size()) {

			for (int idx_chunk = 0; idx_chunk < chunk; idx_chunk++) {

				double value = dpA[dp_source][idx + idx_chunk];
				av += value / chunk;
				avsq += value * value / chunk;
			}

			cstd += (avsq - av * av) / num_chunks;
		}
	}

	if ((dpA[dp_source].size() % chunk) && num_chunks > 1) *std = sqrt(cstd * num_chunks / (num_chunks - 1));
	else *std = sqrt(cstd);

	return error;
}

//get amplitude value (max - min)/2 every pointsPeriod
BError DPArrays::get_amplitude(int dp_source, int pointsPeriod, double* pamplitude)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_source) || dpA[dp_source].size() < pointsPeriod) return error(BERROR_INCORRECTARRAYS);

	double maxval = 0.0;
	double minval = 0.0;

	double amplitude = 0;

	for (int idx = 0; idx < dpA[dp_source].size() / pointsPeriod; idx++) {

		maxval = dpA[dp_source][idx*pointsPeriod];
		minval = dpA[dp_source][idx*pointsPeriod];

		for (int idx_cycle = idx * pointsPeriod; idx_cycle < (idx + 1)*pointsPeriod; idx_cycle++) {

			if (dpA[dp_source][idx_cycle] > maxval) maxval = dpA[dp_source][idx_cycle];
			if (dpA[dp_source][idx_cycle] < minval) minval = dpA[dp_source][idx_cycle];
		}

		amplitude = (amplitude > (maxval - minval) / 2 ? amplitude : (maxval - minval) / 2);
	}

	*pamplitude = amplitude;

	return error;
}

//------------------------------------------------------------------------------------------ linreg

BError DPArrays::linreg(int dp_x, int dp_y, int dp_z, int dp_out, DBL2* pgradient, DBL2* pintercept)
{
	BError error(__FUNCTION__);

	//Formulas for unweighted linear regression :
	// g = (<xy> - <x><y>) / (<x^2> - <x>^2)
	// c = <y> - g<x>

	//sig_g^2 = (1/(N-2)) * <(y-y_fit)^2> / (<x^2> - <x>^2)
	//sig_c^2 = sig_g^2 * <x^2>

	if (!GoodArrays_Unique(dp_x, dp_y)) return error(BERROR_INCORRECTARRAYS);
	
	if (dp_z < 0 || dp_out < 0) {

		//dp_z not used, so just a simple linear regression with no dp_out output
		std::pair<DBL2, DBL2> result = linear_regression(dpA[dp_x], dpA[dp_y]);

		*pgradient = result.first;
		*pintercept = result.second;
	}
	else {

		if (!GoodArrays_Unique(dp_x, dp_y, dp_z, dp_out, dp_out + 1, dp_out + 2, dp_out + 3, dp_out + 4)) return error(BERROR_INCORRECTARRAYS);

		int N = minimum(dpA[dp_x].size(), dpA[dp_y].size(), dpA[dp_z].size());
		if (N < 2) return error(BERROR_INCORRECTARRAYS);

		clear(dp_out, 5);

		int startIdx = 0, lastIdx = 1;
		double zValue = dpA[dp_z][0];

		while (lastIdx < N) {

			//test lastIdx
			if (IsE(dpA[dp_z][lastIdx], zValue)) {

				//same zValue so increment and test next...
				lastIdx++;
				
				//...unless we reached the end
				if (lastIdx < N) continue;
			}

			//not same zValue or lastIdx is equal to N. Do we have enough points for a linear regression?
			if (lastIdx - startIdx >= 2) {

				std::pair<DBL2, DBL2> result = linear_regression(dpA[dp_x], dpA[dp_y], startIdx, lastIdx - 1);

				*pgradient = result.first;
				*pintercept = result.second;

				dpA[dp_out + 0].push_back(zValue);
				dpA[dp_out + 1].push_back(result.first.i);
				dpA[dp_out + 2].push_back(result.first.j);
				dpA[dp_out + 3].push_back(result.second.i);
				dpA[dp_out + 4].push_back(result.second.j);
			}

			//reached the end
			if (lastIdx == N) break;

			//Not the end. start getting indexes for next regression
			startIdx = lastIdx;
			zValue = dpA[dp_z][startIdx];
			lastIdx++;
		}
	}

	return error;
}

//------------------------------------------------------------------------------------------ get_coercivity, get_remanence

BError DPArrays::get_coercivity(int dp_x, int dp_y, DBL3 *pHc_up, DBL3 *pHc_dn)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || !dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	auto find_up_switch = [&](int dp_x, int dp_y, int &idx) -> DBL3 {

		int N = (int)(dpA[dp_x].size() > dpA[dp_y].size() ? dpA[dp_y].size() : dpA[dp_x].size());

		DBL3 Hc = DBL3();

		while (idx < N) {

			//it is guaranteed that dpA[dp_y][idx - 1] < 0.
			if (dpA[dp_y][idx] >= 0 && dpA[dp_x][idx] > dpA[dp_x][idx - 1]) {

				//found switching
				
				//linear interpolation to find switching field
				double alpha = -dpA[dp_y][idx - 1] / (dpA[dp_y][idx] - dpA[dp_y][idx - 1]);
				Hc.x = dpA[dp_x][idx - 1] * (1 - alpha) + dpA[dp_x][idx] * alpha;

				//uncertainty range : Hc (-Hc.y, +Hc.z)
				Hc.y = Hc.x - dpA[dp_x][idx - 1];
				Hc.z = dpA[dp_x][idx] - Hc.x;

				break;
			}
			else idx++;
		}

		return Hc;
	};

	auto find_dn_switch = [&](int dp_x, int dp_y, int &idx) -> DBL3 {

		int N = (int)(dpA[dp_x].size() > dpA[dp_y].size() ? dpA[dp_y].size() : dpA[dp_x].size());

		DBL3 Hc = DBL3();

		while (idx < N) {

			//it is guaranteed that dpA[dp_y][idx - 1] >= 0.
			if (dpA[dp_y][idx] < 0 && dpA[dp_x][idx] < dpA[dp_x][idx - 1]) {

				//found switching

				//linear interpolation to find switching field
				double alpha = dpA[dp_y][idx - 1] / (dpA[dp_y][idx - 1] - dpA[dp_y][idx]);
				Hc.x = dpA[dp_x][idx - 1] * (1 - alpha) + dpA[dp_x][idx] * alpha;

				//uncertainty range : Hc (-Hc.y, +Hc.z)
				Hc.y = Hc.x - dpA[dp_x][idx];
				Hc.z = dpA[dp_x][idx - 1] - Hc.x;

				break;
			}
			else idx++;
		}

		return Hc;
	};

	int idx = 0;

	if (dpA[dp_y][0] < 0) {

		*pHc_up = find_up_switch(dp_x, dp_y, idx);
		*pHc_dn = find_dn_switch(dp_x, dp_y, idx);
	}
	else {

		*pHc_dn = find_dn_switch(dp_x, dp_y, idx);
		*pHc_up = find_up_switch(dp_x, dp_y, idx);
	}

	return error;
}

BError DPArrays::get_remanence(int dp_x, int dp_y, double *pMr_up, double *pMr_dn)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || !dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	bool found_first = false;

	for (int idx = 0; idx < (int)(dpA[dp_x].size() > dpA[dp_y].size() ? dpA[dp_y].size() : dpA[dp_x].size()); idx++) {

		if (IsZ(dpA[dp_x][idx])) {

			if (dpA[dp_y][idx] < 0) *pMr_up = dpA[dp_y][idx];
			else *pMr_dn = dpA[dp_y][idx];

			if (found_first) break;
			else found_first = true;
		}
	}

	return error;
}

//For a hysteresis loop with only one branch continue it by constructing the other direction branch (invert both x and y data and add it in continuation)
BError DPArrays::complete_hysteresis(int dp_x, int dp_y)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || !dpA[dp_y].size() || dpA[dp_x].size() != dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	std::vector<double> dp_x_new(dpA[dp_x].size()), dp_y_new(dpA[dp_y].size());

#pragma omp parallel for
	for (int idx = 0; idx < dpA[dp_x].size(); idx++) {

		dp_x_new[idx] = dpA[dp_x][idx] * -1.0;
		dp_y_new[idx] = dpA[dp_y][idx] * -1.0;
	}

	JoinToVector(dpA[dp_x], dp_x_new);
	JoinToVector(dpA[dp_y], dp_y_new);

	return error;
}


//--------------------- histograms

//From input x - y data build a histogram of number of times x - y data crosses a given line (up or down). 
//The line varies between minimum and maximum of y data in given number of steps. 
//Output the line values in dp_lev with corresponding number of crossings in dp_counts.
BError DPArrays::crossings_histogram(int dp_x, int dp_y, int dp_lev, int dp_counts, int steps)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y, dp_lev, dp_counts) || !dpA[dp_y].size() || dpA[dp_x].size() != dpA[dp_y].size() || steps < 1) return error(BERROR_INCORRECTARRAYS);

	resize(dp_lev, steps + 1);
	resize(dp_counts, steps + 1);

	auto count_crossings = [&](double level) -> int {

		double last_value = dpA[dp_y][0];

		int num_crossings = 0;

		for (int idx = 1; idx < dpA[dp_x].size(); idx++) {

			double current_value = dpA[dp_y][idx];

			if ((last_value <= level && current_value >= level) ||
				(last_value >= level && current_value <= level)) num_crossings++;

			last_value = current_value;
		}

		return num_crossings;
	};

	DBL2 minmax;
	get_min_max(dp_y, &minmax);

#pragma omp parallel for
	for (int idx = 0; idx <= steps; idx++) {

		dpA[dp_lev][idx] = minmax.i + idx * (minmax.j - minmax.i) / steps;
		dpA[dp_counts][idx] = count_crossings(dpA[dp_lev][idx]);
	}

	return error;
}

//From input x-y data build a histogram of average frequency the x-y data crosses a given line (up and down, separated). 
//The line varies between minimum and maximum of y data in given number of steps (100 by default). 
//Output the line values in dp_lev with corresponding crossgins frequencies in dp_freq_up and dp_freq_dn.
BError DPArrays::crossings_frequencies(int dp_x, int dp_y, int dp_lev, int dp_freq_up, int dp_freq_dn, int steps)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y, dp_lev, dp_freq_up, dp_freq_dn) || !dpA[dp_y].size() || dpA[dp_x].size() != dpA[dp_y].size() || steps < 1) return error(BERROR_INCORRECTARRAYS);

	resize(dp_lev, steps + 1);
	resize(dp_freq_up, steps + 1);
	resize(dp_freq_dn, steps + 1);

	auto crossings_frequencies = [&](double level) -> DBL2 {

		double last_value = dpA[dp_y][0];
		
		double last_crossing_up = dpA[dp_x][0];
		double last_crossing_dn = dpA[dp_x][0];

		double average_period_up = 0.0;
		double average_period_dn = 0.0;

		int num_crossings_up = 0;
		int num_crossings_dn = 0;

		for (int idx = 1; idx < dpA[dp_x].size(); idx++) {

			double current_value = dpA[dp_y][idx];

			if (last_value <= level && current_value >= level) {

				//crossing up
				double current_crossing_up = dpA[dp_x][idx];

				average_period_up += (current_crossing_up - last_crossing_up);
				num_crossings_up++;

				last_crossing_up = current_crossing_up;
			}

			if (last_value >= level && current_value <= level) {

				//crossing down
				double current_crossing_dn = dpA[dp_x][idx];

				average_period_dn += (current_crossing_dn - last_crossing_dn);
				num_crossings_dn++;

				last_crossing_dn = current_crossing_dn;
			}

			last_value = current_value;
		}

		DBL2 frequencies;

		if (num_crossings_up && average_period_up) frequencies.i = 1.0 / (average_period_up / num_crossings_up);
		if (num_crossings_dn && average_period_dn) frequencies.j = 1.0 / (average_period_dn / num_crossings_dn);

		return frequencies;
	};

	DBL2 minmax;
	get_min_max(dp_y, &minmax);

#pragma omp parallel for
	for (int idx = 0; idx <= steps; idx++) {

		dpA[dp_lev][idx] = minmax.i + idx * (minmax.j - minmax.i) / steps;
		DBL2 frequencies = crossings_frequencies(dpA[dp_lev][idx]);

		dpA[dp_freq_up][idx] = frequencies.i;
		dpA[dp_freq_dn][idx] = frequencies.j;
	}

	return error;
}

//From input x-y data build a histogram of average frequency of peaks in the x-y data in bands given by the number of steps. 
//The bands vary between minimum and maximum of y data in given number of steps (100 by default). 
//Output the line values in dp_level with corresponding peak frequencies in dp_freq.
BError  DPArrays::peaks_frequencies(int dp_x, int dp_y, int dp_lev, int dp_freq, int steps)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y, dp_lev, dp_freq) || !dpA[dp_y].size() || dpA[dp_x].size() != dpA[dp_y].size() || steps < 1) return error(BERROR_INCORRECTARRAYS);

	resize(dp_lev, steps);
	resize(dp_freq, steps);

	auto peaks_frequency = [&](double band_start, double band_end) -> double {

		double last_value = dpA[dp_y][0];

		double last_peak = dpA[dp_x][0];

		double average_period = 0.0;

		int num_peaks = 0;

		for (int idx = 1; idx < dpA[dp_x].size() - 1; idx++) {

			double current_value = dpA[dp_y][idx];
			double next_value = dpA[dp_y][idx + 1];

			if (band_start <= current_value && current_value <= band_end) {

				if (last_value < current_value && next_value < current_value) {

					double current_peak = dpA[dp_x][idx];

					average_period += (current_peak - last_peak);
					num_peaks++;

					last_peak = current_peak;
				}
			}

			last_value = current_value;
		}

		if (num_peaks && average_period) return 1.0 / (average_period / num_peaks);
		else return 0.0;
	};

	DBL2 minmax;
	get_min_max(dp_y, &minmax);

#pragma omp parallel for
	for (int idx = 0; idx < steps; idx++) {

		dpA[dp_lev][idx] = minmax.i + (idx + 1) * (minmax.j - minmax.i) / steps;
		dpA[dp_freq][idx] = peaks_frequency(dpA[dp_lev][idx] - (minmax.j - minmax.i) / steps, dpA[dp_lev][idx]);
	}

	return error;
}

//--------------------- curve fitting

//fit f(x) = y0 + S dH / (4(x-H0)^2 + dH^2). Return fitting parameters with their standard deviations.
BError DPArrays::fit_lorentz(int dp_x, int dp_y, DBL2 *pS, DBL2 *pH0, DBL2 *pdH, DBL2 *py0)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || dpA[dp_x].size() != dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	CurveFitting curve_fit;

	std::vector<double> params, params_std;

	std::vector<DBL2> xy_data;

	xy_data.resize(dpA[dp_x].size());

#pragma omp parallel for
	for (int idx = 0; idx < xy_data.size(); idx++) {

		xy_data[idx] = DBL2(dpA[dp_x][idx], dpA[dp_y][idx]);
	}

	curve_fit.FitLorentz_LMA(xy_data, params, params_std);

	if (params.size() >= 4 && params_std.size() >= 4) {

		*pS = DBL2(params[0], params_std[0]);
		*pH0 = DBL2(params[1], params_std[1]);
		*pdH = DBL2(params[2], params_std[2]);
		*py0 = DBL2(params[3], params_std[3]);
	}
	else return error(BERROR_INCORRECTARRAYS);

	return error;
}

//fit f(x) = y0 + S (dH + A*(x-H0)) / (4(x-H0)^2 + dH^2). Return fitting parameters with their standard deviations.
BError DPArrays::fit_lorentz2(int dp_x, int dp_y, DBL2 *pS, DBL2 *pA, DBL2 *pH0, DBL2 *pdH, DBL2 *py0)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || dpA[dp_x].size() != dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	CurveFitting curve_fit;

	std::vector<double> params, params_std;

	std::vector<DBL2> xy_data;

	xy_data.resize(dpA[dp_x].size());

#pragma omp parallel for
	for (int idx = 0; idx < xy_data.size(); idx++) {

		xy_data[idx] = DBL2(dpA[dp_x][idx], dpA[dp_y][idx]);
	}

	curve_fit.FitLorentzSA_LMA(xy_data, params, params_std);

	if (params.size() >= 5 && params_std.size() >= 5) {

		*pS = DBL2(params[0], params_std[0]);
		*pA = DBL2(params[1], params_std[1]);
		*pH0 = DBL2(params[2], params_std[2]);
		*pdH = DBL2(params[3], params_std[3]);
		*py0 = DBL2(params[4], params_std[4]);
	}
	else return error(BERROR_INCORRECTARRAYS);

	return error;
}

//fit Mz(x) = Ms * cos(2*arctan(sinh(R/w)/sinh((x-x0)/w))). Return fitting parameters with their standard deviations.
BError DPArrays::fit_skyrmion(int dp_x, int dp_y, DBL2 *pR, DBL2 *px0, DBL2 *pMs, DBL2 *pw)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || dpA[dp_x].size() != dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	CurveFitting curve_fit;

	std::vector<double> params, params_std;

	std::vector<DBL2> xy_data;

	xy_data.resize(dpA[dp_x].size());

#pragma omp parallel for
	for (int idx = 0; idx < xy_data.size(); idx++) {

		xy_data[idx] = DBL2(dpA[dp_x][idx], dpA[dp_y][idx]);
	}

	curve_fit.FitSkyrmion_LMA(xy_data, params, params_std);

	if (params.size() >= 4 && params_std.size() >= 4) {

		*pR = DBL2(params[0], params_std[0]);
		*px0 = DBL2(params[1], params_std[1]);
		*pMs = DBL2(params[2], params_std[2]);
		*pw = DBL2(params[3], params_std[3]);
	}
	else return error(BERROR_INCORRECTARRAYS);

	return error;
}

//fit M(x) = A * tanh((x - x0) / (D / 2)). Return fitting parameters with their standard deviations.
BError DPArrays::fit_domainwall(int dp_x, int dp_y, DBL2 *pA, DBL2 *px0, DBL2 *pD)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_x, dp_y) || dpA[dp_x].size() != dpA[dp_y].size()) return error(BERROR_INCORRECTARRAYS);

	CurveFitting curve_fit;

	std::vector<double> params, params_std;

	std::vector<DBL2> xy_data;

	xy_data.resize(dpA[dp_x].size());

#pragma omp parallel for
	for (int idx = 0; idx < xy_data.size(); idx++) {

		xy_data[idx] = DBL2(dpA[dp_x][idx], dpA[dp_y][idx]);
	}

	curve_fit.FitDW_LMA(xy_data, params, params_std);

	if (params.size() >= 3 && params_std.size() >= 3) {

		*pA = DBL2(params[0], params_std[0]);
		*px0 = DBL2(params[1], params_std[1]);
		*pD = DBL2(params[2], params_std[2]);
	}
	else return error(BERROR_INCORRECTARRAYS);

	return error;
}

//fit STT function : 
BError DPArrays::fit_stt(VEC<DBL3>& T, VEC_VC<DBL3>& M, VEC_VC<DBL3>& J, Rect rectangle, DBL2* pP, DBL2 *pbeta, double *pRsq)
{
	BError error(__FUNCTION__);

	Box cells_box = M.box_from_rect_max(rectangle + M.rect.s);
	INT3 box_sizes = cells_box.size();
	int points = box_sizes.dim();

	std::vector<DBL2> xy_data_adiabatic_y(points, DBL2()), xy_data_nonadiabatic_y(points, DBL2());
	std::vector<double> adiabatic_x(points, 0.0), nonadiabatic_x(points, 0.0);

	//to extract P and beta we separate the orthogonal adiabatic and non-adiabatic terms and fit them separately
	//Thus first fit (J.del)m * Ts = P * MUB_E * [(J.del)m]^2 -> this gives P
	//Next, using the obtained P, fit (m^(J.del)m) * Ts = -P * MUB_E * beta * [m^(J.del)m]^2 -> this gives beta

#pragma omp parallel for
	for (int i = cells_box.s.i; i < cells_box.e.i; i++) {
		for (int j = cells_box.s.j; j < cells_box.e.j; j++) {
			for (int k = cells_box.s.k; k < cells_box.e.k; k++) {

				DBL3 position = M.cellidx_to_position(INT3(i, j, k));
				int idx_M = M.position_to_cellidx(position);
				int idx_J = J.position_to_cellidx(position);
				int idx = (i - cells_box.s.i) + (j - cells_box.s.j) * box_sizes.x + (k - cells_box.s.k) * box_sizes.x * box_sizes.y;

				if (M.is_not_empty(idx_M) && J.is_not_empty(idx_J)) {

					DBL3 m = normalize(M[idx_M]);
					DBL33 grad_m = normalize(M.grad_neu(idx_M), M[idx_M]);
					DBL3 J_dot_del_m = grad_m.x * J[idx_J].x + grad_m.y * J[idx_J].y + grad_m.z * J[idx_J].z;
					DBL3 m_cross_J_dot_del_m = m ^ J_dot_del_m;

					xy_data_adiabatic_y[idx] = DBL2(idx, J_dot_del_m * T[idx_M]);
					xy_data_nonadiabatic_y[idx] = DBL2(idx, m_cross_J_dot_del_m * T[idx_M]);

					adiabatic_x[idx] = J_dot_del_m * J_dot_del_m;
					nonadiabatic_x[idx] = m_cross_J_dot_del_m * m_cross_J_dot_del_m;
				}
			}
		}
	}

	CurveFitting curve_fit;

	std::vector<double> params, params_std;

	curve_fit.FitSTT_LMA(xy_data_adiabatic_y, xy_data_nonadiabatic_y, adiabatic_x, nonadiabatic_x, params, params_std, pRsq);

	if (params.size() >= 2 && params_std.size() >= 2) {

		*pP = DBL2(params[0], params_std[0]);
		*pbeta = DBL2(params[1], params_std[1]);
	}
	else return error(BERROR_INCORRECTARRAYS);

	return error;
}

//fit STT function using a stencil to obtain spatial dependence of P parameter
BError DPArrays::fit_stt_variation(VEC<DBL3>& T, VEC_VC<DBL3>& M, VEC_VC<DBL3>& J, VEC<double>& output, bool adiabatic, double absolute_error_threshold, double Rsq_threshold, double Torque_ratio_threshold, int stencil_size)
{
	BError error(__FUNCTION__);

	//to extract P and beta we separate the orthogonal adiabatic and non-adiabatic terms and fit them separately
	//Thus first fit (J.del)m * Ts = P * MUB_E * [(J.del)m]^2 -> this gives P
	//Next, using the obtained P, fit (m^(J.del)m) * Ts = -P * MUB_E * beta * [m^(J.del)m]^2 -> this gives beta

	INT3 n = M.n;
	int points = n.dim();

	if (!output.assign(M.h, M.rect, 0.0)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//1. Find maximum torque to use when calculating cutoff
	double T_max = 0.0;

	for (int idx = 0; idx < T.linear_size(); idx++) {

		if (T_max < T[idx].norm()) T_max = T[idx].norm();
	}

	//2. Calculate fitting vectors for the entire mesh space
	std::vector<double> xy_data_adiabatic_y(points, 0.0), xy_data_nonadiabatic_y(points, 0.0);
	std::vector<double> adiabatic_x(points, 0.0), nonadiabatic_x(points, 0.0);

#pragma omp parallel for
	for (int i = 0; i < n.i; i++) {
		for (int j = 0; j < n.j; j++) {
			for (int k = 0; k < n.k; k++) {

				DBL3 position = M.cellidx_to_position(INT3(i, j, k));
				int idx_M = M.position_to_cellidx(position);
				int idx_J = J.position_to_cellidx(position);
				int idx = i + j * n.i + k * n.i * n.j;

				if (M.is_not_empty(idx_M) && J.is_not_empty(idx_J)) {

					DBL3 m = normalize(M[idx_M]);
					DBL33 grad_m = normalize(M.grad_neu(idx_M), M[idx_M]);
					DBL3 J_dot_del_m = grad_m.x * J[idx_J].x + grad_m.y * J[idx_J].y + grad_m.z * J[idx_J].z;
					DBL3 m_cross_J_dot_del_m = m ^ J_dot_del_m;

					xy_data_adiabatic_y[idx] = J_dot_del_m * T[idx_M];
					xy_data_nonadiabatic_y[idx] = m_cross_J_dot_del_m * T[idx_M];

					adiabatic_x[idx] = J_dot_del_m * J_dot_del_m;
					nonadiabatic_x[idx] = m_cross_J_dot_del_m * m_cross_J_dot_del_m;
				}
			}
		}
	}

	//3. next use a nearest neighbor stencil centered around a point to fit for beta at that point

	//stencil with odd dimensions
	stencil_size += int((stencil_size % 2) == 0);

	INT3 n_stencil = INT3(stencil_size, stencil_size, 1);

	CurveFitting curve_fit;
	std::vector<double> params, params_std;
	double Rsq;

	std::vector<DBL2> stencil_xy_data_adiabatic_y(n_stencil.dim(), DBL2()), stencil_xy_data_nonadiabatic_y(n_stencil.dim(), DBL2());
	std::vector<double> stencil_adiabatic_x(n_stencil.dim(), 0.0), stencil_nonadiabatic_x(n_stencil.dim(), 0.0);

	//setup indexes in xy_data to fit : these are fixed for the stencil
	for (int idx = 0; idx < n_stencil.dim(); idx++) {

		stencil_xy_data_adiabatic_y[idx].x = idx;
		stencil_xy_data_nonadiabatic_y[idx].x = idx;
	}

	//now fit point by point
	for (int i = (n_stencil.x - 1) / 2; i < n.i - (n_stencil.x - 1) / 2; i++) {
		for (int j = (n_stencil.y - 1) / 2; j < n.j - (n_stencil.y - 1) / 2; j++) {
			for (int k = (n_stencil.z - 1) / 2; k < n.k - (n_stencil.z - 1) / 2; k++) {

				if (M.is_not_empty(INT3(i, j, k)) && ((T[INT3(i, j, k)].norm() / T_max > Torque_ratio_threshold) || IsZ(Torque_ratio_threshold))) {

					//setup stencil data
					for (int i_s = -(n_stencil.x - 1) / 2; i_s <= (n_stencil.x - 1) / 2; i_s++) {
						for (int j_s = -(n_stencil.y - 1) / 2; j_s <= (n_stencil.y - 1) / 2; j_s++) {

							//index to read from main arrays
							int idx = (i + i_s) + (j + j_s) * n.x + k * n.x * n.y;

							//index in stencil arrays
							int idx_stencil = (i_s + (n_stencil.x - 1) / 2) + (j_s + (n_stencil.y - 1) / 2) * n_stencil.x;

							stencil_xy_data_adiabatic_y[idx_stencil].y = xy_data_adiabatic_y[idx];
							stencil_xy_data_nonadiabatic_y[idx_stencil].y = xy_data_nonadiabatic_y[idx];
							stencil_adiabatic_x[idx_stencil] = adiabatic_x[idx];
							stencil_nonadiabatic_x[idx_stencil] = nonadiabatic_x[idx];
						}
					}

					//fit stencil
					curve_fit.FitSTT_LMA(stencil_xy_data_adiabatic_y, stencil_xy_data_nonadiabatic_y, stencil_adiabatic_x, stencil_nonadiabatic_x, params, params_std, &Rsq);

					//set output if threshold criteria are met only
					if (
						((params[1] && abs(params_std[1] / params[1]) < absolute_error_threshold) || IsZ(absolute_error_threshold)) &&
						((Rsq > Rsq_threshold) || IsZ(Rsq_threshold)))
					{
						if (adiabatic) output[INT3(i, j, k)] = params[0];
						else output[INT3(i, j, k)] = params[1];
					}
				}
			}
		}
	}

	return error;
}

//fit SOT function :
BError DPArrays::fit_sot(VEC<DBL3>& T, VEC_VC<DBL3>& M, VEC_VC<DBL3>& J, Rect rectangle, DBL2* pSHAeff, DBL2 *pflST, double *pRsq)
{
	BError error(__FUNCTION__);

	//change the rectangle to absolute coordinates
	rectangle += M.rect.s;
	
	//get contacting rectangle between J and M, as well as intersection with the user specified rectangle
	Rect contact_rect = M.rect.get_intersection(J.rect).get_intersection(rectangle);

	if (contact_rect.IsNull()) {

		//ERROR: heavy metal mesh must be in contact with the focused ferromagnetic mesh and any specified rectangle must intersect the contact.
		return error(BERROR_INCORRECTCONFIG);
	}

	if (!contact_rect.IsPlane() || IsNZ(contact_rect.height())) {

		//ERROR: contacting meshes must be stacked in the z direction; special case where M and J have same rectangle (same mesh) is valid
		if (M.rect != J.rect) return error(BERROR_INCORRECTCONFIG);
	}
	
	double dh = M.h.z;

	Box cells_box = M.box_from_rect_max(contact_rect);
	INT3 box_sizes = cells_box.size();
	int points = box_sizes.dim();

	std::vector<DBL2> xy_data_dampinglike_y(points, DBL2()), xy_data_fieldlike_y(points, DBL2());
	std::vector<double> dampinglike_x(points, 0.0), fieldlike_x(points, 0.0);

#pragma omp parallel for
	for (int i = cells_box.s.i; i < cells_box.e.i; i++) {
		for (int j = cells_box.s.j; j < cells_box.e.j; j++) {
			for (int k = cells_box.s.k; k < cells_box.e.k; k++) {

				DBL3 position = M.cellidx_to_position(INT3(i, j, k));
				DBL3 position_J;

				if (M.rect != J.rect) {

					if (M.rect.s.z > J.rect.s.z) {

						position_J = position + M.rect.s - J.rect.s - DBL3(0, 0, M.h.z / 2) - DBL3(0, 0, J.h.z / 2);
					}
					else {

						position_J = position + M.rect.s - J.rect.s + DBL3(0, 0, M.h.z / 2) + DBL3(0, 0, J.h.z / 2);
					}
				}
				else {

					position_J = position;
				}

				int idx_M = M.position_to_cellidx(position);
				int idx_J = J.position_to_cellidx(position_J);
				int idx = (i - cells_box.s.i) + (j - cells_box.s.j) * box_sizes.x + (k - cells_box.s.k) * box_sizes.x * box_sizes.y;

				if (M.is_not_empty(idx_M) && J.is_not_empty(idx_J)) {

					DBL3 m = normalize(M[idx_M]);
					DBL3 p = (DBL3(0, 0, 1) ^ J[idx_J]) / dh;

					DBL3 mxp = m ^ p;
					DBL3 mxmxp = m ^ mxp;

					xy_data_dampinglike_y[idx] = DBL2(idx, mxmxp * T[idx_M]);
					xy_data_fieldlike_y[idx] = DBL2(idx, mxp * T[idx_M]);

					dampinglike_x[idx] = mxmxp * mxmxp;
					fieldlike_x[idx] = mxp * mxp;
				}
			}
		}
	}

	CurveFitting curve_fit;

	std::vector<double> params, params_std;

	curve_fit.FitSOT_LMA(xy_data_dampinglike_y, xy_data_fieldlike_y, dampinglike_x, fieldlike_x, params, params_std, pRsq);

	if (params.size() >= 2 && params_std.size() >= 2) {

		*pSHAeff = DBL2(params[0], params_std[0]);
		*pflST = DBL2(params[1], params_std[1]);
	}
	else return error(BERROR_INCORRECTARRAYS);

	return error;
}

//fit simultaneously STT and SOT when the self-consistent interfacial spin torque contains both
BError DPArrays::fit_sotstt(VEC<DBL3>& T, VEC_VC<DBL3>& M, VEC_VC<DBL3>& J_hm, VEC_VC<DBL3>& J_fm, Rect rectangle, DBL2* pSHAeff, DBL2 *pflST, DBL2* pP, DBL2 *pbeta, double *pRsq)
{
	BError error(__FUNCTION__);

	//change the rectangle to absolute coordinates
	rectangle += M.rect.s;

	//get contacting rectangle between J and M, as well as intersection with the user specified rectangle
	Rect contact_rect = M.rect.get_intersection(J_hm.rect).get_intersection(rectangle);

	if (contact_rect.IsNull()) {

		//ERROR: heavy metal mesh must be in contact with the focused ferromagnetic mesh and any specified rectangle must intersect the contact.
		return error(BERROR_INCORRECTCONFIG);
	}

	if (!contact_rect.IsPlane() || IsNZ(contact_rect.height())) {

		//ERROR: contacting meshes must be stacked in the z direction; special case where M and J have same rectangle (same mesh) is valid
		if (M.rect != J_hm.rect) return error(BERROR_INCORRECTCONFIG);
	}

	double dh = M.h.z;

	Box cells_box = M.box_from_rect_max(contact_rect);
	INT3 box_sizes = cells_box.size();
	int points = box_sizes.dim();

	std::vector<DBL2> xy_data_adiabatic_y(points, DBL2()), xy_data_nonadiabatic_y(points, DBL2());
	std::vector<DBL2> xy_data_dampinglike_y(points, DBL2()), xy_data_fieldlike_y(points, DBL2());

	std::vector<double> adiabatic_fit_main_x(points, 0.0), adiabatic_fit_dampinglike_x(points, 0.0), adiabatic_fit_fieldlike_x(points, 0.0);
	std::vector<double> nonadiabatic_fit_main_x(points, 0.0), nonadiabatic_fit_dampinglike_x(points, 0.0), nonadiabatic_fit_fieldlike_x(points, 0.0);
	std::vector<double> dampinglike_fit_main_x(points, 0.0), dampinglike_fit_adiabatic_x(points, 0.0), dampinglike_fit_nonadiabatic_x(points, 0.0);
	std::vector<double> fieldlike_fit_main_x(points, 0.0), fieldlike_fit_adiabatic_x(points, 0.0), fieldlike_fit_nonadiabatic_x(points, 0.0);

#pragma omp parallel for
	for (int i = cells_box.s.i; i < cells_box.e.i; i++) {
		for (int j = cells_box.s.j; j < cells_box.e.j; j++) {
			for (int k = cells_box.s.k; k < cells_box.e.k; k++) {

				DBL3 position = M.cellidx_to_position(INT3(i, j, k));
				DBL3 position_J_hm;
				
				if (M.rect != J_hm.rect) {

					if (M.rect.s.z > J_hm.rect.s.z) {

						position_J_hm = position + M.rect.s - J_hm.rect.s - DBL3(0, 0, M.h.z / 2) - DBL3(0, 0, J_hm.h.z / 2);
					}
					else {

						position_J_hm = position + M.rect.s - J_hm.rect.s + DBL3(0, 0, M.h.z / 2) + DBL3(0, 0, J_hm.h.z / 2);
					}
				}
				else {

					position_J_hm = position;
				}

				int idx_M = M.position_to_cellidx(position);
				int idx_J_hm = J_hm.position_to_cellidx(position_J_hm);
				int idx_J_fm = J_fm.position_to_cellidx(position);
				int idx = (i - cells_box.s.i) + (j - cells_box.s.j) * box_sizes.x + (k - cells_box.s.k) * box_sizes.x * box_sizes.y;

				if (M.is_not_empty(idx_M) && J_hm.is_not_empty(idx_J_hm) && J_fm.is_not_empty(idx_J_fm)) {

					DBL3 m = normalize(M[idx_M]);
					DBL3 p = (DBL3(0, 0, 1) ^ J_hm[idx_J_hm]) / dh;

					DBL3 mxp = m ^ p;
					DBL3 mxmxp = m ^ mxp;

					DBL33 grad_m = normalize(M.grad_neu(idx_M), M[idx_M]);
					DBL3 J_dot_del_m = grad_m.x * J_fm[idx_J_fm].x + grad_m.y * J_fm[idx_J_fm].y + grad_m.z * J_fm[idx_J_fm].z;
					DBL3 m_cross_J_dot_del_m = m ^ J_dot_del_m;

					xy_data_dampinglike_y[idx] = DBL2(idx, mxmxp * T[idx_M]);
					xy_data_fieldlike_y[idx] = DBL2(idx, mxp * T[idx_M]);
					xy_data_adiabatic_y[idx] = DBL2(idx, J_dot_del_m * T[idx_M]);
					xy_data_nonadiabatic_y[idx] = DBL2(idx, m_cross_J_dot_del_m * T[idx_M]);

					dampinglike_fit_main_x[idx] = mxmxp * mxmxp;
					dampinglike_fit_adiabatic_x[idx] = mxmxp * J_dot_del_m;
					dampinglike_fit_nonadiabatic_x[idx] = mxmxp * m_cross_J_dot_del_m;

					fieldlike_fit_main_x[idx] = mxp * mxp;
					fieldlike_fit_adiabatic_x[idx] = mxp * J_dot_del_m;
					fieldlike_fit_nonadiabatic_x[idx] = mxp * m_cross_J_dot_del_m;

					adiabatic_fit_main_x[idx] = J_dot_del_m * J_dot_del_m;
					adiabatic_fit_dampinglike_x[idx] = J_dot_del_m * mxmxp;
					adiabatic_fit_fieldlike_x[idx] = J_dot_del_m * mxp;

					nonadiabatic_fit_main_x[idx] = m_cross_J_dot_del_m * m_cross_J_dot_del_m;
					nonadiabatic_fit_dampinglike_x[idx] = m_cross_J_dot_del_m * mxmxp;
					nonadiabatic_fit_fieldlike_x[idx] = m_cross_J_dot_del_m * mxp;
				}
			}
		}
	}

	CurveFitting curve_fit;

	std::vector<double> params, params_std;
	
	curve_fit.FitSOTSTT_LMA(xy_data_adiabatic_y, xy_data_nonadiabatic_y, xy_data_dampinglike_y, xy_data_fieldlike_y,
		adiabatic_fit_main_x, adiabatic_fit_dampinglike_x, adiabatic_fit_fieldlike_x,
		nonadiabatic_fit_main_x, nonadiabatic_fit_dampinglike_x, nonadiabatic_fit_fieldlike_x,
		dampinglike_fit_main_x, dampinglike_fit_adiabatic_x, dampinglike_fit_nonadiabatic_x,
		fieldlike_fit_main_x, fieldlike_fit_adiabatic_x, fieldlike_fit_nonadiabatic_x, 
		params, params_std, pRsq);

	if (params.size() >= 4 && params_std.size() >= 4) {

		*pP = DBL2(params[0], params_std[0]);
		*pbeta = DBL2(params[1], params_std[1]);
		*pSHAeff = DBL2(params[2], params_std[2]);
		*pflST = DBL2(params[3], params_std[3]);
	}
	else return error(BERROR_INCORRECTARRAYS);
	
	return error;
}

//--------------------- data processing

//Replace repeated points from using linear interpolation: if two adjacent sets of repeated points found, replace repeats between the mid-points of the sets
BError DPArrays::replace_repeats(int dp_in, int dp_out)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_in, dp_out) || !dpA[dp_in].size()) return error(BERROR_INCORRECTARRAYS);

	if (dp_in != dp_out) resize(dp_out, dpA[dp_in].size());

	//first remove all "blips" : a blip is a single point which has at least 2 points either side with the same value. 
	//Exceptions: 
	//The point at the end is not considered a blip.
	//The point at the start is considered a blip if it has at least 2 points to the right with the same value.
	//The points next-to-last and next to start are blips if the points either side have the same value.
	//
	//Blips are removed by setting them to the value of the surrounding points.

	//First 2 points treated differently
	if (dpA[dp_in].size() >= 3) {

		//for idx = 0
		double value = dpA[dp_in][1];
		
		if (dpA[dp_in][2] == value && dpA[dp_in][0] != value) {

			dpA[dp_out][0] = value;
		}
		else dpA[dp_out][0] = dpA[dp_in][0];

		//for idx = 1
		value = dpA[dp_in][2];

		if (dpA[dp_in][0] == value && dpA[dp_in][1] != value) {

			dpA[dp_out][1] = value;
		}
		else dpA[dp_out][1] = dpA[dp_in][1];
	}

	//The middle points
	for (int idx = 2; idx < dpA[dp_in].size() - 2; idx++) {

		double value = dpA[dp_in][idx - 2];

		if (dpA[dp_in][idx - 1] == value && dpA[dp_in][idx + 1] == value && dpA[dp_in][idx + 2] == value && dpA[dp_in][idx] != value) {

			//found blip
			dpA[dp_out][idx] = value;
		}
		else dpA[dp_out][idx] = dpA[dp_in][idx];
	}

	//Next-to-last point : idx = dpA[dp_in].size() - 2
	if (dpA[dp_in].size() >= 3) {

		int index = dpA[dp_in].size() - 2;

		double value = dpA[dp_in][index + 1];

		if (dpA[dp_in][index - 1] == value && dpA[dp_in][index] != value) {

			dpA[dp_out][index] = value;
		}
		else dpA[dp_out][index] = dpA[dp_in][index];
	}

	//last point
	dpA[dp_out].back() = dpA[dp_in].back();

	//Now that blips have been removed find the first point that differs from a sequence of points with same value, then replace all points with same value using linear interpolation

	double value_start = dpA[dp_out][0];
	int index_start = 0;
	int index_previous = 0;

	bool first = true;
	int first_index_start = 0;
	int first_index_end = 0;
	
	bool got_gradient = false;

	for (int idx = 1; idx < dpA[dp_out].size(); idx++) {

		if (dpA[dp_out][idx] != value_start) {

			//found a point that differs : replace with interpolation

			if (!first) {

				if (!got_gradient) {

					got_gradient = true;

					double gradient = (dpA[dp_out][idx] - value_start) / (idx - index_start);

					//replace first set using the gradient from this second sequence
					for (int i_idx = first_index_end - 1; i_idx >= first_index_start; i_idx--) {

						dpA[dp_out][i_idx] = (i_idx - first_index_end) * gradient + dpA[dp_out][first_index_end];
					}
				}

				//interpolation if not first time

				for (int i_idx = index_start; i_idx < idx; i_idx++) {

					dpA[dp_out][i_idx] = (dpA[dp_out][idx] * (i_idx - index_start) + value_start * (idx - i_idx)) / (idx - index_start);
				}
			}
			else {

				//first time : save for later. We'll replace the first set with interpolation after we get the gradient from the next sequence

				first = false;

				first_index_start = index_start;
				first_index_end = idx;
			}

			//prepare search for next one
			index_previous = index_start;
			index_start = idx;
			value_start = dpA[dp_out][idx];
		}
	}

	if (!first && got_gradient) {

		//reached the end, but was there a sequence of constant value points at the end we didn't replace? Should we do anything about it?
		int last_sequence_length = dpA[dp_out].size() - 1 - index_start;
		int prev_sequence_length = index_start - index_previous;

		//if the last sequence is at least 1.5 times longer than the previous one then do nothing - 
		//in this case it's considered likely the constant values are real and would not have incremented if the data sequence was captured for longer
		//Otherwise replace with previous gradient

		if (prev_sequence_length > 0 && (double)last_sequence_length / prev_sequence_length < 1.5) {

			double gradient = (dpA[dp_out][index_start] - dpA[dp_out][index_previous]) / prev_sequence_length;

			for (int i_idx = index_start; i_idx < dpA[dp_out].size(); i_idx++) {

				dpA[dp_out][i_idx] = dpA[dp_out][index_start] + gradient * (i_idx - index_start);
			}
		}
	}
	else {

		//just replace start to finish with linear gradient
		//replace first set using the gradient from this second sequence
		for (int i_idx = 0; i_idx < dpA[dp_out].size(); i_idx++) {

			dpA[dp_out][i_idx] = (dpA[dp_out].back() * i_idx + dpA[dp_out].front() * (dpA[dp_out].size() - 1 - i_idx)) / (dpA[dp_out].size() - 1);
		}
	}

	return error;
}

//Subtract the first point (the offset) from all the points in dp_in
BError DPArrays::remove_offset(int dp_in, int dp_out)
{
	BError error(__FUNCTION__);

	if (!GoodArrays(dp_in, dp_out) || !dpA[dp_in].size()) return error(BERROR_INCORRECTARRAYS);

	if (dp_in != dp_out) resize(dp_out, dpA[dp_in].size());

	double offset = dpA[dp_in][0];

#pragma omp parallel for
	for (int idx = 0; idx < dpA[dp_in].size(); idx++) {

		dpA[dp_out][idx] = dpA[dp_in][idx] - offset;
	}

	return error;
}

//smooth data in dp_in using nearest-neighbor averaging with given window size, and place result in dp_out (must be different)
BError DPArrays::adjacent_averaging(int dp_in, int dp_out, int window_size)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_in, dp_out) || !dpA[dp_in].size() || window_size / 2 > dpA[dp_in].size() || window_size < 2 || !dpA[dp_in].size()) return error(BERROR_INCORRECTARRAYS);

	resize(dp_out, dpA[dp_in].size());

	//build starting average
	double average = 0.0;

	for (int idx = 0; idx < window_size / 2; idx++) {

		average += dpA[dp_in][idx];
	}

	average /= (window_size / 2);
	dpA[dp_out][0] = average;

	//the starting tail with partial window
	for (int idx = 1; idx < window_size / 2; idx++) {

		average = (average * (window_size / 2 + idx - 1) + dpA[dp_in][idx + window_size / 2 - 1]) / (window_size / 2 + idx);

		dpA[dp_out][idx] = average;
	}

	//the middle part with full window
	for (int idx = window_size / 2; idx < dpA[dp_in].size() - window_size / 2; idx++) {

		average = average - (dpA[dp_in][idx - window_size / 2] - dpA[dp_in][idx + window_size / 2 - 1]) / window_size;

		dpA[dp_out][idx] = average;
	}

	//the ending tail with partial window
	for (int idx = dpA[dp_in].size() - window_size / 2; idx < dpA[dp_in].size(); idx++) {

		average = (average * (window_size / 2 + dpA[dp_in].size() - idx) - dpA[dp_in][idx]) / (window_size / 2 + dpA[dp_in].size() - idx - 1);

		dpA[dp_out][idx] = average;
	}

	return error;
}

//extract monotonic sequence form input arrays and place it in output arrays. Primary array (dp_in_pri) decides ordering and dependent array (dp_in_dep) must have same size
BError DPArrays::extract_monotonic(int dp_in_pri, int dp_in_dep, int dp_out_pri, int dp_out_dep)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_in_pri, dp_in_dep, dp_out_pri, dp_out_dep) || !dpA[dp_in_pri].size() || dpA[dp_in_pri].size() != dpA[dp_in_dep].size()) return error(BERROR_INCORRECTARRAYS);

	resize(dp_out_pri, 0);
	resize(dp_out_dep, 0);

	bool increasing = dpA[dp_in_pri].back() > dpA[dp_in_pri].front();

	double last_value = dpA[dp_in_pri][0];
	dpA[dp_out_pri].push_back(last_value);
	dpA[dp_out_dep].push_back(dpA[dp_in_dep][0]);

	for (int idx = 1; idx < dpA[dp_in_pri].size(); idx++) {

		if (dpA[dp_in_pri][idx] >= last_value) {

			last_value = dpA[dp_in_pri][idx];
			dpA[dp_out_pri].push_back(last_value);
			dpA[dp_out_dep].push_back(dpA[dp_in_dep][idx]);
		}
	}

	return error;
}

//convert from Cartesian to polar (as r, theta)
BError DPArrays::Cartesian_to_Polar(int dp_in_x, int dp_in_y, int dp_out_r, int dp_out_theta)
{
	BError error(__FUNCTION__);

	if (!GoodArrays_Unique(dp_in_x, dp_in_y, dp_out_r, dp_out_theta)) return error(BERROR_INCORRECTARRAYS);

	resize(dp_out_r, dpA[dp_in_x].size());
	resize(dp_out_theta, dpA[dp_in_y].size());

#pragma omp parallel for
	for (int idx = 0; idx < (dpA[dp_in_x].size() < dpA[dp_in_y].size() ? dpA[dp_in_x].size() : dpA[dp_in_y].size()); idx++) {

		DBL2 polar = ::Cartesian_to_Polar<DBL2>(DBL2(dpA[dp_in_x][idx], dpA[dp_in_y][idx]));

		dpA[dp_out_r][idx] = polar.x;
		dpA[dp_out_theta][idx] = polar.y;
	}

	return error;
}