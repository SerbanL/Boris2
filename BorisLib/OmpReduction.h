#pragma once

#include <omp.h>
#include <vector>

#include "Types.h"

//Implements reductions for minimum, maximum and running average within OpernMP parallel for loops

//Under Visual Studio (currently working with VS2017) only OMP 2.0 is supported !!! (crazy isn't it?). This means min and max reduction needs to be done manually.
//Use the OmpReduction class to make this process easier.

//Example usage:

//Before a parallel loop (e.g. reducing double values) declare :

//OmpReduction<double> omp_reduction;

//You could also have this as a data member in the class where the parallel loop is used. In this case before the loop call:

//omp_reduction.new_minmax_reduction();

//Inside the loop use this on each value you want to include in the reduction (e.g. reducing both minimum and maximum):
//omp_reduction.reduce_minmax(loop_value);

//after the loop obtain the final result as:

//omp_reduction.minmax();


template <typename Type>
class OmpReduction {

	//number of omp threads
	int OmpThreads;
	
	//vectors of size OmpThreads for minimum and maximum reductions
	std::vector<Type> minimum_reduction;
	std::vector<Type> maximum_reduction;
	
	//int not bool : don't use vector<bool> as it's not a real bool array and you'll have some very unpleasant bugs (e.g. #pragma omp parallel for doesn't work with it!!!!)
	std::vector<int> reduction_value_not_set;

	//vectors of size OmpThreads for average reduction : running average on each thread then combine at the end : need to keep track of running averages as well as number of points in averages
	std::vector<Type> average_reduction;
	std::vector<int> average_reduction_points;

public:

	//values available after reduction
	Type min, max, av;

public:

	OmpReduction(void)
	{
		OmpThreads = omp_get_num_procs();

		minimum_reduction.assign(OmpThreads, Type());
		maximum_reduction.assign(OmpThreads, Type());

		reduction_value_not_set.assign(OmpThreads, true);

		average_reduction.assign(OmpThreads, Type());
		average_reduction_points.assign(OmpThreads, 0);

		min = Type();
		max = Type();
		av = Type();
	}
	
	//--------------------- MINIMUM and MAXIMUM

	//call this before starting a new reduction - place before the #pragma omp parallel for
	void new_minmax_reduction(void)
	{
#pragma omp parallel for
		for (int idx = 0; idx < OmpThreads; idx++) {

			reduction_value_not_set[idx] = true;
		}
	}

	//within the parallel loop reduce values by calling this : minimum only
	void reduce_min(Type value)
	{
		int tn = omp_get_thread_num();

		if (reduction_value_not_set[tn]) {

			reduction_value_not_set[tn] = false;
			minimum_reduction[tn] = value;
		}
		else {

			if (value < minimum_reduction[tn]) minimum_reduction[tn] = value;
		}
	}

	//within the parallel loop reduce values by calling this : maximum only
	void reduce_max(Type value)
	{
		int tn = omp_get_thread_num();

		if (reduction_value_not_set[tn]) {

			reduction_value_not_set[tn] = false;
			maximum_reduction[tn] = value;
		}
		else {

			if (value > maximum_reduction[tn]) maximum_reduction[tn] = value;
		}
	}

	//within the parallel loop reduce values by calling this : both minimum and maximum
	void reduce_minmax(Type value)
	{
		int tn = omp_get_thread_num();

		if (reduction_value_not_set[tn]) {

			reduction_value_not_set[tn] = false;
			minimum_reduction[tn] = value;
			maximum_reduction[tn] = value;
		}
		else {

			if (value > maximum_reduction[tn]) maximum_reduction[tn] = value;
			if (value < minimum_reduction[tn]) minimum_reduction[tn] = value;
		}
	}

	//after reduction call this to find reduced value : minimum only
	Type minimum(void)
	{
		int idx;

		for (idx = 0; idx < OmpThreads; idx++) {

			//it's possible not all threads have collected a reduced value
			if (reduction_value_not_set[idx]) continue;

			//found a starting value
			min = minimum_reduction[idx];

			break;
		}

		//reduce remaining threads
		for (++idx; idx < OmpThreads; idx++) {

			if (!reduction_value_not_set[idx]) {

				if (minimum_reduction[idx] < min) min = minimum_reduction[idx];
			}
		}

		return min;
	}

	//after reduction call this to find reduced value : maximum only
	Type maximum(void)
	{
		int idx;

		for (idx = 0; idx < OmpThreads; idx++) {

			//it's possible not all threads have collected a reduced value
			if (reduction_value_not_set[idx]) continue;

			//found a starting value
			max = maximum_reduction[idx];

			break;
		}

		//reduce remaining threads
		for (++idx; idx < OmpThreads; idx++) {

			if (!reduction_value_not_set[idx]) {

				if (maximum_reduction[idx] > max) max = maximum_reduction[idx];
			}
		}

		return max;
	}

	//after reduction call this to find reduced value : both minimum and maximum stored in a VAL2 in this order
	VAL2<Type> minmax(void)
	{
		int idx;

		for (idx = 0; idx < OmpThreads; idx++) {

			//it's possible not all threads have collected a reduced value
			if (reduction_value_not_set[idx]) continue;

			//found a starting value
			min = minimum_reduction[idx];
			max = maximum_reduction[idx];

			break;
		}

		//reduce remaining threads
		for (++idx; idx < OmpThreads; idx++) {

			if (!reduction_value_not_set[idx]) {

				if (minimum_reduction[idx] < min) min = minimum_reduction[idx];
				if (maximum_reduction[idx] > max) max = maximum_reduction[idx];
			}
		}

		return VAL2<Type>(min, max);
	}

	//--------------------- AVERAGE

	//call this before starting a new average reduction - place before the #pragma omp parallel for
	void new_average_reduction(void)
	{
#pragma omp parallel for
		for (int idx = 0; idx < OmpThreads; idx++) {

			average_reduction[idx] = Type();
			average_reduction_points[idx] = 0;
		}
	}

	//reduce for average during a loop
	void reduce_average(Type value)
	{
		int tn = omp_get_thread_num();

		//get running average on this thread
		average_reduction_points[tn]++;
		average_reduction[tn] = (average_reduction[tn] * (average_reduction_points[tn] - 1) + value) / average_reduction_points[tn];
	}

	//call this after a loop to get the final average
	Type average(void)
	{
		av = Type();

		//count total points
		int points_count = 0;
		for (int idx = 0; idx < OmpThreads; idx++) {

			points_count += average_reduction_points[idx];
		}

		if (points_count) {

			//now combine partial averages into a final average
			for (int idx = 0; idx < OmpThreads; idx++) {

				av += average_reduction[idx] * ((double)average_reduction_points[idx] / points_count);
			}
		}

		return av;
	}
};
