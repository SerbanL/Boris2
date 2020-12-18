#pragma once

#include "Types_VAL.h"
#include "Funcs_Math.h"
#include "Funcs_Vectors.h"
#include "Funcs_Files.h"

////////////////////////////////////////////////////////////////////////////////////////////////// Sequence
//
// Has start and end points, with an integer number of steps between them. Minimum of one step (in this case the sequence has just the start and end points)

template <typename Type>
class Sequence {

	//----------------------------- DATA

private:

	Type start, end;
	int steps;

public:

	//----------------------------- VALUE CONSTRUCTORS

	template <typename Type_ = Type>
	Sequence(void) {

		start = Type();
		end = Type();
		steps = 1;
	}

	template <typename Type_ = Type>
	Sequence(Type_ start_, Type_ end_, int steps_) :
		start(start_), end(end_), steps(steps_)
	{}

	//For the case of SEQ3 (i.e. Sequence<DBL3>) we also want to have a 7 parameter constructor as below, but this will not be available for any other types
	template <typename Type_ = Type, std::enable_if_t<std::is_same<Type_, DBL3>::value>* = nullptr>
	Sequence(double sx, double sy, double sz, double ex, double ey, double ez, int steps_) :
		start(DBL3(sx, sy, sz)), end(DBL3(ex, ey, ez)), steps(steps_)
	{}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	Sequence(const Sequence &copyThis) { start = copyThis.start; end = copyThis.end; steps = copyThis.steps; }

	//assignment operator
	Sequence& operator=(const Sequence &rhs) { start = rhs.start; end = rhs.end; steps = rhs.steps; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const Sequence &rhs) { os << ToString(rhs.start) << "; " << ToString(rhs.end) << "; " << rhs.steps; return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const Sequence &rhs) {

		lhs << rhs.start << std::string("; ") << rhs.end << std::string("; ") << ToString(rhs.steps);
		//steps cannot have units, it's just an integer, so pre-convert it to std::string as above

		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend Sequence& operator>>(const std::stringstream &ss, Sequence &rhs)
	{
		//normally the values in std::string representation are as: "start; end; steps"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "start end steps"
		if (components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings
		switch (components.size()) {

		case 3:			//Type is e.g. a double
			rhs.start = ToNum(trimspaces(components[0]), "");
			rhs.end = ToNum(trimspaces(components[1]), "");
			rhs.steps = ToNum(trimspaces(components[2]), "");
			break;
		case 7:			//Type is a 3-component value
			rhs.start = ToNum(trimspaces(components[0]) + ", " + trimspaces(components[1]) + ", " + trimspaces(components[2]), "");
			rhs.end = ToNum(trimspaces(components[3]) + ", " + trimspaces(components[4]) + ", " + trimspaces(components[5]), "");
			rhs.steps = ToNum(trimspaces(components[6]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	Type value(int stepIdx)
	{
		Type increment = (end - start) / steps;

		return (start + increment * stepIdx);
	}

	int number_of_steps(void) { return steps; }
};

typedef Sequence<double> SEQ;
typedef Sequence<DBL3> SEQ3;

////////////////////////////////////////////////////////////////////////////////////////////////// Sequence Polar
//
// Similar to a SEQ3 but start and end values specified in spherical polar coordinates as : magnitude, theta (polar), phi (azimuthal). The output values are still in Cartesian coordinates

template <typename Type = void>			//SEQP requires a number of methods defined in other headers which cannot be included here. SEQP is only used in files which do include those methods, so fake template this to disable compilation here.
class SequencePolar {

	//----------------------------- DATA

private:

	DBL3 start = DBL3();
	DBL3 end = DBL3();
	int steps = 1;

public:

	//----------------------------- VALUE CONSTRUCTORS

	SequencePolar(void) {}

	SequencePolar(DBL3 start_, DBL3 end_, int steps_) :
		start(start_), end(end_), steps(steps_)
	{}

	SequencePolar(double mag_s, double theta_s, double phi_s, double mag_se, double theta_e, double phi_e, int steps_) :
		start(DBL3(mag_s, theta_s, phi_s)), end(DBL3(mag_se, theta_e, phi_e)), steps(steps_)
	{}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	SequencePolar(const SequencePolar &copyThis) { start = copyThis.start; end = copyThis.end; steps = copyThis.steps; }

	//assignment operator
	SequencePolar& operator=(const SequencePolar &rhs) { start = rhs.start; end = rhs.end; steps = rhs.steps; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const SequencePolar &rhs) { os << ToString(rhs.start) << "; " << ToString(rhs.end) << "; " << rhs.steps; return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const SequencePolar &rhs)
	{
		lhs << rhs.start.i << std::string(", ") << ToString(DBL2(rhs.start.j, rhs.start.k)) << std::string("; ");
		lhs << rhs.end.i << std::string(", ") << ToString(DBL2(rhs.end.j, rhs.end.k)) << std::string("; ");
		lhs << ToString(rhs.steps);

		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend SequencePolar& operator>>(const std::stringstream &ss, SequencePolar &rhs)
	{
		//normally the values in std::string representation are as: "start; end; steps"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "start end steps"
		if (components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings
		switch (components.size()) {

		case 3:			//Type is e.g. a double
			rhs.start = ToNum(trimspaces(components[0]), "");
			rhs.end = ToNum(trimspaces(components[1]), "");
			rhs.steps = ToNum(trimspaces(components[2]), "");
			break;
		case 7:			//Type is a 3-component value
			rhs.start = ToNum(trimspaces(components[0]) + ", " + trimspaces(components[1]) + ", " + trimspaces(components[2]), "");
			rhs.end = ToNum(trimspaces(components[3]) + ", " + trimspaces(components[4]) + ", " + trimspaces(components[5]), "");
			rhs.steps = ToNum(trimspaces(components[6]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	DBL3 value(int stepIdx)
	{
		DBL3 increment = (end - start) / steps;

		return Polar_to_Cartesian(start + increment * stepIdx);
	}

	int number_of_steps(void) { return steps; }
};

typedef SequencePolar<void> SEQP;

////////////////////////////////////////////////////////////////////////////////////////////////// Cos function oscillation.
//
// 

template <typename Type>
class CosOscillation {

	//----------------------------- DATA

private:

	//the oscillation value : specify the starting point of the cos function; the cos sequence will oscillate between +oscillation and -oscillation
	Type oscillation;
	//to adjust the frequency the time taken per step must be defined externally, this is of no concern here (e.g. 1 GHz with 20 steps per cycle : increment steps every 50 ps)
	int steps_per_cycle, cycles;

public:

	//----------------------------- VALUE CONSTRUCTORS

	template <typename Type_ = Type>
	CosOscillation(void) {

		oscillation = Type();
		steps_per_cycle = 1;
		cycles = 1;
	}

	template <typename Type_ = Type>
	CosOscillation(Type_ oscillation_, int steps_per_cycle_, int cycles_) :
		oscillation(oscillation_), steps_per_cycle(steps_per_cycle_), cycles(cycles_)
	{}

	//For the case of COSOSC3 (i.e. CosOscillation<DBL3>) we also want to have a 5 parameter constructor as below, but this will not be available for any other types
	template <typename Type_ = Type, std::enable_if_t<std::is_same<Type_, DBL3>::value>* = nullptr>
	CosOscillation(double ox, double oy, double oz, int steps_per_cycle_, int cycles_) :
		oscillation(DBL3(ox, oy, oz)), steps_per_cycle(steps_per_cycle_), cycles(cycles_)
	{}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	CosOscillation(const CosOscillation &copyThis) { oscillation = copyThis.oscillation; steps_per_cycle = copyThis.steps_per_cycle; cycles = copyThis.cycles; }

	//assignment operator
	CosOscillation& operator=(const CosOscillation &rhs) { oscillation = rhs.oscillation; steps_per_cycle = rhs.steps_per_cycle; cycles = rhs.cycles; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const CosOscillation &rhs) { os << ToString(rhs.oscillation) << "; " << rhs.steps_per_cycle << "; " << rhs.cycles; return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const CosOscillation &rhs) {

		lhs << rhs.oscillation << std::string("; ") << ToString(rhs.steps_per_cycle) << std::string("; ") << ToString(rhs.cycles);

		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend CosOscillation& operator>>(const std::stringstream &ss, CosOscillation &rhs)
	{
		//normally the values in std::string representation are as: "oscillation; steps_per_cycle; cycles"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "oscillation steps_per_cycle cycles"
		if (components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings
		switch (components.size()) {

		case 3:			//Type is e.g. a double
			rhs.oscillation = ToNum(trimspaces(components[0]), "");
			rhs.steps_per_cycle = ToNum(trimspaces(components[1]), "");
			rhs.cycles = ToNum(trimspaces(components[2]), "");
			break;
		case 5:			//Type is e.g. a DBL3
			rhs.oscillation = ToNum(trimspaces(components[0]) + ", " + trimspaces(components[1]) + ", " + trimspaces(components[2]), "");
			rhs.steps_per_cycle = ToNum(trimspaces(components[3]), "");
			rhs.cycles = ToNum(trimspaces(components[4]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	Type value(int stepIdx)
	{
		return oscillation * cos(6.2831853071795864769252867666 * (double)stepIdx / steps_per_cycle);
	}

	int number_of_steps(void) { return steps_per_cycle * cycles; }
};

typedef CosOscillation<double> COSOSC;
typedef CosOscillation<DBL3> COSOSC3;

////////////////////////////////////////////////////////////////////////////////////////////////// Sin function oscillation.
//
// 

template <typename Type>
class SinOscillation {

	//----------------------------- DATA

private:

	//the oscillation value : specify the starting point of the cos function; the cos sequence will oscillate between +oscillation and -oscillation
	Type oscillation;
	//to adjust the frequency the time taken per step must be defined externally, this is of no concern here (e.g. 1 GHz with 20 steps per cycle : increment steps every 50 ps)
	int steps_per_cycle, cycles;

public:

	//----------------------------- VALUE CONSTRUCTORS

	template <typename Type_ = Type>
	SinOscillation(void) {

		oscillation = Type();
		steps_per_cycle = 1;
		cycles = 1;
	}

	template <typename Type_ = Type>
	SinOscillation(Type_ oscillation_, int steps_per_cycle_, int cycles_) :
		oscillation(oscillation_), steps_per_cycle(steps_per_cycle_), cycles(cycles_)
	{}

	//For the case of COSOSC3 (i.e. SinOscillation<DBL3>) we also want to have a 5 parameter constructor as below, but this will not be available for any other types
	template <typename Type_ = Type, std::enable_if_t<std::is_same<Type_, DBL3>::value>* = nullptr>
	SinOscillation(double ox, double oy, double oz, int steps_per_cycle_, int cycles_) :
		oscillation(DBL3(ox, oy, oz)), steps_per_cycle(steps_per_cycle_), cycles(cycles_)
	{}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	SinOscillation(const SinOscillation &copyThis) { oscillation = copyThis.oscillation; steps_per_cycle = copyThis.steps_per_cycle; cycles = copyThis.cycles; }

	//assignment operator
	SinOscillation& operator=(const SinOscillation &rhs) { oscillation = rhs.oscillation; steps_per_cycle = rhs.steps_per_cycle; cycles = rhs.cycles; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const SinOscillation &rhs) { os << ToString(rhs.oscillation) << "; " << rhs.steps_per_cycle << "; " << rhs.cycles; return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const SinOscillation &rhs) {

		lhs << rhs.oscillation << std::string("; ") << ToString(rhs.steps_per_cycle) << std::string("; ") << ToString(rhs.cycles);

		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend SinOscillation& operator>>(const std::stringstream &ss, SinOscillation &rhs)
	{
		//normally the values in std::string representation are as: "oscillation; steps_per_cycle; cycles"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "oscillation steps_per_cycle cycles"
		if (components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings
		switch (components.size()) {

		case 3:			//Type is e.g. a double
			rhs.oscillation = ToNum(trimspaces(components[0]), "");
			rhs.steps_per_cycle = ToNum(trimspaces(components[1]), "");
			rhs.cycles = ToNum(trimspaces(components[2]), "");
			break;
		case 5:			//Type is e.g. a DBL3
			rhs.oscillation = ToNum(trimspaces(components[0]) + ", " + trimspaces(components[1]) + ", " + trimspaces(components[2]), "");
			rhs.steps_per_cycle = ToNum(trimspaces(components[3]), "");
			rhs.cycles = ToNum(trimspaces(components[4]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	Type value(int stepIdx)
	{
		return oscillation * sin(6.2831853071795864769252867666 * (double)stepIdx / steps_per_cycle);
	}

	int number_of_steps(void) { return steps_per_cycle * cycles; }
};

typedef SinOscillation<double> SINOSC;
typedef SinOscillation<DBL3> SINOSC3;

////////////////////////////////////////////////////////////////////////////////////////////////// Cos function oscillation on a bias.
//
// 

template <typename Type>
class CosSequence {

	//----------------------------- DATA

private:

	//the bias value (e.g. Hbias field)
	Type bias;
	//the oscillation value : specify the starting point of the cos function; the cos sequence will oscillate between +oscillation and -oscillation (e.g. this is the Hrf field which should be specified to be orthogonal to the bias).
	Type oscillation;
	//to adjust the frequency the time taken per step must be defined externally, this is of no concern here (e.g. 1 GHz with 20 steps per cycle : increment steps every 50 ps)
	int steps_per_cycle, cycles;

public:

	//----------------------------- VALUE CONSTRUCTORS

	template <typename Type_ = Type>
	CosSequence(void) {

		bias = Type();
		oscillation = Type();
		steps_per_cycle = 1;
		cycles = 1;
	}

	template <typename Type_ = Type>
	CosSequence(Type_ bias_, Type_ oscillation_, int steps_per_cycle_, int cycles_) :
		bias(bias_), oscillation(oscillation_), steps_per_cycle(steps_per_cycle_), cycles(cycles_)
	{}

	//For the case of COSSEQ3 (i.e. CosSequence<DBL3>) we also want to have a 8 parameter constructor as below, but this will not be available for any other types
	template <typename Type_ = Type, std::enable_if_t<std::is_same<Type_, DBL3>::value>* = nullptr>
	CosSequence(double bx, double by, double bz, double ox, double oy, double oz, int steps_per_cycle_, int cycles_) :
		bias(DBL3(bx, by, bz)), oscillation(DBL3(ox, oy, oz)), steps_per_cycle(steps_per_cycle_), cycles(cycles_)
	{}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	CosSequence(const CosSequence &copyThis) { bias = copyThis.bias; oscillation = copyThis.oscillation; steps_per_cycle = copyThis.steps_per_cycle; cycles = copyThis.cycles; }

	//assignment operator
	CosSequence& operator=(const CosSequence &rhs) { bias = rhs.bias; oscillation = rhs.oscillation; steps_per_cycle = rhs.steps_per_cycle; cycles = rhs.cycles; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const CosSequence &rhs) { os << ToString(rhs.bias) << "; " << ToString(rhs.oscillation) << "; " << rhs.steps_per_cycle << "; " << rhs.cycles; return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const CosSequence &rhs) {

		lhs << rhs.bias << std::string("; ") << rhs.oscillation << std::string("; ") << ToString(rhs.steps_per_cycle) << std::string("; ") << ToString(rhs.cycles);

		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend CosSequence& operator>>(const std::stringstream &ss, CosSequence &rhs)
	{
		//normally the values in std::string representation are as: "bias; oscillation; steps_per_cycle; cycles"
		std::vector<std::string> components = split(ss.str(), ",", ";");
		//it could be they are also given as: "bias oscillation steps_per_cycle cycles"
		if (components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings
		switch (components.size()) {

		case 4:			//Type is e.g. a double
			rhs.bias = ToNum(trimspaces(components[0]), "");
			rhs.oscillation = ToNum(trimspaces(components[1]), "");
			rhs.steps_per_cycle = ToNum(trimspaces(components[2]), "");
			rhs.cycles = ToNum(trimspaces(components[3]), "");
			break;
		case 8:			//Type is e.g. a DBL3
			rhs.bias = ToNum(trimspaces(components[0]) + ", " + trimspaces(components[1]) + ", " + trimspaces(components[2]), "");
			rhs.oscillation = ToNum(trimspaces(components[3]) + ", " + trimspaces(components[4]) + ", " + trimspaces(components[5]), "");
			rhs.steps_per_cycle = ToNum(trimspaces(components[6]), "");
			rhs.cycles = ToNum(trimspaces(components[7]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	Type value(int stepIdx)
	{
		return oscillation * cos(6.2831853071795864769252867666 * (double)stepIdx / steps_per_cycle) + bias;
	}

	int number_of_steps(void) { return steps_per_cycle * cycles; }
};

typedef CosSequence<double> COSSEQ;
typedef CosSequence<DBL3> COSSEQ3;

////////////////////////////////////////////////////////////////////////////////////////////////// String Sequence
//
// Has start and end points, with an integer number of steps between them. Minimum of one step (in this case the sequence has just the start and end points)

class StringSequence
{

private:

	std::string numsteps_and_text, textOnly;
	int steps;

public:

	//----------------------------- VALUE CONSTRUCTORS

	StringSequence(void)
	{
		steps = 1;
	}

	//input std::string specified as "n: ...", where n is the number of sequence outputs ( = steps + 1). Each step output has the same std::string output.
	StringSequence(std::string string_seq)
	{
		size_t pos = string_seq.find_first_of(':');

		numsteps_and_text = string_seq;
		steps = 0;

		if (pos != std::string::npos) {

			steps = (int)ToNum(string_seq.substr(0, pos)) - 1;
			if (steps < 0) steps = 0;

			textOnly = string_seq.substr(pos + 1);
		}
	}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	StringSequence(const StringSequence &copyThis)
	{
		numsteps_and_text = copyThis.numsteps_and_text;
		textOnly = copyThis.textOnly;
		steps = copyThis.steps;
	}

	//assignment operator
	StringSequence& operator=(const StringSequence &rhs)
	{
		numsteps_and_text = rhs.numsteps_and_text;
		textOnly = rhs.textOnly;
		steps = rhs.steps;
		return *this;
	}

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const StringSequence &rhs)
	{
		os << rhs.numsteps_and_text;
		return os;
	}

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const StringSequence &rhs)
	{
		lhs << rhs.numsteps_and_text;
		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend StringSequence& operator>>(const std::stringstream &ss, StringSequence &rhs)
	{
		size_t pos = ss.str().find_first_of(':');

		rhs.numsteps_and_text = ss.str();
		rhs.steps = 0;

		if (pos != std::string::npos) {

			rhs.steps = (int)ToNum(ss.str().substr(0, pos)) - 1;
			if (rhs.steps < 0) rhs.steps = 0;

			rhs.textOnly = ss.str().substr(pos + 1);
		}

		return rhs;
	}

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	std::string value(int stepIdx)
	{
		return textOnly;
	}

	int number_of_steps(void) { return steps; }
};

////////////////////////////////////////////////////////////////////////////////////////////////// File-loaded Sequence
//
// Data for the sequence is loaded from a file. The file must have columns as follows:
// First column : time (s). Following columns must contain the actual data : one column for scalar quantities, 3 columns for vector data.
//

template <typename Type>
class FileSequence
{

private:

	//the directory location of fileName
	std::string directory;

	//the file name, but not including a directory, as this could change.
	//the directory must be provided externally when creating or reloading the object.
	std::string fileName;

	//the time resolution at which to load data from file
	double time_resolution;

	//output values : number of values here is the number of steps, and these are loaded form the file at the given time resolution
	std::vector<Type> values;

private:

	bool load_values_from_file(void)
	{
		if (!time_resolution) return false;

		std::vector<std::vector<double>> load_arrays;

		if (std::is_same<Type, double>::value) {

			if (ReadDataColumns(directory + fileName, "\t", load_arrays, { 0, 1 })) {

				//all required file columns must have same number of entries
				if (load_arrays[0].size() != load_arrays[1].size()) return false;

				double time_duration = load_arrays[0].back() - load_arrays[0].front();

				int steps = time_duration / time_resolution;
				if (steps <= 0) return false;

				if (!malloc_vector(values, steps)) return false;

				int file_array_index = 0;

				for (int step = 0; step < steps; step++) {

					double time = load_arrays[0].front() + step * time_resolution;

					//set file_array_index so the current time is between points at file_array_index and file_array_index + 1
					if (file_array_index + 1 < load_arrays[0].size()) {

						while (time >= load_arrays[0][file_array_index + 1]) {

							file_array_index++;
							if (file_array_index + 1 >= load_arrays[0].size()) break;
						}
					}

					if (file_array_index + 1 < load_arrays[0].size()) {

						//use interpolation if possible
						double dT = load_arrays[0][file_array_index + 1] - load_arrays[0][file_array_index];
						if (dT) {

							double parameter = (time - load_arrays[0][file_array_index]) / dT;
							double start = load_arrays[1][file_array_index];
							double end = load_arrays[1][file_array_index + 1];
							double value = parametric_interpolation(start, end, parameter);
							values[step] = *reinterpret_cast<Type*>(&value);
						}
						else {
							double value = load_arrays[1][file_array_index];
							values[step] = *reinterpret_cast<Type*>(&value);
						}
					}
					else {
						double value = load_arrays[1].back();
						values[step] = *reinterpret_cast<Type*>(&value);
					}
				}
			}
		}
		else if (std::is_same<Type, DBL3>::value) {
		
			if (ReadDataColumns(directory + fileName, "\t", load_arrays, { 0, 1, 2, 3 })) {

				//all required file columns must have same number of entries
				if (load_arrays[0].size() != load_arrays[1].size() ||
					load_arrays[0].size() != load_arrays[2].size() ||
					load_arrays[0].size() != load_arrays[3].size()) return false;

				double time_duration = load_arrays[0].back() - load_arrays[0].front();

				int steps = time_duration / time_resolution;
				if (steps <= 0) return false;

				if (!malloc_vector(values, steps)) return false;

				int file_array_index = 0;

				for (int step = 0; step < steps; step++) {

					double time = load_arrays[0].front() + step * time_resolution;

					//set file_array_index so the current time is between points at file_array_index and file_array_index + 1
					if (file_array_index + 1 < load_arrays[0].size()) {

						while (time >= load_arrays[0][file_array_index + 1]) {

							file_array_index++;
							if (file_array_index + 1 >= load_arrays[0].size()) break;
						}
					}

					if (file_array_index + 1 < load_arrays[0].size()) {

						//use interpolation if possible
						double dT = load_arrays[0][file_array_index + 1] - load_arrays[0][file_array_index];
						if (dT) {
						
							double parameter = (time - load_arrays[0][file_array_index]) / dT;
							DBL3 start = DBL3(load_arrays[1][file_array_index], load_arrays[2][file_array_index], load_arrays[3][file_array_index]);
							DBL3 end = DBL3(load_arrays[1][file_array_index + 1], load_arrays[2][file_array_index + 1], load_arrays[3][file_array_index + 1]);
							DBL3 value = parametric_interpolation(start, end, parameter);
							values[step] = *reinterpret_cast<Type*>(&value);
						}
						else {

							DBL3 value = DBL3(load_arrays[1][file_array_index], load_arrays[2][file_array_index], load_arrays[3][file_array_index]);
							values[step] = *reinterpret_cast<Type*>(&value);
						}
					}
					else {
						DBL3 value = DBL3(load_arrays[1].back(), load_arrays[2].back(), load_arrays[3].back());
						values[step] = *reinterpret_cast<Type*>(&value);
					}
				}
			}
		}

		return true;
	}

public:


	//----------------------------- VALUE CONSTRUCTORS

	FileSequence(void) {}

	FileSequence(std::string directory, std::string fileName, double time_resolution)
	{
		this->directory = directory;
		this->fileName = fileName;
		this->time_resolution = time_resolution;

		//now attempt to load values
		load_values_from_file();
	}

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	FileSequence(const FileSequence &copyThis) { *this = copyThis; }

	//assignment operator
	FileSequence& operator=(const FileSequence &rhs)
	{
		fileName = rhs.fileName;
		time_resolution = rhs.time_resolution;
		values = rhs.values;

		return *this;
	}

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const FileSequence &rhs)
	{
		os << rhs.fileName;
		return os;
	}

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const FileSequence &rhs)
	{
		lhs << rhs.fileName;
		return lhs;
	}

	//allows conversions from std::string to Sequence
	friend FileSequence& operator>>(const std::stringstream &ss, FileSequence &rhs)
	{
		//time step and directory are not set here. instead these should either already be set, or else to set values you must set them later with a dedicated call.
		rhs.fileName = ss.str();
		rhs.load_values_from_file();

		return rhs;
	}

	//----------------------------- SETTERS

	bool set_filename(std::string directory, std::string fileName)
	{
		this->directory = directory;
		this->fileName = fileName;

		//now attempt to load values
		return load_values_from_file();
	}

	bool set_time_resolution(double time_resolution)
	{
		this->time_resolution = time_resolution;

		//now attempt to load values
		return load_values_from_file();
	}

	bool set_filename_and_time_resolution(std::string directory, std::string fileName, double time_resolution)
	{
		this->directory = directory;
		this->fileName = fileName;
		this->time_resolution = time_resolution;

		//now attempt to load values
		return load_values_from_file();
	}

	//----------------------------- GETTERS

	std::string get_fileName(void) { return fileName; }

	//----------------------------- OTHERS

	//return value for given step index (ranges from 0 to steps)
	Type value(int stepIdx) 
	{ 
		if (stepIdx < values.size()) return values[stepIdx];
		else return Type();
	}

	int number_of_steps(void) { return values.size() - 1; }
};

typedef FileSequence<double> FILESEQ;
typedef FileSequence<DBL3> FILESEQ3;

////////////////////////////////////////////////////////////////////////////////////////////////// exclusions
//
//

template <typename Type>
struct exclusions {

private:

	std::vector< std::vector<Type> > sets;

	std::vector<Type> empty_set;

public:

	exclusions(void) {}

	std::vector<Type>& operator[](const Type& value)
	{
		for (int idx = 0; idx < (int)sets.size(); idx++) {

			if (vector_contains(sets[idx], value)) return sets[idx];
		}

		return empty_set;
	}

	void clear(void) { sets.resize(0); }

	void storeset(std::vector<Type> new_set) { sets.push_back(new_set); }
	template <typename ... VType> void storeset(VType ... values) { sets.push_back(make_vector(values...)); }
};
