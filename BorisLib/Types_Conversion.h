#pragma once

#include <string>
#include <sstream>
#include <fstream>

#include "Introspection_base.h"
#include "Funcs_Aux_base.h"
#include "Funcs_Math_base.h"

////////////////////////////////////////////////////////////////////////////////////////////////// CONVERSIONS - used for converting types to/from strings, allowing units to be used
//
//

namespace Conversion {

	template <typename Type>
	std::string ToString_convertible(const Type& value, std::true_type)
	{
		//called if the operation ss << value is possible (Type has is_streamable_out trait for std::stringstream)
		std::stringstream ss;
		ss << value;

		return ss.str();
	}

	template <typename Type>
	std::string ToString_convertible(const Type& value, std::false_type)
	{
		return "";
	}

	//-----------------------------------------

	class tostringconversion {

	private:

		std::string text;
		std::string unit;

	private:

		//this is called if value can be streamed to a std::stringstream
		//
		template <typename Type> 
		std::string convert(const Type& value, std::true_type)
		{
			std::stringstream ss;

			if (!unit.length()) { ss << value; return ss.str(); }

			//Note! Unit should only be set for types that can be converted to a double. Need the reinterpret_cast to stop compilation errors.
			int decexp;
			double value_adjusted = decexp_eng(value, &decexp);

			std::string unitmagnitude;

			switch (decexp) {

			case -18:
				unitmagnitude = "a";
				break;
			case -15:
				unitmagnitude = "f";
				break;
			case -12:
				unitmagnitude = "p";
				break;
			case -9:
				unitmagnitude = "n";
				break;
			case -6:
				unitmagnitude = "u";
				break;
			case -3:
				unitmagnitude = "m";
				break;
			case 0:
				unitmagnitude = "";
				break;
			case 3:
				unitmagnitude = "k";
				break;
			case 6:
				unitmagnitude = "M";
				break;
			case 9:
				unitmagnitude = "G";
				break;
			case 12:
				unitmagnitude = "T";
				break;
			case 15:
				unitmagnitude = "P";
				break;
			default:
				unitmagnitude = "";
				break;
			}

			ss << value_adjusted << unitmagnitude + unit;

			return ss.str();
		}

		std::string convert(const std::string& rhs, std::true_type) { return rhs; }

		template <typename Type>
		std::string convert(const Type& value, std::false_type) { return ""; }

	public:

		tostringconversion(void) {}
	
		void set_unit(const std::string& unit_) { unit = unit_; }

		template <typename Type> 
		tostringconversion& operator<<(const Type& rhs)
		{ 
			text += convert(rhs, is_streamable_out<std::stringstream, Type>());

			return *this; 
		}

		std::string str(void) { return text; }
	};

	//-----------------------------------------

	class tonumberconversion {

	private:

		std::string text;
		std::string unit;

	private:

		template <typename Type>
		Type convertible(std::true_type)
		{
			Type value;

			//if unit specified (do not leave space between number of unit) then replace it accordingly. Note, str may contain multiple numbers, hence multiple conversions required (e.g. a DBL3)
			if (unit.length()) {

				replaceall(text, "a" + unit, "e-18");
				replaceall(text, "f" + unit, "e-15");
				replaceall(text, "p" + unit, "e-12");
				replaceall(text, "n" + unit, "e-9");
				replaceall(text, "u" + unit, "e-6");
				replaceall(text, "m" + unit, "e-3");
				replaceall(text, "k" + unit, "e+3");
				replaceall(text, "M" + unit, "e+6");
				replaceall(text, "G" + unit, "e+9");
				replaceall(text, "T" + unit, "e+12");
				replaceall(text, "P" + unit, "e+15");
				replaceall(text, unit, "");
			}

			std::stringstream ss(text);
			//>> operator is overloaded for any class that needs conversion from std::string as friend Type& operator>>(const std::stringstream &ss, Type &rhs) {...}
			ss >> value;
			return value;
		}

		template <typename Type>
		Type convertible(std::false_type) { return Type(); }

	public:

		tonumberconversion(const std::string& text_, const std::string& unit_ = "") : text(text_), unit(unit_) {}
		tonumberconversion(std::string&& text_, const std::string& unit_ = "") : text( move(text_) ), unit(unit_) {}

		//conversion operator from Conversion to Type
		template <typename Type> 
		operator Type() 
		{
			//can only convert if Type accepts >> operator from strinstream
			return convertible<Type>(is_streamable_in<std::stringstream, Type>());
		}

		//operator std::string&() { return move(text); }
		operator std::string&() { return text; }
	};
};