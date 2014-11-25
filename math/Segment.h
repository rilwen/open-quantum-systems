#pragma once

#include <ostream>
#include "MathCore.h"

namespace rql { namespace math {

	struct Segment
	{
	public:
		RQL_MATH_API_DEBUG Segment(double lower, double upper, bool lower_inclusive, bool upper_inclusive);
		RQL_MATH_API_DEBUG Segment(double value); // create a segment with delta support
		RQL_MATH_API_DEBUG Segment(const Segment& other);
	
		RQL_MATH_API_DEBUG double lower() const
		{
			return l_;
		}
	
		RQL_MATH_API_DEBUG double upper() const
		{
			return u_;
		}
	
		RQL_MATH_API_DEBUG bool lower_inclusive() const
		{
			return li_;
		}
	
		RQL_MATH_API_DEBUG bool upper_inclusive() const
		{
			return ui_;
		}
	
		RQL_MATH_API_DEBUG bool is_empty() const
		{
			return l_ == u_ && !(li_ || ui_);
		}

		// liminf of the distance between points from this and other segment
		RQL_MATH_API_DEBUG double min_distance(const Segment& other) const;

		RQL_MATH_API_DEBUG Segment& operator=(const Segment& other);
	
		RQL_MATH_API_DEBUG bool operator==(const Segment& other) const;
		RQL_MATH_API_DEBUG bool operator!=(const Segment& other) const;
	
		// "a < b" <==> a.lower() < b.lower() OR (a.lower() == b.lower() AND a.upper() < b.upper())
		// OR (a.lower() == b.lower() AND a.upper() == b.upper() AND a.lower_inclusive() AND !b.lower_inclusive())
		// OR (a.lower() == b.lower() AND a.upper() == b.upper() AND a.lower_inclusive() AND b.lower_inclusive() AND !a.lower_inclusive() AND b.upper_inclusive())
		// empty segment is lower than any other segment
		RQL_MATH_API_DEBUG bool operator<(const Segment& other) const;
		RQL_MATH_API_DEBUG bool operator<=(const Segment& other) const;

		friend RQL_MATH_API_DEBUG std::ostream& operator<<(std::ostream& os, const Segment& s);

		// is this segment totally below the other segment (any point from this segment is below all points from the other segment)?
		// ill-defined if any of the segments is empty
		RQL_MATH_API_DEBUG bool is_below(const Segment& other) const;
	
		// Can this segment be merged with other segment into a new Segment
		// ignore_missing_point -- return true even if there is a singe point between the segments which belongs to neither of them
		RQL_MATH_API_DEBUG bool can_be_merged_with(const Segment& other, bool ignore_missing_point = false) const;

		// Merge with other segment, assuming it's possible.
		// ignore_missing_point -- return true even if there is a singe point between the segments which belongs to neither of them
		RQL_MATH_API_DEBUG Segment merge_with(const Segment& other, bool ignore_missing_point = false) const;

		// return a completion (inclusive at both ends)
		RQL_MATH_API_DEBUG Segment complete() const;

		// set product
		friend RQL_MATH_API_DEBUG Segment operator*(const Segment& s1, const Segment& s2);

		RQL_MATH_API_DEBUG const static Segment EMPTY;
		RQL_MATH_API_DEBUG const static double INF;
	private:
		Segment(); // create empty segment
		static Segment create(double new_lower, double new_upper, bool new_lower_inclusive, bool new_upper_inclusive);
		double l_;
		double u_;
		bool li_;
		bool ui_;
	};

	RQL_MATH_API_DEBUG Segment operator*(const Segment& s1, const Segment& s2);

}}