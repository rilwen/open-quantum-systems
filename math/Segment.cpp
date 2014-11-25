#include "Segment.h"
#include <algorithm>
#include <cassert>
#include <vector>
#include <limits>
#include <cstdio>

namespace rql { namespace math {

	const Segment Segment::EMPTY = Segment();
	const double Segment::INF = std::numeric_limits<double>::infinity();

	Segment::Segment(double lower, double upper, bool lower_inclusive, bool upper_inclusive)
	: l_(lower), u_(upper), li_(lower_inclusive), ui_(upper_inclusive)
	{
	assert( lower <= upper );
	assert( lower < upper || (lower_inclusive == upper_inclusive) );
	}

	Segment::Segment(double value)
	: l_(value), u_(value), li_(true), ui_(true)
	{
	}

	Segment::Segment()
	: l_(0.0), u_(0.0), li_(false), ui_(false)
	{
	}

	Segment::Segment(const Segment& o)
	: l_(o.l_), u_(o.u_), li_(o.li_), ui_(o.ui_)
	{
	}

	Segment& Segment::operator=(const Segment& o)
	{
		if (this != &o) {
			l_ = o.l_;
			u_ = o.u_;
			li_ = o.li_;
			ui_ = o.ui_;
		}
		return *this;
	}

	bool Segment::operator==(const Segment& o) const
	{
		if (this == &o)
			return true;
		if (is_empty()) {
			return o.is_empty();
		} else {
			return l_ == o.l_ && u_ == o.u_ && li_ == o.li_ && ui_ == o.ui_;
		}
	}

	bool Segment::operator!=(const Segment& o) const
	{
		return !(*this==o);
	}

	bool less_than_helper(const std::vector<double>& diffs)
	{
		for (std::vector<double>::const_iterator i = diffs.begin(); i != diffs.end(); ++i) {
			if (*i < 0)
				return true;
			else if (*i > 0)
				return false;
		}
		return false;
	}

	double flag2double(bool flag)
	{
		return flag ? 1.0 : 0.0;
	}

	bool Segment::operator<(const Segment& other) const
	{
		if (is_empty())
			return true;
		else if (other.is_empty())
			return false;
		std::vector<double> diffs(4);
		diffs[0] = l_ - other.l_;
		diffs[1] = u_ - other.u_;
		diffs[2] = flag2double(other.li_) - flag2double(li_);
		diffs[3] = flag2double(ui_) - flag2double(other.ui_);
		return less_than_helper(diffs);
	}

	bool Segment::operator<=(const Segment& other) const
	{
		return !(other < *this);
	}

	double Segment::min_distance(const Segment& other) const
	{
		if (u_ <= other.l_)
			return other.l_ - u_;
		else if (other.u_ <= l_)
			return l_ - other.u_;
		else
			return 0.0;
	}

	bool Segment::is_below(const Segment& other) const
	{
		return u_ < other.l_ || (u_ == other.l_ && !(ui_ && other.li_));
	}

	bool Segment::can_be_merged_with(const Segment& other, bool ignore_missing_point) const
	{
		if (is_empty() || other.is_empty())
			return true;
		if (this->min_distance(other) != 0.0)
			return false;
		if (ignore_missing_point)
			return true;
		if (u_ == other.l_)
			return ui_ || other.li_;
		if (l_ == other.u_)
			return li_ || other.ui_;
		return true;
	}

	Segment Segment::merge_with(const Segment& o, bool /*ignore_missing_point*/) const
	{
		//assert( can_be_merged_with(o, ignore_missing_point) );
		if (is_empty())
			return o;
		else if (o.is_empty())
			return *this;
		else {
			double new_lower, new_upper;
			bool new_lower_inclusive, new_upper_inclusive;
			if (l_ < o.l_) {
				new_lower = l_;
				new_lower_inclusive = li_;
			} else if (o.l_ < l_) {
				new_lower = o.l_;
				new_lower_inclusive = o.li_;
			} else {
				new_lower = l_;
				new_lower_inclusive = li_ || o.li_;
			}
			if (u_ > o.u_) {
				new_upper = u_;
				new_upper_inclusive = ui_;
			} else if (o.u_ > u_) {
				new_upper = o.u_;
				new_upper_inclusive = o.ui_;
			} else {
				new_upper = u_;
				new_upper_inclusive = ui_ || o.ui_;
			}
			return Segment::create(new_lower, new_upper, new_lower_inclusive, new_upper_inclusive);
		}
	}

	Segment operator*(const Segment& s1, const Segment& s2)
	{
		double new_lower, new_upper;
		bool new_lower_inclusive, new_upper_inclusive;
		if (s1.lower() < s2.lower()) {
			new_lower = s2.lower();
			new_lower_inclusive = s2.lower_inclusive();
		} else if (s1.lower() > s2.lower()) {
			new_lower = s1.lower();
			new_lower_inclusive = s1.lower_inclusive();
		} else {
			new_lower = s1.lower();
			new_lower_inclusive = s1.lower_inclusive() && s2.lower_inclusive();
		}
		if (s1.upper() > s2.upper()) {
			new_upper = s2.upper();
			new_upper_inclusive = s2.upper_inclusive();
		} else if (s1.upper() < s2.upper()) {
			new_upper = s1.upper();
			new_upper_inclusive = s1.upper_inclusive();
		} else {
			new_upper = s1.upper();
			new_upper_inclusive = s1.upper_inclusive() && s2.upper_inclusive();
		}
		return Segment::create(new_lower, new_upper, new_lower_inclusive, new_upper_inclusive);
	}

	Segment Segment::create(double new_lower, double new_upper, bool new_lower_inclusive, bool new_upper_inclusive)
	{
		if (new_lower < new_upper) {
			return Segment(new_lower, new_upper, new_lower_inclusive, new_upper_inclusive);
		} else if (new_lower == new_upper && new_lower_inclusive && new_upper_inclusive) {
			return Segment(new_lower);
		} else {
			return Segment();
		}
	}

	const char* bool2str(bool value)
	{
		return value ? "TRUE" : "FALSE";
	}

	std::ostream& operator<<(std::ostream& os, const Segment& s)
	{
		if (s.is_empty())
			os << "[]";
		else
			os << "[" << s.lower() << "," << s.upper() << "," << bool2str(s.lower_inclusive()) << "," << bool2str(s.upper_inclusive()) << "]";
		return os;
	}

	Segment Segment::complete() const
	{
		return is_empty() ? *this : Segment::create(l_, u_, true, true);
	}

}}
