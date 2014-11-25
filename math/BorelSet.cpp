#include "BorelSet.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

namespace rql { namespace math {

	BorelSet::BorelSet(const Segment& segment)
	{
		segments_.push_back(segment);
	}

	BorelSet BorelSet::complete() const
	{
		return BorelSet(segments_.begin(), segments_.end(), true);
	}

	bool BorelSet::is_empty() const
	{
		assert(segments_.size() > 0);
		return segments_[0].is_empty();
	}

	BorelSet BorelSet::operator +(const BorelSet& o) const
	{
		std::vector<Segment> tmp;
		tmp.reserve(nbr_segments() + o.nbr_segments());
		for (std::vector<Segment>::const_iterator i = segments_.begin(); i != segments_.end(); ++i)
			tmp.push_back(*i);
		for (std::vector<Segment>::const_iterator i = o.segments_.begin(); i != o.segments_.end(); ++i)
			tmp.push_back(*i);
		std::sort(tmp.begin(), tmp.end());
		return BorelSet(tmp.begin(), tmp.end());
	}

	BorelSet BorelSet::operator *(const BorelSet& o) const
	{
		std::vector<Segment> tmp;
		tmp.reserve(nbr_segments() * o.nbr_segments());
		for (std::vector<Segment>::const_iterator i = segments_.begin(); i != segments_.end(); ++i)
			for (std::vector<Segment>::const_iterator j = o.segments_.begin(); j != o.segments_.end(); ++j)
				tmp.push_back((*i) * (*j));
		std::sort(tmp.begin(), tmp.end());
		return BorelSet(tmp.begin(), tmp.end());
	}

	BorelSet operator!(const BorelSet& bs)
	{
		std::vector<Segment> tmp;
		tmp.reserve(2);
		bool frist = true;
		BorelSet result;
		for (std::vector<Segment>::const_iterator j = bs.segments_.begin(); j != bs.segments_.end(); ++j) {
			tmp.clear();
			if (j->lower() > -std::numeric_limits<double>::infinity() || (!j->lower_inclusive()))
				tmp.push_back(Segment(-std::numeric_limits<double>::infinity(), j->lower(), true, !j->lower_inclusive()));
			if (j->upper() < std::numeric_limits<double>::infinity() || (!j->upper_inclusive()))
				tmp.push_back(Segment(j->upper(), std::numeric_limits<double>::infinity(), !j->upper_inclusive(), true));
			result = frist ? BorelSet(tmp.begin(), tmp.end()) : result * BorelSet(tmp.begin(), tmp.end());
		}
		return result;
	}

	bool BorelSet::operator ==(const BorelSet& other) const
	{
		if (this == &other)
			return true;
		if (nbr_segments() != other.nbr_segments())
			return false;
		for (size_t i = 0; i < nbr_segments(); ++i) {
			if (segments_[i] != other.segments_[i])
				return false;
		}
		return true;
	}

	bool BorelSet::operator !=(const BorelSet& other) const
	{
		return !(*this == other);
	}

	std::ostream& operator<<(std::ostream& os, const BorelSet& bs)
	{
		os << "[" << bs.segment(0);
		for (size_t i = 1; i < bs.nbr_segments(); ++i)
			os << "," << bs.segment(i);
		os << "]";
		return os;
	}

	double BorelSet::lim_inf() const
	{
		assert(segments_.size() > 0);
		return segments_[0].lower();
	}

	double BorelSet::lim_sup() const
	{
		assert(segments_.size() > 0);
		return segments_[segments_.size() - 1].upper();
	}

}}