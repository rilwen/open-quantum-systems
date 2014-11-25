#pragma once
#include <vector>
#include <cassert>
#include <ostream>
#include "Segment.h"
#include "MathCore.h"

namespace rql { namespace math {

	class BorelSet
	{
	public:
		template <class Iter> 
		BorelSet(Iter begin, const Iter end, bool complete = false);
		RQL_MATH_API_DEBUG BorelSet(const Segment& segment = Segment::EMPTY);

		RQL_MATH_API_DEBUG size_t nbr_segments() const
		{
			return segments_.size();
		}
	
		RQL_MATH_API_DEBUG const Segment& segment(size_t i) const
		{
	#ifdef NDEBUG
			return segments_[i];
	#else
			return segments_.at(i);
	#endif
		}

		RQL_MATH_API_DEBUG BorelSet complete() const;

		RQL_MATH_API_DEBUG bool is_empty() const;

		// set product
		RQL_MATH_API_DEBUG BorelSet operator*(const BorelSet& other) const;
		// set union
		RQL_MATH_API_DEBUG BorelSet operator+(const BorelSet& other) const;
		// adjoint set
		RQL_MATH_API_DEBUG friend BorelSet operator!(const BorelSet& other);

		RQL_MATH_API_DEBUG bool operator==(const BorelSet& other) const;
		RQL_MATH_API_DEBUG bool operator!=(const BorelSet& other) const;

		RQL_MATH_API_DEBUG friend std::ostream& operator<<(std::ostream& os, const BorelSet& bs);

		RQL_MATH_API_DEBUG double lim_inf() const;
		RQL_MATH_API_DEBUG double lim_sup() const;
	private:
		std::vector<Segment> segments_;
	};

	template <class Iter>
	BorelSet::BorelSet(Iter i, const Iter end, bool complete)
	{
		segments_.push_back(Segment::EMPTY); // ensures we have something
		for (; i != end; ++i) {
			const Segment s = complete ? i->complete() : *i;
			if (segments_.rbegin()->can_be_merged_with(s, complete))
				*segments_.rbegin() = segments_.rbegin()->merge_with(s, complete);
			else
				segments_.push_back(s);
		}
		assert(segments_.size() > 0);
	}

}}