#include "InterpolatorImplPiecewise.h"
#include "../SegmentSearch.h"
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cassert>


namespace rql {
	namespace interp {		
		InterpolatorImplPiecewise::InterpolatorImplPiecewise(const std::vector<double>& x, bool leftInclusive)
			: m_x(x), m_left_inclusive(leftInclusive)
		{
			setup();
		}

		double InterpolatorImplPiecewise::evaluate(double x) const
		{
			if (x < m_x.front() || x > m_x.back())
				throw std::out_of_range("X outside range");

			int segment_idx;
			if (m_left_inclusive)
			{
				segment_idx = std::min<int>(SegmentSearch::binary_search_left_inclusive(m_x, x), m_nbr_segments - 1);
				assert( segment_idx >= 0 );
			}
			else
			{
				segment_idx = std::max<int>(SegmentSearch::binary_search_right_inclusive(m_x, x), 0);
			}
			return evaluate(x, segment_idx);
		}

		double InterpolatorImplPiecewise::lowerBound() const
		{
			return m_x.front();
		}

		double InterpolatorImplPiecewise::upperBound() const
		{
			return m_x.back();
		}		

		bool InterpolatorImplPiecewise::left_inclusive() const
		{
			return m_left_inclusive;
		}
	}
}