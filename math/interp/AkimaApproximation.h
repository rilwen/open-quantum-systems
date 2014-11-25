#ifndef __MATH_INTERP_AKIMA_APPROXIMATION_H
#define __MATH_INTERP_AKIMA_APPROXIMATION_H

#include <vector>
#include <stdexcept>
#include "../IndexShifter.h"
#include <cassert>

namespace rql {
	namespace interp {
		/** Approximates dy/dx given x and y(x) using the Akima algorithm
		@tparam ValueType Function value type
		*/
		template <class ValueType>
		struct AkimaApproximation
		{
			//! @tparam V1 Vector class with operator[] and .size()
			//! @tparam V2 Vector class with operator[] and .size()
			//! @tparam V3 Vector class with operator[] and .size()
			template <class V1, class V2, class V3>
			static void calculate(const V1& x, const V2& y, V3& dy);
		};

		template <class ValueType>
		template <class V1, class V2, class V3>
		void AkimaApproximation<ValueType>::calculate(const V1& x, const V2& y, V3& dy)
		{
			const size_t size = x.size();
			if (y.size() != size)
				throw std::domain_error("AkimaApproximation: x and y size mismatch");
			if (dy.size() != size)
				throw std::domain_error("AkimaApproximation: x and dy size mismatch");
			if (size < 2)
				throw std::domain_error("AkimaApproximation: needs at least 1 segment");
			const size_t nbr_segs = size - 1;

			ShiftedVector<std::vector<ValueType>,ValueType> d(std::vector<ValueType>(size+3), -2);
			for (size_t i = 0; i < nbr_segs; ++i)
				d[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
			d[-1] = 2*d[0]-d[1];
			d[-2] = 2*d[-1]-d[0];
			d[nbr_segs] = 2*d[nbr_segs-1]-d[nbr_segs-2];
			d[nbr_segs+1] = 2*d[nbr_segs] - d[nbr_segs-2];

			// approximate the derivatives of y
			for (size_t i = 0; i < size; ++i)
			{
				const ValueType w_l = std::abs(d[i+1] - d[i]);
				const ValueType w_r = std::abs(d[i-1] - d[i-2]);
				if (w_l != 0 || w_r != 0) {
					dy[i] = (w_l*d[i-1] + w_r*d[i])/(w_l + w_r);
				} else {
					dy[i] = 0.5*(d[i - 1] + d[i]);
				}
			}
		}
	}
}

#endif // __MATH_INTERP_AKIMA_APPROXIMATION_H
