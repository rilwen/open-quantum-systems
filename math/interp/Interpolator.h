#ifndef __MATH_INTERPOLATOR_H
#define __MATH_INTERPOLATOR_H

#include <memory>
#include "../MathCore.h"
#include "../Pimpl.h"
#include "InterpolatorImpl.h"

namespace rql {
	namespace interp {

		// PIMPL idiom enables the use of Interpolator as a functor.
		// Relies on auto-generated copy operators/constructors.
		class Interpolator: public Pimpl<InterpolatorImpl>
		{
		public:
			RQL_MATH_API_DEBUG Interpolator(std::shared_ptr<InterpolatorImpl> impl);
			inline double operator()(double x) const;
			inline double lowerBound() const;
			inline double upperBound() const;
			inline Interpolator operator+(double x) const;
			inline Interpolator operator-(double x) const;
			inline Interpolator operator*(double x) const;
			inline Interpolator operator/(double x) const;
			inline Interpolator& operator+=(double x);
			inline Interpolator& operator-=(double x);
			inline Interpolator& operator*=(double x);
			inline Interpolator& operator/=(double x);
		};

		inline double Interpolator::operator()(double x) const
		{
			return impl().evaluate(x);
		}
		
		inline double Interpolator::lowerBound() const
		{
			return impl().lowerBound();
		}

		inline double Interpolator::upperBound() const
		{
			return impl().upperBound();
		}

		inline Interpolator Interpolator::operator+(double x) const
		{
			std::shared_ptr<InterpolatorImpl> clone(impl().clone());
			*clone += x;
			return Interpolator(clone);
		}

		inline Interpolator Interpolator::operator-(double x) const
		{
			std::shared_ptr<InterpolatorImpl> clone(impl().clone());
			*clone -= x;
			return Interpolator(clone);
		}

		inline Interpolator Interpolator::operator*(double x) const
		{
			std::shared_ptr<InterpolatorImpl> clone(impl().clone());
			*clone *= x;
			return Interpolator(clone);
		}

		inline Interpolator Interpolator::operator/(double x) const
		{
			std::shared_ptr<InterpolatorImpl> clone(impl().clone());
			*clone /= x;
			return Interpolator(clone);
		}

		inline Interpolator& Interpolator::operator+=(double x)
		{
			impl() += x;
			return *this;
		}

		inline Interpolator& Interpolator::operator-=(double x)
		{
			impl() -= x;
			return *this;
		}

		inline Interpolator& Interpolator::operator*=(double x)
		{
			impl() *= x;
			return *this;
		}

		inline Interpolator& Interpolator::operator/=(double x)
		{
			impl() /= x;
			return *this;
		}
	}
}

#endif // __MATH_INTERPOLATOR_H 