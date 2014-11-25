#include "InterpolatorImplConstant.h"
#include <limits>

namespace rql {
	namespace interp {
		InterpolatorImplConstant::InterpolatorImplConstant(double constant)
			: m_const(constant), m_lower_bound(-std::numeric_limits<double>::infinity())
		{
		}

		InterpolatorImplConstant::InterpolatorImplConstant(double constant, double lowerBound)
			: m_const(constant), m_lower_bound(lowerBound)
		{
		}

		double InterpolatorImplConstant::evaluate(double) const
		{
			return m_const;
		}

		double InterpolatorImplConstant::lowerBound() const
		{
			return m_lower_bound;
		}

		double InterpolatorImplConstant::upperBound() const
		{
			return std::numeric_limits<double>::infinity();
		}

		std::shared_ptr<InterpolatorImpl> InterpolatorImplConstant::clone() const
		{
			return std::shared_ptr<InterpolatorImplConstant>(new InterpolatorImplConstant(m_const));
		}

		InterpolatorImpl& InterpolatorImplConstant::operator+=(double x)
		{
			m_const += x;
			return *this;
		}

		InterpolatorImpl& InterpolatorImplConstant::operator-=(double x)
		{
			m_const -= x;
			return *this;
		}

		InterpolatorImpl& InterpolatorImplConstant::operator*=(double x)
		{
			m_const *= x;
			return *this;
		}

		InterpolatorImpl& InterpolatorImplConstant::operator/=(double x)
		{
			m_const /= x;
			return *this;
		}
	}
}