#include "InterpolatorImplFactory.h"
#include "InterpolatorImplPiecewiseConstant.h"
#include "InterpolatorImplPiecewisePolynomial.h"
#include "InterpolatorImplConstant.h"
#include "AkimaApproximation.h"
#include <stdexcept>
#include <vector>

namespace rql {
	namespace interp {

		std::shared_ptr<InterpolatorImpl> InterpolatorImplFactory::constant(double y)
		{
			return std::shared_ptr<InterpolatorImpl>(new InterpolatorImplConstant(y));
		}

		std::shared_ptr<InterpolatorImpl> InterpolatorImplFactory::constant(double y, double lb)
		{
			return std::shared_ptr<InterpolatorImpl>(new InterpolatorImplConstant(y, lb));
		}

		std::shared_ptr<InterpolatorImpl> InterpolatorImplFactory::piecewiseConstant(const std::vector<double>& x, const std::vector<double>& y, bool leftInclusive)
		{
			if (x.size() != y.size() + 1)
				throw std::domain_error("InterpolatorImplFactory::piecewiseConstant: x and y size mismatch");
			return std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewiseConstant(x, y, leftInclusive));
		}

		std::shared_ptr<InterpolatorImpl> InterpolatorImplFactory::piecewiseLinear(const std::vector<double>& x, const std::vector<double>& y)
		{
			const unsigned int size = x.size();
			if (y.size() != size)
				throw std::domain_error("InterpolatorImplFactory: x and y size mismatch");
			std::vector<typename InterpolatorImplPiecewisePolynomial<1>::DataNode> data(size);
			for (unsigned int i = 0; i < size; ++i)
			{
				data[i].x() = x[i];
				data[i].y()[0] = y[i];
			}
			return std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewisePolynomial<1>(data));
		}

		std::shared_ptr<InterpolatorImpl> InterpolatorImplFactory::piecewiseCubic(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& dy)
		{
			const unsigned int size = x.size();
			if (y.size() != size)
				throw std::domain_error("InterpolatorImplFactory: x and y size mismatch");
			if (dy.size() != size)
				throw std::domain_error("InterpolatorImplFactory: x and dy size mismatch");
			std::vector<typename InterpolatorImplPiecewisePolynomial<2>::DataNode> data(size);
			for (unsigned int i = 0; i < size; ++i)
			{
				data[i].x() = x[i];
				data[i].y()[0] = y[i];
				data[i].y()[1] = dy[i];
			}
			return std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewisePolynomial<2>(data));
		}

		std::shared_ptr<InterpolatorImpl> InterpolatorImplFactory::akima(const std::vector<double>& x, const std::vector<double>& y)
		{
			const unsigned int size = x.size();
			if (y.size() != size)
				throw std::domain_error("InterpolatorImplFactory: x and y size mismatch");
			std::vector<double> dy(size);
			AkimaApproximation<double>::calculate(x, y, dy);
			return piecewiseCubic(x, y, dy);
		}
	}
}