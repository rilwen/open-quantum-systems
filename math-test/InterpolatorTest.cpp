#include <algorithm>
#include <limits>
#include <list>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <gtest/gtest.h>
#include "math/interp/Interpolator.h"
#include "math/interp/InterpolatorImplConstant.h"
#include "math/interp/InterpolatorImplPiecewiseConstant.h"
#include "math/interp/InterpolatorImplFactory.h"

namespace rql {
	namespace interp {

		class InterpolatorTest: public testing::Test
		{
		public:
			static void testOperations(Interpolator i);
		};

		TEST_F(InterpolatorTest,Constant)
		{
			Interpolator interpolator(std::shared_ptr<InterpolatorImpl>(new InterpolatorImplConstant(0.1)));
			EXPECT_NEAR(0.1, interpolator(45.12), 1E-20);
			EXPECT_EQ(-std::numeric_limits<double>::infinity(), interpolator.lowerBound());
			EXPECT_EQ(std::numeric_limits<double>::infinity(), interpolator.upperBound());
			testOperations(interpolator);
		}

		TEST_F(InterpolatorTest, PiecewiseConstant)
		{
			std::vector<double> x(3);
			std::vector<double> y(2);
			x[0] = 0;
			x[1] = 2;
			x[2] = 3;
			y[0] = -1;
			y[1] = 1;
			Interpolator i1(std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewiseConstant(x, y, false)));
			testOperations(i1);
			Interpolator i2(std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewiseConstant(x, y, true)));
			testOperations(i2);
			EXPECT_EQ(0, i1.lowerBound());
			EXPECT_EQ(0, i2.lowerBound());
			EXPECT_EQ(3, i1.upperBound());
			EXPECT_EQ(3, i2.upperBound());
			EXPECT_EQ(y[0], i1(x[0]));
			EXPECT_EQ(y[0], i2(x[0]));
			EXPECT_EQ(y[1], i1(x[2]));
			EXPECT_EQ(y[1], i2(x[2]));
			EXPECT_ANY_THROW(i1(x[0]-1));
			EXPECT_ANY_THROW(i2(x[0]-1));
			EXPECT_ANY_THROW(i1(x[2]+1));
			EXPECT_ANY_THROW(i2(x[2]+1));
			EXPECT_EQ(y[0], i1(x[1]));
			EXPECT_NEAR(2*y[0], (i1*2)(x[1]), 1E-16);
			EXPECT_EQ(y[1], i2(x[1]));
			EXPECT_NEAR(2*y[1], (i2*2)(x[1]), 1E-16);
			EXPECT_EQ(y[0], i1(0.5*(x[0]+x[1])));
			EXPECT_EQ(y[0], i2(0.5*(x[0]+x[1])));
			EXPECT_EQ(y[1], i1(0.5*(x[1]+x[2])));
			EXPECT_EQ(y[1], i2(0.5*(x[1]+x[2])));
			
			std::list<double> x2;
			std::copy(x.begin(), x.end(), std::back_inserter(x2));
			Interpolator i4(std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewiseConstant(x2.begin(), x2.end(), y.begin(), y.end(), true)));

			Eigen::VectorXd x3(x.size());
			Interpolator i3(std::shared_ptr<InterpolatorImpl>(new InterpolatorImplPiecewiseConstant(x3, y, false)));
		}

		TEST_F(InterpolatorTest, PiecewiseLinear)
		{
			std::vector<double> x(3);
			std::vector<double> y(3);
			x[0] = 0;
			x[1] = 2;
			x[2] = 3;
			y[0] = -1;
			y[1] = 1;
			y[2] = 0;
			Interpolator i1(std::shared_ptr<InterpolatorImpl>( InterpolatorImplFactory::piecewiseLinear(x, y) ));
			testOperations(i1);
			EXPECT_EQ(0, i1.lowerBound());
			EXPECT_EQ(3, i1.upperBound());
			EXPECT_EQ(y[0], i1(x[0]));
			EXPECT_EQ(y[1], i1(x[1]));
			EXPECT_EQ(y[2], i1(x[2]));
			EXPECT_ANY_THROW(i1(x[0]-1));
			EXPECT_ANY_THROW(i1(x[2]+1));
			EXPECT_NEAR(0, i1(0.5*(x[0]+x[1])), 1E-16);
			EXPECT_NEAR(0.5, i1(0.5*(x[1]+x[2])), 1E-16);
		}

		TEST_F(InterpolatorTest, PiecewiseCubic)
		{
			std::vector<double> x(3);
			std::vector<double> y(3);
			std::vector<double> dy(3);
			x[0] = 0;
			x[1] = 2;
			x[2] = 3;
			y[0] = -1;
			y[1] = 1;
			y[2] = 0;
			dy[0] = 0.1;
			dy[1] = -0.2;
			dy[2] = 0.5;
			Interpolator i1(std::shared_ptr<InterpolatorImpl>( InterpolatorImplFactory::piecewiseCubic(x, y, dy) ));
			testOperations(i1);
			EXPECT_EQ(0, i1.lowerBound());
			EXPECT_EQ(3, i1.upperBound());
			EXPECT_EQ(y[0], i1(x[0]));
			EXPECT_EQ(y[1], i1(x[1]));
			EXPECT_NEAR(y[2], i1(x[2]), 4E-16);
			EXPECT_ANY_THROW(i1(x[0]-1));
			EXPECT_ANY_THROW(i1(x[2]+1));			
		}

		double fun(double x)
		{
			return cos(x);
		}

		TEST_F(InterpolatorTest, Akima)
		{
			const unsigned int n = 10;
			std::vector<double> x(n);
			std::vector<double> y(n);
			for (unsigned int i = 0; i < n; ++i)
			{
				x[i] = i;
				y[i] = fun(x[i]);
			}
			Interpolator interp(std::shared_ptr<InterpolatorImpl>( InterpolatorImplFactory::akima(x, y) ));
			for (unsigned int i = 0; i < n; ++i)
			{
				EXPECT_NEAR(y[i], interp(x[i]), 2E-16);
			}
			for (double v = x.front(); v <= x.back(); v += 0.01)
			{
				EXPECT_NEAR(fun(v), interp(v), 0.1);
			}
		}

		void InterpolatorTest::testOperations(Interpolator i)
		{
			const double a = 0.5;
			Interpolator i_mul(i*a);
			Interpolator i_div(i/a);
			Interpolator i_add(i+a);
			Interpolator i_sub(i-a);
			EXPECT_EQ(i_mul.lowerBound(), i.lowerBound());
			EXPECT_EQ(i_mul.upperBound(), i.upperBound());
			EXPECT_EQ(i_div.lowerBound(), i.lowerBound());
			EXPECT_EQ(i_div.upperBound(), i.upperBound());
			EXPECT_EQ(i_add.lowerBound(), i.lowerBound());
			EXPECT_EQ(i_add.upperBound(), i.upperBound());
			EXPECT_EQ(i_sub.lowerBound(), i.lowerBound());
			EXPECT_EQ(i_sub.upperBound(), i.upperBound());

			double x0, x1;

			if (i.lowerBound() > -std::numeric_limits<double>::infinity())
				x0 = i.lowerBound();
			else
				x0 = -1;
			if (i.upperBound() < std::numeric_limits<double>::infinity())
				x1 = i.upperBound();
			else
				x1 = 1;
			const double dx = (x1 - x0) / 100;

			for (double x = x0; x <= x1;  x += dx)
			{
				const double y = i(x);
				const double tol = 2E-16 * (1 + std::abs(y));
				EXPECT_NEAR(y*a, i_mul(x), tol);
				EXPECT_NEAR(y/a, i_div(x), tol);
				EXPECT_NEAR(y+a, i_add(x), tol);
				EXPECT_NEAR(y-a, i_sub(x), tol);
			}

			i_mul /= a;
			i_div *= a;
			i_add -= a;
			i_sub += a;
			for (double x = x0; x <= x1;  x += dx)
			{
				const double y = i(x);
				const double tol = 2E-16 * (1 + std::abs(y));
				EXPECT_NEAR(y, i_mul(x), tol);
				EXPECT_NEAR(y, i_div(x), tol);
				EXPECT_NEAR(y, i_add(x), tol);
				EXPECT_NEAR(y, i_sub(x), tol);
			}
		}
	}
}
