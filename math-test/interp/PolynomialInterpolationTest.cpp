#include <gtest/gtest.h>
#include "math/interp/PolynomialInterpolation.h"
#include <boost/array.hpp>
#include <array>

using namespace rql::interp;

class PolynomialInterpolationTest: public testing::Test
{
};

TEST_F(PolynomialInterpolationTest,Linear)
{
	const double x = 1;
	// f(x) = 0.5 + x
	const boost::array<double,1> y0 = {0.5};
	const boost::array<double,1> yx = {1.5};
	boost::array<double,2> a;
	PolynomialInterpolation::interpolate(x, y0, yx, a);
	EXPECT_NEAR(y0[0], a[0], 1E-16);
	EXPECT_NEAR(yx[0]-y0[0], a[1], 1E-16);
}

TEST_F(PolynomialInterpolationTest,Cubic)
{
	const double x = 1;
	// f(x) = 0.1*x^3 + x^2 - 2*x + 0.5
	const std::array<double,2> y0 = {0.5,-2.0};
	const std::array<double,2> yx = {-0.4,0.3};
	std::array<double,4> a;
	PolynomialInterpolation::interpolate(x, y0, yx, a);
	EXPECT_NEAR(0.5, a[0], 1E-16);
	EXPECT_NEAR(-2, a[1], 1E-16);
	EXPECT_NEAR(1, a[2], 1E-16);
	EXPECT_NEAR(0.1, a[3], 1E-16);
}
