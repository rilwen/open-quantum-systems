#include <gtest/gtest.h>
#include <Eigen/Core>
#include <cmath>
#include <array>
#include "math/interp/AkimaApproximation.h"

class AkimaApproximationTest: public testing::Test
{
};

TEST_F(AkimaApproximationTest, Test)
{
	const size_t n = 10;
	std::vector<double> x(n);
	Eigen::VectorXd y(x.size());
	std::array<double,n> dy;
	std::vector<double> exact_dy(n);
	for (unsigned int i = 0; i < n; ++i)
	{
		x[i] = i;
		y[i] = cos(i*0.5);
		exact_dy[i] = -0.5*sin(0.5*i);
	}
	rql::interp::AkimaApproximation<double>::calculate(x, y, dy);
	for (unsigned int i = 0; i < n; ++i)
		EXPECT_NEAR(exact_dy[i], dy[i], 0.05);
	for (size_t i = 0; i < n; ++i) {
		y[i] = 2;
	}
	rql::interp::AkimaApproximation<double>::calculate(x, y, dy);
	for (size_t i = 0; i < n; ++i)
		EXPECT_NEAR(0, dy[i], 1E-15);
}
