#include <gtest/gtest.h>
#include "math/integr/FourierIntegral.h"
#include <cmath>

class FourierIntegralTest: public testing::Test
{
};

using namespace rql::integr;

static std::complex<double> function(double x)
{
	return std::exp(std::complex<double>(0, -x));
}

TEST_F(FourierIntegralTest, Test)
{
	FourierIntegral fi(2, 0, 1, 100);
	const std::complex<double> actual = fi.integrate(function);
	const std::complex<double> expected = (std::exp(std::complex<double>(0, 1)) - std::complex<double>(1,0)) / std::complex<double>(0,1);
	ASSERT_NEAR(expected.real(), actual.real(), 1E-8) << "real";
	ASSERT_NEAR(expected.imag(), actual.imag(), 1E-7) << "imag";
}
