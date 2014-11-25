#include <gtest/gtest.h>
#include "math/integr/GSLFourierIntegral.h"

class GSLFourierIntegralTest: public testing::Test
{
};

using namespace rql::math::integr;

static std::complex<double> function(double x)
{
	return std::exp(std::complex<double>(0, -x));
}

TEST_F(GSLFourierIntegralTest, Test)
{
	const size_t n = 20;
	GSLFourierIntegral fi(2, 1, n);
	std::complex<double> actual;
	std::complex<double> abserr;
	fi.integrate(function, 0, 1E-10, 1E-10, n, actual, abserr);
	const std::complex<double> expected = (std::exp(std::complex<double>(0, 1)) - std::complex<double>(1,0)) / std::complex<double>(0,1);
	ASSERT_NEAR(expected.real(), actual.real(), 1E-10) << "real";
	ASSERT_NEAR(expected.imag(), actual.imag(), 1E-10) << "imag";
	ASSERT_TRUE(std::abs(abserr) < 1E-10);
}

