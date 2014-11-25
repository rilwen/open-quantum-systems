#include <gtest/gtest.h>
#include "math/integr/WeightedIntegral.h"
#include "math/integr/IntegrationDensities.h"
#include <stdexcept>
#include <cmath>

class WeightedIntegralTest: public testing::Test
{
};

TEST_F(WeightedIntegralTest,Test)
{
	rql::integr::WeightedIntegral wi(rql::integr::WeightedIntegral::gen_xs(0, 1, 100), rql::integr::IntegrationDensities::Constant(1.0));
	const double integral = wi.integrate<double(*)(double)>(cos);
	const double sin1 = sin(1.0);
	ASSERT_NEAR(sin1, integral, 1E-5);
	rql::integr::WeightedIntegral ws(rql::integr::WeightedIntegral::gen_xs(0, 1, 100), rql::integr::IntegrationDensities::Sine(2.0));
	const double integral2 = ws.integrate<double(*)(double)>(cos);
	const double cos1 = cos(1.0);
	const double expected2 = (1-pow(cos1, 3))*2.0/3;
	ASSERT_NEAR(expected2, integral2, 1E-8);
	rql::integr::WeightedIntegral wc(rql::integr::WeightedIntegral::gen_xs(0, 1, 100), rql::integr::IntegrationDensities::Cosine(2.0));
	const double integral3 = wc.integrate<double(*)(double)>(sin);
	const double expected3 = (3*cos1 - cos(3.0)-2)/6.0;
	ASSERT_NEAR(expected3, integral3, 1E-8);
}

double constant_one(double x)
{
	return 1;
}

double linear(double x)
{
	return x;
}

TEST_F(WeightedIntegralTest,Sine)
{
	rql::integr::WeightedIntegral wu(rql::integr::WeightedIntegral::gen_xs(0, 1, 100), rql::integr::IntegrationDensities::Constant(1.0));
	const double integral_u = wu.integrate<double(*)(double)>(sin);
	ASSERT_NEAR(1 - cos(1.0), integral_u, 1E-5) << "unity";
	rql::integr::WeightedIntegral ws(rql::integr::WeightedIntegral::gen_xs(0, 1, 100), rql::integr::IntegrationDensities::Sine(1.0));
	const double integral_s = ws.integrate(constant_one);
	ASSERT_NEAR(1 - cos(1.0), integral_s, 1E-10) << "sine_weight";
	const double integral_linear = ws.integrate(linear);
	ASSERT_NEAR(sin(1.0) - cos(1.0), integral_linear, 1E-10);
}
