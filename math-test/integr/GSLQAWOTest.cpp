#include <gtest/gtest.h>
#include "math/integr/GSLQAWO.h"
#include <cmath>

class GSLQAWOTest: public testing::Test
{
};

using namespace rql::math::integr;

TEST_F(GSLQAWOTest,Test)
{
	GSLQAWO wc(2.0, 1, GSL_INTEG_COSINE, 20);
	double result;
	double abserr;
	wc.integrate<double(*)(double)>(sin, 0, 1E-10, 1E-10, 20, result, abserr);
	const double cos1 = cos(1.0);
	double expected = (3*cos1 - cos(3.0)-2)/6.0;
	ASSERT_NEAR(expected, result, 1E-10);
	GSLQAWO ws(2.0, 1, GSL_INTEG_SINE, 20);
	ws.integrate<double(*)(double)>(cos, 0, 1E-10, 1E-10, 20, result, abserr);
	expected = (1-pow(cos1, 3))*2.0/3;
	ASSERT_NEAR(expected, result, 1E-10);
	wc.setType(GSL_INTEG_SINE);
	wc.integrate<double(*)(double)>(cos, 0, 1E-10, 1E-10, 20, result, abserr);
	ASSERT_NEAR(expected, result, 1E-10);
}
