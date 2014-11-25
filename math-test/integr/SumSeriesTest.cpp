#include <gtest/gtest.h>
#include <math/integr/SumSeries.h>
#include <cmath>

class SumSeriesTest: public testing::Test
{
};

struct PowerSeries
{
	double q;
	double operator()(int k)
	{
		return pow(q, k);
	}
};

using namespace rql::math::integr;

TEST_F(SumSeriesTest, Power)
{	
	SumSeries<double> ss;
	PowerSeries power_series;
	power_series.q = 0.99;
	const double expected = 1 / (1 - power_series.q);
	const double actual = ss.sum(power_series, 0., 0, 1);
	ASSERT_NEAR(expected, actual, 2E-11);
	ASSERT_NEAR(expected, 1.0 + ss.sum(power_series, 0, 1, 1), 2E-11);
	ASSERT_NEAR(expected + 1.0, ss.sum(power_series, 1.0, 0, 1), 2E-11);

}

double alternating_series(int n)
{	
	if (n % 2 == 0) {
		return -1.0 / n;
	} else {
		return 1.0 / n;
	}
}

TEST_F(SumSeriesTest, Alternating)
{
	SumSeries<double> ss(1E-4);
	const double expected = log(2.0);
	const double actual = ss.sum(alternating_series, 0.0, 1, 1);
	ASSERT_NEAR(expected, actual, 1E-4);
}
