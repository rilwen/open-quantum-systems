#include <gtest/gtest.h>
#include "double_exponential_quadrature.h"
#include <cmath>

class DoubleExponentialQuadratureTest: public testing::Test
{
};

TEST_F(DoubleExponentialQuadratureTest,Cosine)
{
	std::vector<double> samplingpnts;
	std::vector<double> weights;

	DoubleExponentialQuadrature::build(1,2,25, samplingpnts, weights);

	double sum = 0;
	for (unsigned int k = 0; k < samplingpnts.size(); ++k) {
		sum += weights[k]*cos(samplingpnts[k]);
	}
	ASSERT_NEAR(sin(2.0)-sin(1.0), sum, 1e-10);	
}

TEST_F(DoubleExponentialQuadratureTest,NonAnalytic)
{
	std::vector<double> samplingpnts;
	std::vector<double> weights;

	DoubleExponentialQuadrature::build(0,3,50, samplingpnts, weights);

	double sum = 0;
	for (unsigned int k = 0; k < samplingpnts.size(); ++k) {
		if (samplingpnts[k] < 2) sum += weights[k]*1;
		else sum += weights[k]*(samplingpnts[k] - 1);
	}
	ASSERT_NEAR(3.5, sum, 5e-3);
}

TEST_F(DoubleExponentialQuadratureTest,Discontinuous)
{
	std::vector<double> samplingpnts;
	std::vector<double> weights;

	DoubleExponentialQuadrature::build(0,3,500, samplingpnts, weights);

	double sum = 0;
	for (unsigned int k = 0; k < samplingpnts.size(); ++k) {
		if (samplingpnts[k] < 2) sum += weights[k]*1;
		else sum += weights[k]*(samplingpnts[k] - 2);
	}
	ASSERT_NEAR(2.5, sum, 5e-3);
}