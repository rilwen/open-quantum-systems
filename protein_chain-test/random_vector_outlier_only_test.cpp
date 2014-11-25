#include <gtest/gtest.h>
#include <cmath>
#include "random_vector_outlier_only.h"

class RandomVectorOutlierOnlyTest: public testing::Test
{
};

TEST_F(RandomVectorOutlierOnlyTest, Distr)
{
	const unsigned int dim = 4;
	const double outlierProb = 0.1;
	const double outlierSize = 4;	
	std::vector<unsigned int> cntrs_neg(dim, 0u);
	std::vector<unsigned int> cntrs_pos(dim, 0u);	
	std::vector<double> epsilons(dim);
	RandomVectorOutlierOnly rvo(dim, outlierProb, outlierSize);
	ASSERT_EQ(dim, rvo.dim());
	const unsigned int iters = 90000;
	for (unsigned int i = 0; i < iters; ++i) {
		rvo.draw(epsilons);
		for (unsigned int j = 0; j < dim; ++j) {
			const double eps = epsilons[j];
			ASSERT_TRUE(eps == 0 || std::abs(eps) == outlierSize) << i << " " << j;
			if (eps == -outlierSize) {
				++cntrs_neg[j];
			} else if (eps == outlierSize) {
				++cntrs_pos[j];
			}
		}
	}
	const double tol = 0.5/sqrt(static_cast<double>(iters));
	for (unsigned int j = 0; j < dim; ++j) {
		const double ratio_neg = static_cast<double>(cntrs_neg[j]) / iters;
		EXPECT_NEAR(outlierProb/2, ratio_neg, tol) << "ratio_neg: " << j;
		const double ratio_pos = static_cast<double>(cntrs_pos[j]) / iters;
		EXPECT_NEAR(outlierProb/2, ratio_pos, tol) << "ratio_pos: " << j;
	}
}

