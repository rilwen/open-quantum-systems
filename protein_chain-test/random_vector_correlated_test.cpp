#include <gtest/gtest.h>
#include "random_vector_correlated.h"
#include "random_vector_gaussian.h"
#include <math/stats/covariance.h>
#include <boost/make_shared.hpp>

using namespace rql::math::stats;

class RandomVectorCorrelatedTest: public testing::Test
{
};

TEST_F(RandomVectorCorrelatedTest,Covariance)
{
	const unsigned int dim = 3;
	boost::shared_ptr<RandomVector> iid = boost::make_shared<RandomVectorGaussian>(dim, 1.0);
	Eigen::MatrixXd cov(dim, dim);
	cov(0,0) = 1;
	cov(1,1) = 0.5;
	cov(2,2) = 1.5;
	cov(0,1) = -0.1;
	cov(1,0) = cov(0,1);
	cov(0,2) = 0.21;
	cov(2,0) = cov(0,2);
	cov(1,2) = 0.5;
	cov(2,1) = cov(1,2);
	RandomVectorCorrelated rvg(iid, cov);
	Covariance covcalc(dim);
	std::vector<double> vec(dim);
	const unsigned int nbr_draws = 120000;
	for (unsigned int n = 0; n < nbr_draws; ++n) {
		rvg.draw(vec);
		covcalc.update(vec.begin(), vec.end());
	}

	Covariance::matrix_type measured_covariance(covcalc.covariance());
	//std::cout << norm_diag << " " << norm << std::endl;
	EXPECT_NEAR(0.0, (cov - measured_covariance).norm(), 0.01);
}
