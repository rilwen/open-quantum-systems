#include <gtest/gtest.h>
#include "random_vector_gaussian.h"
#include "random_vector_levy_symmetric.h"
#include <math/stats/covariance.h>
#include <vector>
#include <iostream>
#include <cassert>

using namespace rql::math::stats;

class RandomVectorGaussianTest: public testing::Test
{
protected:
	void testDist(RandomVector& rvg, const std::vector<double>& bukkits, const std::vector<double>& expectedCdf, double tolerance);
};

TEST_F(RandomVectorGaussianTest,Mean)
{
	const unsigned int dim = 30;
	RandomVectorGaussian rvg(dim);
//	RandomVectorLevySymmetric rvg(dim,1,1);	
	std::vector<double> vec(dim);
	std::vector<double> average(dim, 0.0);
	double counter = 0;
	const unsigned int nbr_draws = 20000;
	for (unsigned int n = 0; n < nbr_draws; ++n) {
		rvg.draw(vec);
		++counter;
		for (unsigned int k = 0; k < dim; ++k) {
			const double old = average[k];
			average[k] = old + (vec[k] - old) / counter;
		}
	}
	double norm = 0;
	for (unsigned int k = 0; k < dim; ++k) {
		norm = std::max(std::abs(average[k]), norm);		
	}
	ASSERT_NEAR(0.0, norm, 1e-8);
}

TEST_F(RandomVectorGaussianTest,Covariance)
{
	const unsigned int dim = 10;
	RandomVectorGaussian rvg(dim);
	Covariance covcalc(dim);
	std::vector<double> vec(dim);
	const unsigned int nbr_draws = 60000;
	for (unsigned int n = 0; n < nbr_draws; ++n) {
		rvg.draw(vec);
		covcalc.update(vec.begin(), vec.end());
	}

	Covariance::matrix_type cov(covcalc.covariance());
	double norm = 0;
	double norm_diag = 0;
	for (unsigned int k = 0; k < dim; ++k) {
		norm_diag = std::max(std::abs(cov(k,k)-1), norm_diag);
		for (unsigned int l = 0; l < k; ++l) {
			norm = std::max(std::abs(cov(k,l)), norm);		
		}
	}
	//std::cout << norm_diag << " " << norm << std::endl;
	EXPECT_NEAR(0.0, norm_diag, 0.015);
	EXPECT_NEAR(0.0, norm, 0.015);
}

void RandomVectorGaussianTest::testDist(RandomVector& rvg, const std::vector<double>& bukkits, const std::vector<double>& expectedCdf, double tolerance)
{
	assert(bukkits.size() == expectedCdf.size());
	const unsigned int nbrSamples = 1000000;
	std::vector<unsigned int> counters(bukkits.size(), 0);
	std::vector<double> x(1);
	for (unsigned int i = 0; i < nbrSamples; ++i) {
		rvg.draw(x);
		for (unsigned int j = 0; j < bukkits.size(); ++j) {
			if (x[0] < bukkits[j]) {
				++counters[j];
			}
		}
	}
	for (unsigned int j = 0; j < bukkits.size(); ++j) {
		const double cdf = static_cast<double>(counters[j]) / nbrSamples;
		EXPECT_NEAR(expectedCdf[j], cdf, tolerance/sqrt(static_cast<double>(nbrSamples))) << "bukkits[j] == " << bukkits[j];
	}
}

TEST_F(RandomVectorGaussianTest,Distribution)
{
	const unsigned int dim = 1;
	RandomVectorGaussian rvg(dim);	
	std::vector<double> bukkits;
	bukkits.push_back(-6);
	bukkits.push_back(-5);
	bukkits.push_back(-4);
	bukkits.push_back(-3);
	bukkits.push_back(-2);
	bukkits.push_back(-1);
	bukkits.push_back(0);
	bukkits.push_back(1);
	bukkits.push_back(2);
	bukkits.push_back(3);
	bukkits.push_back(4);
	bukkits.push_back(5);
	bukkits.push_back(6);
	bukkits.push_back(-0.5);
	bukkits.push_back(0.5);
	std::vector<double> expectedCdf;
	expectedCdf.push_back(9.86587644913328e-010);
	expectedCdf.push_back(2.86651571868024e-007);
	expectedCdf.push_back(3.16712418331200e-005);
	expectedCdf.push_back(0.00134989803163010);
	expectedCdf.push_back(0.0227501319481792);
	expectedCdf.push_back(0.158655253931457);
	expectedCdf.push_back(0.5);
	expectedCdf.push_back(0.841344746068543);
	expectedCdf.push_back(0.977249868051821);
	expectedCdf.push_back(0.998650101968370);
	expectedCdf.push_back(0.999968328758167);
	expectedCdf.push_back(0.999999713348428);
	expectedCdf.push_back(0.999999999013412);
	expectedCdf.push_back(0.308537538725987);
	expectedCdf.push_back(0.691462461274013);
	testDist(rvg, bukkits, expectedCdf, 0.61);
}