#include <gtest/gtest.h>
#include "random_vector_levy_symmetric.h"
#include <math/average.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

using namespace rql::math;

class RandomVectorLevySymmetricTest: public testing::Test
{
protected:
	void testMean(double alpha, double tol);
	void testDist(RandomVector& rvg, const std::vector<double>& bukkits, const std::vector<double>& expectedCdf, double tolerance);
	void testIndependenceCovariance(RandomVector& rvg, unsigned int i1, unsigned int i2, double tolerance);
	void testIndependenceKendallTau(RandomVector& rvg, unsigned int i1, unsigned int i2, double tolerance);
};

void RandomVectorLevySymmetricTest::testMean(double alpha, double tol)
{
	const unsigned int dim = 1;
	RandomVectorLevySymmetric rvg(dim,1,alpha);
	std::vector<double> vec(dim);
	std::vector<Average<double> > average(dim);
	double counter = 0;
	const unsigned int nbr_draws = 1000000;
	for (unsigned int n = 0; n < nbr_draws; ++n) {
		rvg.draw(vec);
		++counter;
		for (unsigned int k = 0; k < dim; ++k) {
			average[k].update(vec[k]);
		}
	}
	double norm = 0;
	for (unsigned int k = 0; k < dim; ++k) {
		norm = std::max(std::abs(average[k].value()), norm);
	}
	ASSERT_NEAR(0.0, norm, tol);
}

TEST_F(RandomVectorLevySymmetricTest,MeanFatTails)
{
	testMean(1.5, 3E-2);
}

TEST_F(RandomVectorLevySymmetricTest,MeanGaussian)
{
	testMean(2, 2E-3);
}

void RandomVectorLevySymmetricTest::testDist(RandomVector& rvg, const std::vector<double>& bukkits, const std::vector<double>& expectedCdf, double tolerance)
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

TEST_F(RandomVectorLevySymmetricTest, LevyDist)
{
	const unsigned int dim = 1;
	RandomVectorLevySymmetric rvg(dim,0.1,0.5);	
	std::vector<double> bukkits(6);
	bukkits[0] = 0;
	bukkits[1] = 0.1;
	bukkits[2] = 1;
	bukkits[3] = 10;
	bukkits[4] = 100;
	bukkits[5] = 1000;
	bukkits.push_back(-1);
	bukkits.push_back(-10);
	bukkits.push_back(-100);
	bukkits.push_back(-1000);
	std::vector<double> expectedCdf(6);
	expectedCdf[0] = 0.5;
	expectedCdf[1] = 0.728719687384151;
	expectedCdf[2] = 0.888714609369534;
	expectedCdf[3] = 0.961664099290069;
	expectedCdf[4] = 0.987542439590809;
	expectedCdf[5] = 0.996026457936124;
	expectedCdf.push_back(0.111285390630466);
	expectedCdf.push_back(0.0383359007099306);
	expectedCdf.push_back(0.0124575604091910);
	expectedCdf.push_back(0.00397354206387601);
	testDist(rvg, bukkits, expectedCdf, 0.6);
}

TEST_F(RandomVectorLevySymmetricTest, LorentzDist)
{
	const unsigned int dim = 1;
	RandomVectorLevySymmetric rvg(dim,0.1,1);
	std::vector<double> bukkits(6);
	bukkits[0] = 0;
	bukkits[1] = 0.1;
	bukkits[2] = 0.5;
	bukkits[3] = 1;
	bukkits[4] = 10;
	bukkits[5] = 50;
	bukkits.push_back(-0.1);
	bukkits.push_back(-0.5);
	bukkits.push_back(-10);
	bukkits.push_back(-50);
	std::vector<double> expectedCdf(6);
	expectedCdf[0] = 0.5;
	expectedCdf[1] = 0.75;
	expectedCdf[2] = 0.937167041810999;
	expectedCdf[3] = 0.968274482569446;
	expectedCdf[4] = 0.996817007235092;
	expectedCdf[5] = 0.999363381076457;
	expectedCdf.push_back(0.250000000000000);
	expectedCdf.push_back(0.0628329581890011);
	expectedCdf.push_back(0.00318299276490819);
	expectedCdf.push_back(6.36618923543275e-004);
	testDist(rvg, bukkits, expectedCdf, 0.6);
}

TEST_F(RandomVectorLevySymmetricTest, GaussDist)
{
	const unsigned int dim = 1;
	RandomVectorLevySymmetric rvg(dim,1,2);
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
	testDist(rvg, bukkits, expectedCdf, 0.71);
}

void RandomVectorLevySymmetricTest::testIndependenceCovariance(RandomVector& rv, unsigned int i1, unsigned int i2, double tolerance)
{
	const unsigned int nbrSamples = 1000000;
	Average<double> a1;
	Average<double> a2;
	Average<double> a12;
	std::vector<double> x(rv.dim());
	assert( i1 < rv.dim() );
	assert( i2 < rv.dim() );
	assert( i1 != i2 );
	for (unsigned int i = 0; i < nbrSamples; ++i) {
		rv.draw(x);
		a1.update(x[i1]);
		a2.update(x[i2]);
		a12.update(x[i1]*x[i2]);
	}
	EXPECT_NEAR(0.0, a12.value() - a1.value()*a2.value(), tolerance / sqrt(static_cast<double>(nbrSamples)));
}

void RandomVectorLevySymmetricTest::testIndependenceKendallTau(RandomVector& rv, unsigned int i1, unsigned int i2, double tolerance)
{
	const unsigned int nbrSamples = 20000;	
	std::vector<double> x(rv.dim());
	assert( i1 < rv.dim() );
	assert( i2 < rv.dim() );
	assert( i1 != i2 );
	std::vector<std::pair<double,double> > pairs;
	for (unsigned int i = 0; i < nbrSamples; ++i) {
		rv.draw(x);
		const double x1 = x[i1];
		const double x2 = x[i2];
		pairs.push_back(std::pair<double,double>(x1, x2));
	}
	unsigned int nbrConcordant = 0;
	unsigned int nbrDiscordant = 0;
	for (unsigned int i = 0; i < nbrSamples; ++i) {
		const double xi = pairs[i].first;
		const double yi = pairs[i].second;
		for (unsigned int j = 0; j < i; ++j) {
			const double xj = pairs[j].first;
			const double yj = pairs[j].second;
			if ((xi < xj && yi < yj) || (xi > xj) && (yi > yj))
				++nbrConcordant;
			else if ((xi < xj && yi > yj) || (xi > xj) && (yi < yj))
				++nbrDiscordant;
		}
	}
	const double tau = (nbrConcordant - nbrDiscordant) / (0.5*nbrSamples*(nbrSamples-1.0));
	//std::cout << "Tau: " << tau << std::endl;
	EXPECT_NEAR(0.0, tau, tolerance) << "Kendall Tau";
}

TEST_F(RandomVectorLevySymmetricTest, GaussIndependence)
{
	RandomVectorLevySymmetric rv(5,1,2);
	testIndependenceCovariance(rv, 0, 3, 0.1);
	testIndependenceKendallTau(rv, 0, 3, 0.004);
}

TEST_F(RandomVectorLevySymmetricTest, LevyIndependence)
{
	RandomVectorLevySymmetric rv(5,1,0.5);
	testIndependenceKendallTau(rv, 0, 3, 0.008);
}

TEST_F(RandomVectorLevySymmetricTest, LorentzIndependence)
{
	RandomVectorLevySymmetric rv(5,1,1);
	testIndependenceKendallTau(rv, 0, 3, 0.002);
}
