#include <gtest/gtest.h>
#include <qsd/strunz_calculator_grid_explicit.h>
#include "gaussian_correlation_function.h"
#include <protein_chain/hamiltonian_factory.h>

class StrunzCalculatorGridExplicitTest: public testing::Test
{
};

TEST_F(StrunzCalculatorGridExplicitTest,Test)
{
	boost::shared_ptr<const CorrelationFunction> singleAlpha = boost::shared_ptr<CorrelationFunction>(new Gaussian);
	const size_t nbrSites = 2;
	std::vector<boost::shared_ptr<const CorrelationFunction> > alpha(nbrSites, singleAlpha);
	Eigen::MatrixXcd hamiltonian(nbrSites, nbrSites);
	HamiltonianFactory::buildChain(std::vector<double>(nbrSites, 0.0), -1, 0, hamiltonian);
	const double dt = 0.2;
	const size_t nbrTimePts = 2;
	StrunzCalculatorGridExplicit calc(alpha, hamiltonian, nbrTimePts, dt);
	ASSERT_EQ(dt, calc.dt());
	ASSERT_EQ(nbrSites * nbrTimePts, calc.nbrOperators());
	StrunzCalculatorGridExplicit::Workspace wksp(calc.workspace());
	calc.step(wksp);
	Eigen::MatrixXcd Heff = calc.effectiveHamiltonianTimesMinusI(wksp) * std::complex<double>(0, 1);
	for (size_t m = 0; m < 2; ++m)
		for (size_t n = 0; n < 2; ++n) {
			EXPECT_EQ(m == n ? 0.0 : -1.0, Heff(m,n).real());
			EXPECT_EQ(0.0, Heff(m,n).imag());
		}
	calc.step(wksp);
	Heff = calc.effectiveHamiltonianTimesMinusI(wksp) * std::complex<double>(0, 1);
	for (size_t m = 0; m < 2; ++m)
		for (size_t n = 0; n < 2; ++n) {
			EXPECT_EQ(m == n ? 0.0 : -1.0192157887830464, Heff(m,n).real());
			EXPECT_NEAR(m == n ? -0.196079 : 0.0, Heff(m,n).imag(), 1E-6);
		}
}
