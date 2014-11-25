#include <gtest/gtest.h>
#include <qsd/correlation_functions_decomposition_factory.h>
#include <qsd/strunz_calculator_reduced.h>
#include <boost/make_shared.hpp>
#include <protein_chain/hamiltonian_factory.h>
#include <protein_chain/correlation_function_lorentzians.h>

class StrunzCalculatorReducedTest: public testing::Test
{
};

TEST_F(StrunzCalculatorReducedTest, ZeroAlpha)
{
	const size_t nbr_sites = 3;
	const size_t nbr_peaks = 2;
	Eigen::VectorXd omega(nbr_peaks);
	Eigen::VectorXd hwhm(nbr_peaks);
	Eigen::VectorXd height(nbr_peaks);
	omega[0] = 0.5;
	omega[1] = 1;
	hwhm[0] = 0.1;
	hwhm[1] = 0.25;
	height[0] = 0.0;
	height[1] = 0.0;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(omega, hwhm, height);
	std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> > alphas(nbr_sites, alpha);
	Eigen::MatrixXcd hamiltonian(nbr_sites, nbr_sites);
	HamiltonianFactory::buildChain(std::vector<double>(nbr_sites, 0.0), -1, 0.0, hamiltonian);
	StrunzCalculatorReduced calculator(CorrelationFunctionsDecompositionFactory::decomposeLorentzians(alphas), hamiltonian);
	StrunzCalculatorReduced::Data state(calculator.data());
	StrunzCalculatorReduced::Data f(calculator.data());
	Eigen::MatrixXcd m(nbr_sites, nbr_sites);
	calculator.calculate_effective_hamiltonian_times_minus_i(state, m);
	ASSERT_EQ(hamiltonian*std::complex<double>(0,-1), m);
	calculator.calculate_time_derivative(state, m, f);
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			EXPECT_EQ(Eigen::MatrixXcd::Zero(nbr_sites, nbr_sites), f(m, k));
	calculator.calculate_effective_hamiltonian_times_minus_i_and_time_derivative(state, m, f);
	ASSERT_EQ(hamiltonian*std::complex<double>(0,-1), m);
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			ASSERT_EQ(Eigen::MatrixXcd::Zero(nbr_sites, nbr_sites), f(m, k));
}

TEST_F(StrunzCalculatorReducedTest, Data)
{
	const size_t nbr_sites = 3;
	const size_t nbr_peaks = 2;
	Eigen::VectorXd omega(nbr_peaks);
	Eigen::VectorXd hwhm(nbr_peaks);
	Eigen::VectorXd height(nbr_peaks);
	omega[0] = 0.5;
	omega[1] = 1;
	hwhm[0] = 0.1;
	hwhm[1] = 0.25;
	height[0] = 0.0;
	height[1] = 0.0;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(omega, hwhm, height);
	std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> > alphas(nbr_sites, alpha);
	Eigen::MatrixXcd hamiltonian(nbr_sites, nbr_sites);
	HamiltonianFactory::buildChain(std::vector<double>(nbr_sites, 0.0), -1, 0.0, hamiltonian);
	StrunzCalculatorReduced calculator(CorrelationFunctionsDecompositionFactory::decomposeLorentzians(alphas), hamiltonian);
	StrunzCalculatorReduced::Data d1(calculator.data());
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			EXPECT_EQ(Eigen::MatrixXcd::Zero(nbr_sites, nbr_sites), d1(m, k));
	StrunzCalculatorReduced::Data d2(calculator.data());
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k) {
			d1(m,k) = Eigen::MatrixXcd::Random(nbr_sites, nbr_sites);
			d2(m,k) = Eigen::MatrixXcd::Random(nbr_sites, nbr_sites);
		}
	StrunzCalculatorReduced::Data d3(calculator.data());
	d3.copy_data_from(d1);
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			EXPECT_EQ(d1(m, k), d3(m,k));
	d3 *= 2;
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			ASSERT_EQ(2*d1(m, k), d3(m,k));
	d3 /= 2;
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			ASSERT_EQ(d1(m, k), d3(m,k)) << m << " " << k;
	d3 += d2;
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			ASSERT_NEAR( (d1(m, k) + d2(m,k) - d3(m,k)).norm(), 0.0, 1E-14) << m << " " << k;
	d3.copy_data_from(d1);
	d3.add(d2, -0.21);
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			ASSERT_NEAR( (d1(m, k) - 0.21*d2(m,k) - d3(m,k)).norm(), 0.0, 1E-14) << m << " " << k;
	d3.clear();
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k)
			ASSERT_EQ(Eigen::MatrixXcd::Zero(nbr_sites, nbr_sites), d3(m,k)) << m << " " << k;
}

TEST_F(StrunzCalculatorReducedTest, Monomers)
{
	const size_t nbr_sites = 3;
	const size_t nbr_peaks = 2;
	Eigen::VectorXd omega(nbr_peaks);
	Eigen::VectorXd hwhm(nbr_peaks);
	Eigen::VectorXd height(nbr_peaks);
	omega[0] = 0.5;
	omega[1] = 1;
	hwhm[0] = 0.1;
	hwhm[1] = 0.25;
	height[0] = 0.9;
	height[1] = 1.8;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(omega, hwhm, height);
	std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> > alphas(nbr_sites, alpha);
	Eigen::MatrixXcd hamiltonian(nbr_sites, nbr_sites);
	hamiltonian.setZero();
	hamiltonian(0,0)=-1;
	hamiltonian(2,2)=1;
	StrunzCalculatorReduced calculator(CorrelationFunctionsDecompositionFactory::decomposeLorentzians(alphas), hamiltonian);
	StrunzCalculatorReduced::Data state(calculator.data());
	StrunzCalculatorReduced::Data f(calculator.data());
	const double dt = 0.1;
	Eigen::MatrixXcd m(nbr_sites, nbr_sites);
	for (size_t i = 0; i < 20; ++i) {
		calculator.calculate_effective_hamiltonian_times_minus_i(state, m);
		calculator.calculate_time_derivative(state, m, f);
		state.add(f, dt);
		calculator.calculate_effective_hamiltonian_times_minus_i_and_time_derivative(state, m, f);
		state.add(f, dt);
	}
	for (size_t m = 0; m < nbr_sites; ++m)
		for (size_t k = 0; k < nbr_peaks; ++k) {
			Eigen::MatrixXcd tmp = state(m, k);
			tmp(m,m) = 0;
			ASSERT_EQ(0.0, tmp.norm()) << m << " " << k;
			tmp = f(m,k);
			tmp(m,m) = 0;
			ASSERT_EQ(0.0, tmp.norm()) << m << " " << k;
			ASSERT_EQ(state(0, k)(0,0), state(m, k)(m,m)) << m << " " << k;
			ASSERT_EQ(f(0, k)(0,0), f(m, k)(m,m)) << m << " " << k;
		}
	for (size_t m1 = 0; m1 < nbr_sites; ++m1)
		m(m1,m1) = 0;
	ASSERT_EQ(0.0, m.norm());
}

TEST_F(StrunzCalculatorReducedTest, TraceNormalization)
{
	const size_t nbr_sites = 3;
	const size_t nbr_peaks = 2;
	Eigen::VectorXd omega(nbr_peaks);
	Eigen::VectorXd hwhm(nbr_peaks);
	Eigen::VectorXd height(nbr_peaks);
	omega[0] = 0.5;
	omega[1] = 1;
	hwhm[0] = 0.1;
	hwhm[1] = 0.25;
	height[0] = 0.9;
	height[1] = 1.8;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(omega, hwhm, height);
	std::vector<boost::shared_ptr<const CorrelationFunctionLorentzians> > alphas(nbr_sites, alpha);
	Eigen::MatrixXcd hamiltonian(nbr_sites, nbr_sites);
	hamiltonian.setZero();
	hamiltonian(0,0)=-1;
	hamiltonian(2,2)=1;
	StrunzCalculatorReduced calculator(CorrelationFunctionsDecompositionFactory::decomposeLorentzians(alphas), hamiltonian);
	StrunzCalculatorReduced::Data state(calculator.data());
	for (size_t m = 0; m < nbr_sites; ++m) {
		for (size_t k = 0; k < nbr_peaks; ++k) {
			state(m, k) = Eigen::MatrixXcd::Random(nbr_sites, nbr_sites);			
		}
	}
	const double time = 10.1;
	calculator.normalize_traces(state, time);
	ASSERT_NEAR(0.0, calculator.total_absolute_trace_deviation_from_theoretical(state, time), 1E-10);
}
