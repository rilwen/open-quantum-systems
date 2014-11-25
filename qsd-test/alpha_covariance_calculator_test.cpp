#include <gtest/gtest.h>
#include <qsd/alpha_covariance_calculator.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <boost/make_shared.hpp>

class AlphaCovarianceCalculatorTest: public testing::Test
{
};

TEST_F(AlphaCovarianceCalculatorTest,test)
{
	const size_t nbr_steps = 4;
	const size_t dim = 2;
	
	const double dt = 0.5;
	std::vector<boost::shared_ptr<const CorrelationFunction> > alpha(dim, boost::make_shared<CorrelationFunctionLorentzians>(1, 0.25, 0.5));
	alpha[1] = boost::make_shared<CorrelationFunctionLorentzians>(0.1, 0.5, 0.25);

	AlphaCovarianceCalculator acc(alpha, nbr_steps, dt);
	ASSERT_EQ(dim, acc.dim());
	ASSERT_EQ(nbr_steps, acc.nbr_time_pts());
	ASSERT_EQ(dt, acc.delta());
	for (size_t m = 0; m < dim; ++m) {
		SingleAlphaCovarianceCalculator sacc(acc.single_alpha_calculator(m));
		for (size_t k = 0; k < nbr_steps; ++k)
			for (size_t l = 0; l < nbr_steps; ++l) {
				const std::complex<double> expected = (*alpha[m])(static_cast<int>(k-l)*dt);
				ASSERT_EQ(expected, acc.covariance(m, k, l)) << "(m, k, l) == (" << m << ", " << k << ", " << l << ")";
				ASSERT_EQ(expected, sacc.covariance(k, l)) << "(m, k, l) == (" << m << ", " << k << ", " << l << ")";
				ASSERT_EQ(acc.covariance_RI(m, k, l), sacc.covariance_RI(k, l)) << "(m, k, l) == (" << m << ", " << k << ", " << l << ")";
				ASSERT_EQ(acc.covariance_RR_II(m, k, l), sacc.covariance_RR_II(k, l)) << "(m, k, l) == (" << m << ", " << k << ", " << l << ")";
				ASSERT_EQ(0.5*expected.real(), acc.covariance_RR_II(m, k, l)) << "(m, k, l) == (" << m << ", " << k << ", " << l << ")";
				ASSERT_EQ(0.5*expected.imag(), acc.covariance_RI(m, k, l)) << "(m, k, l) == (" << m << ", " << k << ", " << l << ")";
			}
	}
}

TEST_F(AlphaCovarianceCalculatorTest,test_common)
{
	const size_t nbr_steps = 4;
	const size_t dim = 2;
	
	const double dt = 0.5;
	const boost::shared_ptr<const CorrelationFunction> alpha = boost::make_shared<CorrelationFunctionLorentzians>(1, 0.25, 0.5);

	AlphaCovarianceCalculator acc(alpha, dim, nbr_steps, dt);
	ASSERT_EQ(dim, acc.dim());
	ASSERT_EQ(nbr_steps, acc.nbr_time_pts());
	ASSERT_EQ(dt, acc.delta());
	for (size_t m = 0; m < dim; ++m) {
		SingleAlphaCovarianceCalculator sacc(acc.single_alpha_calculator(m));
		for (size_t k = 0; k < nbr_steps; ++k)
			for (size_t l = 0; l < nbr_steps; ++l) {
				const std::complex<double> expected = (*alpha)(static_cast<int>(k-l)*dt);
				ASSERT_EQ(expected, acc.covariance(m, k, l));
				ASSERT_EQ(expected, sacc.covariance(k, l));
				ASSERT_EQ(acc.covariance_RI(m, k, l), sacc.covariance_RI(k, l));
				ASSERT_EQ(acc.covariance_RR_II(m, k, l), sacc.covariance_RR_II(k, l));
				ASSERT_EQ(0.5*expected.real(), acc.covariance_RR_II(m, k, l));
				ASSERT_EQ(0.5*expected.imag(), acc.covariance_RI(m, k, l));
			}
	}
}

