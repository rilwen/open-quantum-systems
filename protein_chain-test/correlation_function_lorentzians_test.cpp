#include <gtest/gtest.h>
#include <protein_chain/correlation_function_lorentzians.h>

class CorrelationFunctionLorentziansTest: public testing::Test
{
};

TEST_F(CorrelationFunctionLorentziansTest, Test)
{
	const double omega = 0.5;
	const double gamma = 1;
	const double height = 1;
	CorrelationFunctionLorentzians alpha(omega, gamma, height);
	/*const double d = 1E-9;
	const double t = 0.5;
	const double s = 0.22;*/
	/*const size_t k_t = static_cast<size_t>(t/d);
	const size_t k_s = static_cast<size_t>(s/d);*/
	//ASSERT_NEAR(alpha(t - s).real()*d*d, alpha.covariance(k_t - k_s, d).real(), 3*d*d*d);
	//ASSERT_NEAR(alpha(t - s).imag()*d*d, alpha.covariance(k_t - k_s, d).imag(), 3*d*d*d);
	const double tau = 0.34;
	const std::complex<double> expected = 3.1415926535897931*exp(std::complex<double>(-gamma,-omega)*tau);
	ASSERT_NEAR(std::abs(expected - alpha(tau)), 0.0, 1E-10);
	ASSERT_NEAR(3.1415926535897931, alpha(0).real(), 1E-10);
	ASSERT_NEAR(0.0, alpha(0).imag(), 1E-10);
	ASSERT_EQ(alpha(-1).real(), alpha(1).real());
	ASSERT_EQ(alpha(-1).imag(), -alpha(1).imag());

	alpha = CorrelationFunctionLorentzians(0.0, 1, 1);
	ASSERT_TRUE(alpha(-1) == alpha(1));

	const double expected_scale = 3.1415926535897931;
	ASSERT_NEAR(expected_scale, CorrelationFunctionLorentzians::scale(1, 1), 1E-10);
	ASSERT_NEAR(0.456, CorrelationFunctionLorentzians::height(CorrelationFunctionLorentzians::scale(0.456, 0.91), 0.91), 1E-10);
}

TEST_F(CorrelationFunctionLorentziansTest, DISABLED_Temperature)
{
	const size_t pfd_order = 20;

	CorrelationFunctionLorentzians J(1.0, 0.2, 1.0, 0.000001, pfd_order);
	const std::complex<double> v0_vlow = J(0);
	ASSERT_NEAR(v0_vlow.imag(), 0.0, 1E-6);
	ASSERT_TRUE(v0_vlow.real() > 0.0);
	std::complex<double> v1 = J(0.5);
	std::complex<double> v2 = J(-0.5);
	ASSERT_NEAR((v1+v2).imag(), 0.0, 0);
	ASSERT_NEAR(v1.real(), v2.real(), 0);

	J = CorrelationFunctionLorentzians(1.0, 0.2, 1.0, 0.1, pfd_order);
	const std::complex<double> v0_low = J(0);
	ASSERT_NEAR(v0_low.imag(), 0.0, 1E-5);
	ASSERT_TRUE(v0_low.real() > 0.0);
	v1 = J(0.5);
	v2 = J(-0.5);
	ASSERT_NEAR((v1+v2).imag(), 0.0, 0);
	ASSERT_NEAR(v1.real(), v2.real(), 0);

	J = CorrelationFunctionLorentzians(1.0, 0.2, 1.0, 1, pfd_order);
	const std::complex<double> v0_mid = J(0);
	ASSERT_NEAR(v0_mid.imag(), 0.0, 1E-5);
	ASSERT_TRUE(v0_mid.real() > 0.0);
	ASSERT_LT(v0_mid.real(), v0_low.real());
	v1 = J(0.5);
	v2 = J(-0.5);
	ASSERT_NEAR((v1+v2).imag(), 0.0, 0);
	ASSERT_NEAR(v1.real(), v2.real(), 0);

	J = CorrelationFunctionLorentzians(1.0, 0.2, 1.0, 10, pfd_order);
	const std::complex<double> v0_high = J(0);
	ASSERT_NEAR(v0_high.imag(), 0.0, 1E-5);
	ASSERT_TRUE(v0_high.real() > 0.0);
	ASSERT_TRUE(v0_high.real() < v0_mid.real());
	v1 = J(0.5);
	v2 = J(-0.5);
	ASSERT_NEAR((v1+v2).imag(), 0.0, 0);
	ASSERT_NEAR(v1.real(), v2.real(), 0);
}

