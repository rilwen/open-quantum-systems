#include <gtest/gtest.h>
#include <protein_chain/correlation_function_temperature_calculator.h>
#include <protein_chain/complex_analytic_spectral_density.h>
#include <protein_chain/partial_function_decomposition_bose.h>
#include <math/MathUtils.h>
#include <cmath>
#include <stdexcept>

class CorrelationFunctionTemperatureCalculatorTest: public testing::Test
{
protected:
	CorrelationFunctionTemperatureCalculatorTest();
	void test_simple_lorentzian(double gamma, double kBT, bool high_order, bool positive_tau);

	PartialFunctionDecompositionBose m_pfd_lo;
	PartialFunctionDecompositionBose m_pfd_hi;
};

CorrelationFunctionTemperatureCalculatorTest::CorrelationFunctionTemperatureCalculatorTest()
	: m_pfd_lo(10), m_pfd_hi(40)
{
}


class SimpleLorentzian: public ComplexAnalyticSpectralDensity
{
public:
  SimpleLorentzian(double gamma)
  {
    m_gamma = gamma;
  }
  double spectral_density(double w) const
  {
    return spectral_density(std::complex<double>(w)).real();
  }
  virtual std::complex<double> spectral_density(const std::complex<double>& w) const
  {
    return m_gamma / (w*w + m_gamma*m_gamma) / rql::math::PI;
  }
  virtual double spectral_density_derivative(double w) const
  {
    return - 2 * w * m_gamma / pow(w*w + m_gamma*m_gamma, 2) / rql::math::PI;
  }
  virtual size_t nbrPoles() const { return 2; }
  std::complex<double> pole(size_t idx) const 
  { 
    if (idx == 0) {
      return std::complex<double>(0.0, m_gamma);
    } else {
      return std::complex<double>(0.0, -m_gamma);
    }
  }
  std::complex<double> residue(size_t idx) const 
  { 
    if (idx == 0) {
      return std::complex<double>(0.0, -0.5/rql::math::PI);
    } else {
      return std::complex<double>(0.0, 0.5/rql::math::PI);
    }
  }
private:
  double m_gamma;
};

void CorrelationFunctionTemperatureCalculatorTest::test_simple_lorentzian(double gamma, double kBT, bool high_order, bool positive_tau)
{	
	SimpleLorentzian J(gamma);
	const PartialFunctionDecompositionBose& pfd = high_order ? m_pfd_hi : m_pfd_lo;
	//std::cout << "Test: " << gamma << " " << kBT << " " << pfd.order() << " " << positive_tau << "\n";
	CorrelationFunctionTemperatureCalculator calc(J, pfd, kBT);
	for (int i = -10; i <= 10; ++i) {
		const double tau = i*0.3 / kBT;
		const std::complex<double> a0 = calc.evaluate(tau);
		//std::cout << tau << "\t" << a0 << "\n";
		const std::complex<double> a1 = calc.evaluate(-tau);
		EXPECT_NEAR(a0.real(), a1.real(), 1E-12) << gamma << " " << kBT << " " << pfd.order() << " " << tau << " " << positive_tau;
		EXPECT_NEAR(a0.imag(), -a1.imag(), 1E-12) << gamma << " " << kBT << " " << pfd.order() << " " << tau << " " << positive_tau;
	}
}

TEST_F(CorrelationFunctionTemperatureCalculatorTest, DISABLED_SimpleLorentzian)
{
	test_simple_lorentzian(1.0, 1E-14, false, false);
	test_simple_lorentzian(1.0, 1E-14, false, true);
	test_simple_lorentzian(1.0, 0.001, false, false);
	test_simple_lorentzian(1.0, 0.001, true, false);
	test_simple_lorentzian(1.0, 1, false, true);
	test_simple_lorentzian(1.0, 1, true, true);
	test_simple_lorentzian(1.0, 100, false,true);
	test_simple_lorentzian(1.0, 100, true,true);

	const double gamma = 1.0;
	SimpleLorentzian J(1.0);
	const double kBT = 1E-14;
	CorrelationFunctionTemperatureCalculator calc(J, m_pfd_hi, kBT);
	for (int i = -10; i <= 10; ++i) {
		const double tau = i*0.5/gamma;
		const double expected = exp(-std::abs(tau)*gamma);
		const std::complex<double> actual = calc.evaluate(tau);
		EXPECT_NEAR(0.0, actual.imag(), 1E-11) << tau;
		EXPECT_NEAR(expected, actual.real(), 1E-12) << tau;
	}
}

