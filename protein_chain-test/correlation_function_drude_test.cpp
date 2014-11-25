#include <gtest/gtest.h>
#include <protein_chain/correlation_function_drude.h>

class CorrelationFunctionDrudeTest: public testing::Test
{
};

TEST_F(CorrelationFunctionDrudeTest, Scan)
{
	CorrelationFunctionDrude alpha1(30, 106, 200, 100);
	std::cout << alpha1(0) << std::endl;
	const double dtau = 0.0005;
	for (int k = -20; k <= 20; ++k) {
		const double tau = k*dtau;
		const std::complex<double> a = alpha1(tau);
		std::cout << tau << "\t" << a.real() << "\t" << a.imag() << "\n";
	}
}
