#include <protein_chain/partial_function_decomposition_bose.h>
#include <math/MathUtils.h>
#include <gtest/gtest.h>
#include <cmath>

class PartialFunctionDecompositionBoseTest: public testing::Test
{
};

TEST_F(PartialFunctionDecompositionBoseTest,Compare)
{
	const unsigned int order = 20;
	PartialFunctionDecompositionBose pfd(order);
	ASSERT_EQ(order, pfd.order());
	for (double log_x = -5; log_x <= 5; log_x += 0.1) {
		double x = exp(log_x);
		std::complex<double> actual = pfd.evaluate(x);
		ASSERT_NEAR(0.0, actual.imag(), 1E-12) << x;
		ASSERT_TRUE(actual.real() > 0) << x;
		double expected = pfd.full_evaluate(x);
		ASSERT_NEAR(expected, pfd.full_evaluate(std::complex<double>(x)).real(), std::abs(expected)*1E-12);
		ASSERT_NEAR(0.0, pfd.full_evaluate(std::complex<double>(x)).imag(), 0.0);
		double relative_error = std::abs((actual - expected) / expected);
		ASSERT_NEAR(0.0, relative_error, 1E-2*x) << x;
		//std::cout << x << "\t" << actual.real() << "\t" << pfd.full_evaluate(x) << "\t" << actual.imag() << "\n";

		x = -x;
		actual = pfd.evaluate(x);
		ASSERT_NEAR(0.0, actual.imag(), 1E-12) << x;
		//EXPECT_TRUE(actual.real() < 0) << x << " " << actual;
		expected = pfd.full_evaluate(x);
		ASSERT_NEAR(expected, pfd.full_evaluate(std::complex<double>(x)).real(), std::abs(expected)*1E-12);
		ASSERT_NEAR(0.0, pfd.full_evaluate(std::complex<double>(x)).imag(), 0.0);
		/*relative_error = std::abs((actual - expected) / expected);
		ASSERT_NEAR(0.0, relative_error, 1E-2*std::abs(x)) << x;*/
	}
}
