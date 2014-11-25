#include <gtest/gtest.h>
#include <math/stats/covariance_decomposition_calculator.h>
#include <math/stats/covariance_decomposition_calculator_iterative.h>
#include <math/stats/covariance_decomposition_calculator_one_shot.h>

using namespace rql::math::stats;

class CovarianceDecompositionCalculatorTest: public testing::Test
{
protected:
	void test_calculator(const CovarianceDecompositionCalculator& calc, const Eigen::MatrixXd& covariance, Eigen::MatrixXd& decomposition);
	void test_calculator(const CovarianceDecompositionCalculator& calc, Eigen::MatrixXd& decomposition);
};

void CovarianceDecompositionCalculatorTest::test_calculator(const CovarianceDecompositionCalculator& calc, Eigen::MatrixXd& decomposition)
{
	Eigen::MatrixXd covariance(3, 3);
	covariance(0,0) = 1;
	covariance(0,1) = -0.2;
	covariance(0,2) = 0.15;
	covariance(1,1) = 0.9;
	covariance(1,2) = 0.3;
	covariance(2,2) = 1.3;
	for (size_t i = 0; i < 3; ++i)
		for (size_t j = 0; j < i; ++j)
			covariance(i, j) = covariance(j, i);
	test_calculator(calc, covariance, decomposition);
}

void CovarianceDecompositionCalculatorTest::test_calculator(const CovarianceDecompositionCalculator& calc, const Eigen::MatrixXd& covariance, Eigen::MatrixXd& decomposition)
{
	calc.decompose(covariance, decomposition);
	ASSERT_EQ(3, decomposition.rows());
	ASSERT_EQ(3, decomposition.cols());
	Eigen::MatrixXd actual_covariance(decomposition * decomposition.transpose());
	const double diff_norm = (actual_covariance - covariance).norm();
	ASSERT_TRUE(diff_norm < 2E-15) << diff_norm;
}

TEST_F(CovarianceDecompositionCalculatorTest, one_shot)
{
	CovarianceDecompositionCalculatorOneShot one_shot_calc;
	Eigen::MatrixXd decomposition;
	test_calculator(one_shot_calc, decomposition);
}

TEST_F(CovarianceDecompositionCalculatorTest, iterative)
{
	CovarianceDecompositionCalculatorIterative iterative_calc;
	Eigen::MatrixXd transform;
	test_calculator(iterative_calc, transform);
	for (int i = 0; i < transform.rows(); ++i)
		for (int k = 0; k < transform.cols(); ++k)
			if (k > i)
				ASSERT_EQ(0.0, transform(i, k)) << i << " " << k;
}
