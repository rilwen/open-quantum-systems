#include <gtest/gtest.h>
#include <math/linalg/MatrixFunctions.h>

using namespace rql::math::linalg;

class MatrixFunctionsTest: public testing::Test
{
};

static double id(double x)
{
	return x;
}

TEST_F(MatrixFunctionsTest, Id)
{
	const int dim = 4;
	Eigen::MatrixXd A(4, 4);
	for (size_t r = 0; r < dim; ++r) {
		for (size_t c = 0; c < dim; ++c) {
			A(r,c) = r + c;
		}
	}
	Eigen::MatrixXd result;
	MatrixFunctions::hermitian_matrix_function(A, id, result);
	ASSERT_NEAR(0.0, (A - result).norm(), 5E-14);
	MatrixFunctions::hermitian_matrix_function(result, id, result);
	ASSERT_NEAR(0.0, (A - result).norm(), 5E-14);
}
