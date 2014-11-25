#include "diagonalizer_realhermitian_test.h"
#include "diagonalizer_realhermitian_eigen.h"

class DiagonalizerRealHermitianEigenTest: public DiagonalizerRealHermitianTest
{
};

TEST_F(DiagonalizerRealHermitianEigenTest,Test)
{
	DiagonalizerRealHermitianEigen diag(dim());
	test(diag);
}
