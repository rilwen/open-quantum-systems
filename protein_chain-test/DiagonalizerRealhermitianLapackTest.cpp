#include "diagonalizer_realhermitian_test.h"
#include "diagonalizer_realhermitian_lapack.h"

class DiagonalizerRealHermitianLapackTest: public DiagonalizerRealHermitianTest
{
};

TEST_F(DiagonalizerRealHermitianLapackTest,Test)
{
	DiagonalizerRealHermitianLapack diag(dim());
	test(diag);
}
