#include <gtest/gtest.h>

class DiagonalizerRealHermitian;

class DiagonalizerRealHermitianTest: public testing::Test
{
protected:	
	void test(DiagonalizerRealHermitian& diagonalizer);
	static unsigned int dim() { return 5; }
};
