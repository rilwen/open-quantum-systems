#include "AccumulatorTest.h"
#include "math/Accumulator.h"

using namespace rql::math;

TEST_F(AccumulatorTest,Double)
{
	double sum = 0;
	Accumulator<double> accumulator;

	double x = 1;
	int sign = 1;
	for (int i = 0; i < 24; ++i)
	{
		sum += x;
		accumulator += sign + x;
		x *= 0.01;
		sign *= -1;
	}

	EXPECT_TRUE(accumulator.sum() == accumulator);
	EXPECT_TRUE(accumulator == accumulator.sum());
	EXPECT_FALSE(sum == accumulator.sum());
	EXPECT_NEAR(sum, accumulator.sum(), 1E-6);
}
