#include "DomainConverterTest.h"
#include "math/integr/DomainConverter.h"
#include <Eigen/Core>
#include <cmath>
#include <limits>

using namespace rql::math::integr;
using namespace Eigen;

double f1(double x)
{
	return exp(-x*x);
}

double f2(double x)
{
	return x;
}

Vector2d f12(double x)
{
	Vector2d v;
	v[0] = f1(x);
	v[1] = f2(x);
	return v;
}

TEST_F(DomainConverterTest, Scalar)
{
	const DomainConverterObj<double(&)(double),double> dc1(f1);
	EXPECT_NEAR(1, dc1(0), 1E-16);
	EXPECT_EQ(0, dc1(-1));
	EXPECT_EQ(0, dc1(1));
	EXPECT_NEAR(dc1(-0.5), dc1(0.5), 1E-16);

	const DomainConverterObj<double(&)(double),double> dc2(f2);
	EXPECT_NEAR(0, dc2(0), 1E-16);
	EXPECT_NEAR(-dc2(-0.5), dc2(0.5), 1E-16);
}

TEST_F(DomainConverterTest, Vector)
{
	const DomainConverterObj<Vector2d(&)(double),Vector2d> dc(f12);
	Vector2d f = dc(0);
	EXPECT_NEAR(1, f[0], 1E-16);
	EXPECT_NEAR(0, f[1], 1E-16);
}

TEST_F(DomainConverterTest, Limits)
{
	const double y = -4;
	const double x = DomainConverter::to_finite_range(y);
	EXPECT_NEAR(y, DomainConverter::from_finite_range(x), 1E-15);
	EXPECT_NEAR(0, DomainConverter::from_finite_range(0), 1E-16);
	EXPECT_NEAR(0, DomainConverter::to_finite_range(0), 1E-16);
	EXPECT_NEAR(DomainConverter::from_finite_range(x), -DomainConverter::from_finite_range(-x), 1E-16);
	EXPECT_NEAR(DomainConverter::to_finite_range(y), -DomainConverter::to_finite_range(-y), 1E-16);
	EXPECT_NEAR(DomainConverter::lower_finite_range(), DomainConverter::to_finite_range(-std::numeric_limits<double>::infinity()), 1E-16);
	EXPECT_NEAR(DomainConverter::upper_finite_range(), DomainConverter::to_finite_range(std::numeric_limits<double>::infinity()), 1E-16);
	EXPECT_EQ(-std::numeric_limits<double>::infinity(), DomainConverter::from_finite_range(DomainConverter::lower_finite_range()));
	EXPECT_EQ(std::numeric_limits<double>::infinity(), DomainConverter::from_finite_range(DomainConverter::upper_finite_range()));
}
