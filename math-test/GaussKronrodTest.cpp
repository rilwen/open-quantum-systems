#include <Eigen/Core>
#include "GaussKronrodTest.h"
#include "math/MathUtils.h"
#include "math/integr/GaussKronrod.h"

using namespace Eigen;
using namespace rql::math::integr;

double unity(double x)
{
	return 1;
}

double fun1(double x)
{
	return exp(- x);
}

double fun2(double x)
{
	return x;
}

class Function
{
	public:
		double operator()(double x) const
		{
			return exp(-x);
		}
};

const double frequency = 100;

double fun3(double x)
{
	return cos(frequency*x);
}

Vector2d funVec(double x)
{
	return Vector2d(fun1(x), fun3(x));
}

double abs(const Vector2d& v)
{
	return std::abs(v[0]) + std::abs(v[1]);
}

TEST_F(GaussKronrodTest, GaussKronrod21)
{
	GaussKronrod21 igk;
	double result = igk.integrate(fun1, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, 2E-16);
	Function fun1obj;
	result = igk.integrate<Function>(fun1obj, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, 2E-16);
	result = igk.integrate(fun2, 0, 1);
	EXPECT_NEAR(0.5, result, 2E-16);
	result = igk.integrate(unity, -1, 1);
	EXPECT_NEAR(2, result, 1E-15);

	Vector2d resVec(0.0, 0.0);
	igk.integrate(funVec, 0, 1, resVec);
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], 2E-16);

	resVec *= 0.0;
	Vector2d errVec(0.0, 0.0);
	igk.integrate_with_error(funVec, 0, 1, resVec, errVec);
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], 2E-16);

	std::pair<double,double> resWErr = igk.integrate_with_error(fun2, 0, 1);
	EXPECT_NEAR(0.5, resWErr.first, 2E-16);
	EXPECT_NEAR(0, resWErr.second, 1E-1);
}

TEST_F(GaussKronrodTest, GaussKronrod41)
{
	GaussKronrod41 igk;
	double result = igk.integrate(fun1, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, 2E-16);
	Function fun1obj;
	result = igk.integrate(fun1obj, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, 2E-16);
	result = igk.integrate(fun2, 0, 1);
	EXPECT_NEAR(0.5, result, 1E-16);
	result = igk.integrate(fun3, 0, 1);
	EXPECT_NEAR(sin(frequency) / frequency, result, 1E-6);
	result = igk.integrate(unity, -1, 1);
	EXPECT_NEAR(2, result, 1E-15);

	Vector2d resVec(0.0, 0.0);
	igk.integrate(funVec, 0, 1, resVec);
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], 2E-16);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], 1E-6);

	resVec *= 0.0;
	Vector2d errVec(0.0, 0.0);
	igk.integrate_with_error(funVec, 0, 1, resVec, errVec);
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], 2E-16);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], 1E-6);

	std::pair<double,double> resWErr = igk.integrate_with_error(fun2, 0, 1);
	EXPECT_NEAR(0.5, resWErr.first, 2E-16);
	EXPECT_NEAR(0, resWErr.second, 2E-16);
}

TEST_F(GaussKronrodTest, GaussKronrod61)
{
	GaussKronrod61 igk;
	double result = igk.integrate(fun1, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, 2E-16);
	Function fun1obj;
	result = igk.integrate(fun1obj, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, 2E-16);
	result = igk.integrate(fun2, 0, 1);
	EXPECT_NEAR(0.5, result, 1E-16);
	result = igk.integrate(fun3, 0, 1);
	EXPECT_NEAR(sin(frequency) / frequency, result, 1E-15);
	result = igk.integrate(unity, -1, 1);
	EXPECT_NEAR(2, result, 1E-15);

	Vector2d resVec(0.0, 0.0);
	igk.integrate(funVec, 0, 1, resVec);
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], 2E-16);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], 1E-15);

	resVec *= 0.0;
	Vector2d errVec(0.0, 0.0);
	igk.integrate_with_error(funVec, 0, 1, resVec, errVec);
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], 2E-16);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], 1E-15);
	EXPECT_NEAR(0, abs(errVec), 1E-2);

	std::pair<double,double> resWErr = igk.integrate_with_error(fun2, 0, 1);
	EXPECT_NEAR(0.5, resWErr.first, 1E-16);
	EXPECT_NEAR(0, resWErr.second, 2E-16);
}
