#include "AdaptiveGKTest.h"
#include <Eigen/Core>
#include "math/MathUtils.h"
#include "AdaptiveGKTest.h"
#include "math/integr/AdaptiveGK.h"

using namespace Eigen;
using namespace rql::math::integr;

static double unity(double x)
{
	return 1;
}

static double fun1(double x)
{
	return exp(- x);
}

static double fun2(double x)
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

static const double frequency = 100;

static double fun3(double x)
{
	return cos(frequency*x);
}

static Vector2d funVec(double x)
{
	return Vector2d(fun1(x), fun3(x));
}

static double abs(const Vector2d& v)
{
	return std::abs(v[0]) + std::abs(v[1]);
}

struct Vector2dNorm
{
	double operator()(const Vector2d& v) const
	{
		return abs(v);
	}
};

TEST_F(AdaptiveGKTest, Adaptive15)
{
	const double tolerance = 1E-14;
	AdaptiveGK15 integrator(20, tolerance, false);

	double result = integrator.integrate(fun2, -1, 1);
	EXPECT_NEAR(0, result, tolerance);

	result = integrator.integrate(fun1, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, tolerance);

	Function fun1obj;
	result = integrator.integrate(fun1obj, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, tolerance);

	result = integrator.integrate(fun3, 0, 1);
	EXPECT_NEAR(sin(frequency) / frequency, result, tolerance);
	const double lInf_err = std::abs(sin(frequency) / frequency - result);

	Vector2d resVec;
	integrator.integrate(funVec, 0, 1, resVec, Vector2dNorm());
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], tolerance);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], tolerance);


	integrator.integrate(funVec, 0, 1, resVec, Vector2dNorm());
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], tolerance);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], tolerance);
	const double l1_err = std::abs(sin(frequency) / frequency - result);
	EXPECT_TRUE(l1_err <= lInf_err);
}

//TEST_F(AdaptiveGKTest, Adaptive21)
//{
//	const double tolerance = 1E-14;
//	AdaptiveGK21 integrator(20, tolerance, false);
//
//	double result = integrator.integrate(fun2, -1, 1);
//	EXPECT_NEAR(0, result, tolerance);
//
//	result = integrator.integrate(fun1, 0, 1);
//	EXPECT_NEAR(1 - exp(-1.0), result, tolerance);
//
//	Function fun1obj;
//	result = integrator.integrate<Function,double>(fun1obj, 0, 1, std::abs<double>);
//	EXPECT_NEAR(1 - exp(-1.0), result, tolerance);
//
//	result = integrator.integrate(fun3, 0, 1);
//	EXPECT_NEAR(sin(frequency) / frequency, result, tolerance);
//	const double lInf_err = std::abs(sin(frequency) / frequency - result);
//
//	Vector2d resVec = integrator.integrate(funVec, 0, 1, Vector2dNorm());
//	EXPECT_NEAR(1 - exp(-1.0), resVec[0], tolerance);
//	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], tolerance);
//
//
//	resVec = integrator.integrate(funVec, 0, 1, Vector2dNorm());
//	EXPECT_NEAR(1 - exp(-1.0), resVec[0], tolerance);
//	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], tolerance);
//	const double l1_err = std::abs(sin(frequency) / frequency - result);
//	EXPECT_TRUE(l1_err <= lInf_err);
//}

TEST_F(AdaptiveGKTest, Adaptive41)
{
	const double tolerance = 1E-14;
	AdaptiveGK41 integrator(20, tolerance, false);

	double result = integrator.integrate(fun2, -1, 1);
	EXPECT_NEAR(0, result, tolerance);

	result = integrator.integrate(fun1, 0, 1);
	EXPECT_NEAR(1 - exp(-1.0), result, tolerance);

	Function fun1obj;
	result = integrator.integrate(fun1obj, 0., 1.);
	EXPECT_NEAR(1 - exp(-1.0), result, tolerance);

	result = integrator.integrate(fun3, 0, 1);
	EXPECT_NEAR(sin(frequency) / frequency, result, tolerance);
	const double lInf_err = std::abs(sin(frequency) / frequency - result);

	Vector2d resVec;
	integrator.integrate(funVec, 0, 1, resVec, Vector2dNorm());
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], tolerance);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], tolerance);


	integrator.integrate(funVec, 0, 1, resVec, Vector2dNorm());
	EXPECT_NEAR(1 - exp(-1.0), resVec[0], tolerance);
	EXPECT_NEAR(sin(frequency) / frequency, resVec[1], tolerance);
	const double l1_err = std::abs(sin(frequency) / frequency - result);
	EXPECT_TRUE(l1_err <= lInf_err);
}
