#include <gtest/gtest.h>
#include "math/ode/solver_runge_kutta_generic.h"
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <Eigen/Core>
#include "van_der_pol.h"

using namespace rql::math::ode;

class SolverRungeKuttaGenericTest: public testing::Test
{
};


typedef Eigen::Matrix<double,1,1> Double;

struct ODE
{
	void operator()(double t, const Double& y, Double& r)
	{
		r[0] = -y[0]*y[0];
	}
};

double y(double t, double y0)
{
	return y0 / (1 + t*y0);
}

TEST_F(SolverRungeKuttaGenericTest,1D)
{
	ODE ode;
	typedef SolverRungeKuttaGeneric<Double,false> solver_type;
	Double state;
	const double y0 = 1;
	state[0] = y0;	
	solver_type solver;
	double prev_t = 0;
	for (size_t i = 1; i < 20; ++i) {
		const double t = i*0.1;
		solver.solve(ode, state, prev_t, t);
		EXPECT_NEAR(state[0], y(t, y0), 1E-6) << t;
		prev_t = t;
	}
}

TEST_F(SolverRungeKuttaGenericTest, VanDerPol)
{
	VanDerPol vdp;
	typedef SolverRungeKuttaGeneric<VanDerPol::Vec, true> solver_type;
	solver_type solver;
	VanDerPol::Vec state;
	VanDerPol::initialConditions(state);
	SolverRungeKuttaGenericNorm<VanDerPol::Vec,true>::norm(state);
	double prev_t = 0;
	const size_t n = 10000;
	for (size_t i = 1; i <= n; ++i) {
		const double t = i*(VanDerPol::defaultTime/n);
		solver.solve(vdp, state, prev_t, t);
		prev_t = t;
	}
	EXPECT_NEAR(VanDerPol::defaultResult, state[0], 1E-5);
}

