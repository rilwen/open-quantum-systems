#include <math/ode/GSLODESystem.h>
#include <math/ode/GSLODEDriver.h>
#include <math/ode/GSLODESolver.h>
#include <gtest/gtest.h>
#include <cmath>
#include "van_der_pol.h"

using namespace rql::math::ode;

class GSLODETest: public testing::Test
{
};




TEST_F(GSLODETest, VanDerPol)
{
	boost::shared_ptr<IGSLODESystem> system = GSLODESystem<VanDerPol>::build(2, VanDerPol());
	GSLODEDriver driver = GSLODEDriver::y(system, gsl_odeiv2_step_rk8pd, 1E-6, 1E-6, 0.0);
	const size_t n = 101;
	std::vector<double> t(n);
	for (size_t i = 0; i < n; ++i)
		t[i] = i * (VanDerPol::defaultTime / 100.0);
	std::vector<Eigen::VectorXd> x(n);
	x[0].resize(2);
	VanDerPol::initialConditions(x[0]);
	GSLODESolver::solve(driver, t, x);
	EXPECT_NEAR(VanDerPol::defaultResult, x[100][0],1E-5);
}
