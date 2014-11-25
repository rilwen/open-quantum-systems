#include <gtest/gtest.h>
#include <math/stats/monte_carlo_result.h>
#include <math/sf/NormalDistribution.h>

using namespace rql::math::stats;

class MonteCarloResultTest: public testing::Test
{
};

TEST_F(MonteCarloResultTest,Gaussian)
{
	const unsigned int n = 10000;
	MonteCarloResult mcr;
	mcr.update(0);
	for (unsigned int i = 1; i <= n; ++i) {
		const double p = 0.5 + (0.5*i) / static_cast<double>(n+1);
		const double x = rql::math::sf::normsinv(p);
		mcr.update(x);
		mcr.update(-x);
	}
	const unsigned int nbr_samples = 2*n + 1;
	const double expected_error = sqrt(1.0 / n);
	//std::cout << "Error: " << mcr.error() << std::endl;
	ASSERT_NEAR(expected_error, mcr.error(), 3E-3);
	ASSERT_NEAR(0.0, mcr.value(), expected_error);
}
