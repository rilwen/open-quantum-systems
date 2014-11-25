#include <gtest/gtest.h>
#include <boost/make_shared.hpp>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/hamiltonian_factory.h>
#include <qsd/strunz_simulator_reduced_stochastic_linear.h>

class StrunzSimulatorReducedStochasticLinearTest: public testing::Test
{
};

TEST_F(StrunzSimulatorReducedStochasticLinearTest,Regression)
{
	Eigen::VectorXd omega0(7);
	Eigen::VectorXd hwhm(7);
	Eigen::VectorXd height(7);
	omega0[0] = 0.16; hwhm[0] = 0.06; height[0] = 0.5;
	omega0[1] = 0.43; hwhm[1] = 0.05; height[1] = 0.43;
	omega0[2] = 0.86; hwhm[2] = 0.05; height[2] = 0.7;
	omega0[3] = 1.21; hwhm[3] = 0.05; height[3] = 0.7;
	omega0[4] = 1.38; hwhm[4] = 0.06; height[4] = 1.9;
	omega0[5] = 1.56; hwhm[5] = 0.07; height[5] = 1.6;
	omega0[6] = 0.05; hwhm[6] = 0.04; height[6] = 0.12;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(omega0, hwhm, height);
	const size_t N = 4;
	std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> > alphas(N, alpha);
	Eigen::MatrixXcd hamiltonian(N, N);
	HamiltonianFactory::buildChain(std::vector<double>(N, 0.0), -1.0, 0, hamiltonian);
	const double dt = 0.01;
	const size_t nbr_steps = 50;
	StrunzSimulatorReducedStochasticLinear simulator(alphas, hamiltonian, dt, nbr_steps);

	Eigen::VectorXcd initial(N);
	initial.setZero();
	initial[0] = 1;

	StrunzSimulatorReducedStochasticLinear::Workspace wksp(simulator.workspace());
	std::vector<Eigen::VectorXcd> trajectory;
	simulator.simulate(wksp, initial, trajectory);
	EXPECT_EQ(nbr_steps + 1, trajectory.size());
	EXPECT_EQ(initial, trajectory[0]);

	std::vector<double> expected;	
	expected.push_back(0.7869459914838699);
	expected.push_back(0.4312159323797339);
	expected.push_back(-0.105193661531058);
	expected.push_back(0.3680525704509895);
	expected.push_back(-0.08741801794028017);
	expected.push_back(-0.01940657538643234);
	expected.push_back(0.006137791049392543);
	expected.push_back(-0.01653135153995525);

	for (size_t i = 0; i < N; ++i) {
		/*std::cout << std::setprecision(16) << trajectory.back()[i].real() << std::endl;
		std::cout << std::setprecision(16) << trajectory.back()[i].imag() << std::endl;*/
		ASSERT_NEAR(expected[2*i], trajectory.back()[i].real(), 5E-8);
		ASSERT_NEAR(expected[2*i + 1], trajectory.back()[i].imag(), 5E-8);
	}
}
