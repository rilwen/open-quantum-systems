#include <gtest/gtest.h>
#include <boost/make_shared.hpp>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/hamiltonian_factory.h>
#include <qsd/strunz_simulator_reduced_stochastic_nonlinear.h>

class StrunzSimulatorReducedStochasticNonlinearTest: public testing::Test
{
};

TEST_F(StrunzSimulatorReducedStochasticNonlinearTest,Regression)
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
	StrunzSimulatorReducedStochasticNonlinear simulator(alphas, hamiltonian, dt, nbr_steps);

	Eigen::VectorXcd initial(N);
	initial.setZero();
	initial[0] = 1;

	StrunzSimulatorReducedStochasticNonlinear::Workspace wksp(simulator.workspace());
	std::vector<Eigen::VectorXcd> trajectory;
	simulator.simulate(wksp, initial, trajectory);
	EXPECT_EQ(nbr_steps + 1, trajectory.size());
	EXPECT_EQ(initial, trajectory[0]);
	for (size_t i = 0; i < nbr_steps; ++i) {
		ASSERT_NEAR(1.0, trajectory[i+1].norm(), 1E-7) << i << " " << trajectory[i+1];
	}

	std::vector<double> expected;
	expected.push_back(0.7857467335385239);
	expected.push_back(0.4851328080063778);
	expected.push_back(-0.112505934494126);
	expected.push_back(0.3560416043255957);
	expected.push_back(-0.08437996314761742);
	expected.push_back(-0.02021385288987206);
	expected.push_back(0.006216117389288343);
	expected.push_back(-0.01605971269794333);

	for (size_t i = 0; i < N; ++i) {
		/*std::cout << std::setprecision(16) << trajectory.back()[i].real() << std::endl;
		std::cout << std::setprecision(16) << trajectory.back()[i].imag() << std::endl;*/
		ASSERT_NEAR(expected[2*i], trajectory.back()[i].real(), 5E-8);
		ASSERT_NEAR(expected[2*i + 1], trajectory.back()[i].imag(), 5E-8);
	}

	std::vector<double> expected_occupation_probabilities;
	expected_occupation_probabilities.push_back(0.9999010899337275);
	expected_occupation_probabilities.push_back(0.9966277840987565);
	expected_occupation_probabilities.push_back(0.9892701471581221);
	expected_occupation_probabilities.push_back(0.9784924633402122);
	expected_occupation_probabilities.push_back(0.9649355007763267);
	expected_occupation_probabilities.push_back(0.9489967205959219);
	expected_occupation_probabilities.push_back(0.9316330433060124);
	expected_occupation_probabilities.push_back(0.9125430927949354);
	expected_occupation_probabilities.push_back(0.8921505318858678);
	expected_occupation_probabilities.push_back(0.8707946979925719);

	for (size_t i = 0; i < nbr_steps; i += 5) {
		//std::cout << std::setprecision(16) << pow(std::abs(trajectory[i+1][0]), 2) << std::endl;
		ASSERT_NEAR(expected_occupation_probabilities[i/5], pow(std::abs(trajectory[i+1][0]), 2), 1E-7);
	}

	// test that simulations are random
	std::vector<Eigen::VectorXcd> trajectory2;
	simulator.simulate(wksp, initial, trajectory2);
	EXPECT_EQ(nbr_steps + 1, trajectory2.size());
	EXPECT_EQ(initial, trajectory2[0]);
	for (size_t i = 0; i < nbr_steps; ++i) {
		ASSERT_NEAR(1.0, trajectory2[i+1].norm(), 1E-7) << i << " " << trajectory2[i+1];
		ASSERT_NE(trajectory[i+1], trajectory2[i+1]);
	}
}
