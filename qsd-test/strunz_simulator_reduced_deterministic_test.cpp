#include <gtest/gtest.h>
#include <boost/make_shared.hpp>
#include <qsd/strunz_simulator_reduced_deterministic.h>
#include <protein_chain/absorption_spectrum_calculator_time_dependent.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <protein_chain/hamiltonian_factory.h>

class StrunzSimulatorReducedDeterministicTest: public testing::Test
{
protected:
	void test_outputs(bool psi_norm_correction);
	void test_partial_simulations(bool psi_norm_correction);
};

void StrunzSimulatorReducedDeterministicTest::test_outputs(bool psi_norm_correction)
{
	const double hwhm = 0.25;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(1, hwhm, CorrelationFunctionLorentzians::height(0.64,hwhm));
	const size_t N = 2;
	std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> > alphas(N, alpha);
	Eigen::MatrixXcd hamiltonian(N, N);
	HamiltonianFactory::buildChain(std::vector<double>(N, 0.0), 0, 0, hamiltonian);
	StrunzSimulatorReducedDeterministic simulator(alphas, hamiltonian, psi_norm_correction);
	Eigen::VectorXcd initial(N);
	initial.fill(1/sqrt(static_cast<double>(N)));
	std::vector<Eigen::VectorXcd> states;
	StrunzSimulatorReducedDeterministic::Workspace wksp(simulator.workspace());
	const double dt = 0.1;
	const size_t nbrSteps = 50;
	simulator.simulate(wksp, initial, states, dt, nbrSteps);
	ASSERT_EQ(nbrSteps + 1, states.size());
	const double initial_norm = initial.norm();
	for (std::vector<Eigen::VectorXcd>::const_iterator it = states.begin(); it != states.end(); ++it) {
		ASSERT_EQ(N, it->size());
		ASSERT_TRUE(it->norm() <= initial_norm) << (*it);// print the vector contents on screen if we fail the norm test
	}

	std::vector<double> scal_prod_times(nbrSteps + 1);
	std::vector<std::complex<double> > scal_prods(nbrSteps + 1);
	for (size_t i = 0; i <= nbrSteps; ++i) {
		scal_prod_times[i] = i*dt;
		scal_prods[i] = initial.dot(states[i]);
	}

	std::vector<std::complex<double> > scal_prods_from_simulator;
	simulator.simulate(wksp, initial, scal_prods_from_simulator, dt, nbrSteps);
	ASSERT_EQ(nbrSteps + 1, scal_prods_from_simulator.size());
	
	for (size_t i = 0; i <= nbrSteps; ++i)
		ASSERT_EQ(scal_prods[i], scal_prods_from_simulator[i]) << i;
}

TEST_F(StrunzSimulatorReducedDeterministicTest, Outputs_NormCorrectionOn)
{
	test_outputs(true);
}

TEST_F(StrunzSimulatorReducedDeterministicTest, Outputs_NormCorrectionOff)
{
	test_outputs(false);
}

TEST_F(StrunzSimulatorReducedDeterministicTest, Regression)
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
	const size_t N = 3;
	std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> > alphas(N, alpha);
	Eigen::MatrixXcd hamiltonian(N, N);
	HamiltonianFactory::buildRing(std::vector<double>(N, 0.0), -1.3, hamiltonian);
	StrunzSimulatorReducedDeterministic simulator(alphas, hamiltonian, true);
	Eigen::VectorXcd initial(N);
	initial.fill(1/sqrt(static_cast<double>(N)));
	StrunzSimulatorReducedDeterministic::Workspace wksp(simulator.workspace());
	const double dt = 0.02;
	const size_t nbrSteps = 2500;
	std::vector<std::complex<double> > scal_prods;
	simulator.simulate(wksp, initial, scal_prods, dt, nbrSteps);
	ASSERT_EQ(nbrSteps + 1, scal_prods.size());

	std::vector<double> scal_prod_times(nbrSteps + 1);	
	for (size_t i = 0; i <= nbrSteps; ++i)
		scal_prod_times[i] = i*dt;

	const double w1 = 4;
	const double w0 = -4;
	const size_t nbrOmegas = 10;
	const double dw = (w1 - w0) / (nbrOmegas - 1);
	std::vector<double> sigma;
	AbsorptionSpectrumCalculatorTimeDependent abs_calc(w0, w1, dw);
	abs_calc.calculateAbsorptionSpectrum(scal_prod_times, scal_prods, sigma);
	std::vector<double> expected_sigma;

	expected_sigma.push_back(0.07970674228575431);
	expected_sigma.push_back(60.65302237133663);
	expected_sigma.push_back(5.52264690724893);
	expected_sigma.push_back(3.590675158433397);
	expected_sigma.push_back(0.5746663354258182);
	expected_sigma.push_back(0.1436565224050128);
	expected_sigma.push_back(0.1939409392475076);
	expected_sigma.push_back(0.5448375694355211);
	expected_sigma.push_back(0.2389281988493277);
	expected_sigma.push_back(0.1298693961176897);

	for (size_t i = 0; i < nbrOmegas; ++i)
		//std::cout << std::setprecision(16) << sigma[i] << std::endl;
		EXPECT_NEAR(expected_sigma[i], sigma[i], 1E-12*std::abs(expected_sigma[i])) << "w == " << w0 + i*dw;
}

void StrunzSimulatorReducedDeterministicTest::test_partial_simulations(bool psi_norm_correction)
{
	const double hwhm = 0.25;
	boost::shared_ptr<const CorrelationFunctionLorentzians> alpha = boost::make_shared<CorrelationFunctionLorentzians>(1, hwhm, CorrelationFunctionLorentzians::height(0.64,hwhm));
	const size_t N = 2;
	std::vector<boost::shared_ptr<const CorrelationFunctionDecomposable> > alphas(N, alpha);
	Eigen::MatrixXcd hamiltonian(N, N);
	HamiltonianFactory::buildChain(std::vector<double>(N, 0.0), 0, 0, hamiltonian);
	StrunzSimulatorReducedDeterministic simulator(alphas, hamiltonian, psi_norm_correction);

	const double dt = 0.01;
	const size_t nbr_steps = 100;
	std::vector<Eigen::MatrixXcd> hamiltonians;
	StrunzSimulatorReducedDeterministic::WorkspaceHamiltonian wksp_ham(simulator.workspace_hamiltonian());
	simulator.simulate_hamiltonians_times_minus_i_nonadaptive(wksp_ham, hamiltonians, dt, nbr_steps, true);
	ASSERT_EQ(4*nbr_steps + 1, hamiltonians.size());
	for (size_t i = 0; i < hamiltonians.size(); ++i) {
		ASSERT_EQ(N, hamiltonians[i].rows()) << i;
		ASSERT_EQ(N, hamiltonians[i].cols()) << i;
	}
	std::vector<Eigen::MatrixXcd> hamiltonians_not_all;
	simulator.simulate_hamiltonians_times_minus_i_nonadaptive(wksp_ham, hamiltonians_not_all, dt, nbr_steps, false);
	ASSERT_EQ(nbr_steps + 1, hamiltonians_not_all.size());
	ASSERT_EQ(hamiltonians.back(), hamiltonians_not_all.back());
	for (size_t i = 0; i < nbr_steps; ++i)
		ASSERT_EQ(hamiltonians[4*i], hamiltonians_not_all[i]) << i;	

	Eigen::VectorXcd initial(N);
	initial.fill(1/sqrt(static_cast<double>(N)));
	StrunzSimulatorReducedDeterministic::Workspace wksp_all(simulator.workspace());
	std::vector<Eigen::VectorXcd> trajectory1;
	simulator.simulate(wksp_all, initial, trajectory1, dt, nbr_steps);
	std::vector<Eigen::VectorXcd> trajectory2;
	StrunzSimulatorReducedDeterministic::WorkspaceState wksp_st(simulator.workspace_state());
	simulator.simulate_psi(wksp_st, hamiltonians, initial, trajectory2, dt, nbr_steps); // use pre-calculated hamiltonians to evolve psi

	// check if psi trajectory calculated in 2 different ways is approx. the same
	ASSERT_EQ(trajectory1.size(), trajectory2.size());
	ASSERT_EQ(nbr_steps + 1, trajectory1.size());
	for (size_t i = 0; i <= nbr_steps; ++i) {
		ASSERT_NEAR(0.0, (trajectory1[i] - trajectory2[i]).norm(), 1E-3) << i;
	}

	std::vector<Eigen::MatrixXcd> hamiltonians_adaptive;
	simulator.simulate_hamiltonians_times_minus_i_adaptive(wksp_ham, hamiltonians_adaptive, dt/2, 2*nbr_steps);
	ASSERT_EQ(2*nbr_steps + 1, hamiltonians_adaptive.size());
	std::vector<Eigen::MatrixXcd> hamiltonians_adaptive2;
	hamiltonians_adaptive2.resize(4*nbr_steps);
	for (size_t i = 0; i < nbr_steps; ++i) {
		hamiltonians_adaptive2[4*i] = hamiltonians_adaptive[2*i];
		hamiltonians_adaptive2[4*i + 1] = hamiltonians_adaptive[2*i + 1];
		hamiltonians_adaptive2[4*i + 2] = hamiltonians_adaptive[2*i + 1];
		hamiltonians_adaptive2[4*i + 3] = hamiltonians_adaptive[2*i + 2];
	}
	std::vector<Eigen::VectorXcd> trajectory3;
	simulator.simulate_psi(wksp_st, hamiltonians_adaptive2, initial, trajectory3, dt, nbr_steps); // use pre-calculated hamiltonians to evolve psi
	for (size_t i = 0; i <= nbr_steps; ++i) {
		ASSERT_NEAR(0.0, (trajectory1[i] - trajectory3[i]).norm(), 2E-3) << i;
	}
}

TEST_F(StrunzSimulatorReducedDeterministicTest, PartialSimulations_NormCorrectionOn)
{
	test_partial_simulations(true);	
}

TEST_F(StrunzSimulatorReducedDeterministicTest, PartialSimulations_NormCorrectionOff)
{
	test_partial_simulations(false);
}
