#include <gtest/gtest.h>
#include <hamiltonian_factory.h>
#include <protein_chain/absorption_spectrum_calculator_time_dependent.h>
#include <protein_chain/correlation_function_lorentzians.h>
#include <qsd/strunz_simulator_grid.h>
#include <cmath>
#include <iomanip>
#include <iostream>

class StrunzSimulatorGridTest: public testing::Test
{
};

//! Never made it work
//TEST_F(StrunzSimulatorTest,GridStochastic)
//{
//	boost::shared_ptr<CorrelationFunction> singleAlpha = boost::shared_ptr<CorrelationFunction>(new CorrelationFunctionLorentzians(0, 0.25, 0.5));
//	const size_t nbrSites = 2;
//	Eigen::MatrixXcd hamiltonian(nbrSites, nbrSites);
//	std::vector<boost::shared_ptr<const CorrelationFunction> > alpha(nbrSites, singleAlpha);
//	HamiltonianFactory::buildChain(std::vector<double>(nbrSites, 0.0), -1, 0, hamiltonian);
//	const size_t nbr_steps = 20;
//	StrunzSimulatorGrid simulator(alpha, hamiltonian, 0.01, nbr_steps, 1, true);
//	Eigen::VectorXcd initial(2);
//	initial[0] = 1;
//	initial[1] = 0;
//	std::vector<Eigen::VectorXcd> trajectory;
//	StrunzSimulatorGrid::Workspace wksp(simulator.workspace());
//	simulator.simulateTrajectory(wksp, initial, trajectory);
//	const double tol = 1E-15;
//	EXPECT_EQ(nbr_steps + 1, trajectory.size());
//	EXPECT_EQ(initial, trajectory[0]);
//	EXPECT_NEAR(1.024355909405471, trajectory.back()[0].real(), tol);
//	EXPECT_NEAR(-0.08215468589202808, trajectory.back()[0].imag(), tol);
//	EXPECT_NEAR(0.009627224891992399, trajectory.back()[1].real(), tol);
//	EXPECT_NEAR(0.2100467452367024, trajectory.back()[1].imag(), tol);
//}

TEST_F(StrunzSimulatorGridTest,GridCommon)
{
	boost::shared_ptr<CorrelationFunction> singleAlpha = boost::shared_ptr<CorrelationFunction>(new CorrelationFunctionLorentzians(0, 0.25, 0.5));
	const size_t nbrSites = 2;
	Eigen::MatrixXcd hamiltonian(nbrSites, nbrSites);
	HamiltonianFactory::buildChain(std::vector<double>(nbrSites, 0.0), -1, 0, hamiltonian);
	StrunzSimulatorGrid simulator(singleAlpha, nbrSites, hamiltonian, 0.01, 20, 1, false);
	Eigen::VectorXcd initial(2);
	initial[0] = 1;
	initial[1] = 0;
	std::vector<Eigen::VectorXcd> trajectory;
	StrunzSimulatorGrid::Workspace wksp(simulator.workspace());
	simulator.simulateTrajectory(wksp, initial, trajectory);
	const double tol = 1E-15;
	EXPECT_NEAR(0.97355971494807791, trajectory.back()[0].real(), tol);
	EXPECT_NEAR(0.0, trajectory.back()[0].imag(), tol);
	EXPECT_NEAR(0.0, trajectory.back()[1].real(), tol);
	EXPECT_NEAR(0.19775911747980654, trajectory.back()[1].imag(), tol);
}

TEST_F(StrunzSimulatorGridTest,GridNonStochastic)
{
	boost::shared_ptr<CorrelationFunction> singleAlpha = boost::shared_ptr<CorrelationFunction>(new CorrelationFunctionLorentzians(0, 0.25, 0.5));
	const size_t nbrSites = 2;
	Eigen::MatrixXcd hamiltonian(nbrSites, nbrSites);
	HamiltonianFactory::buildChain(std::vector<double>(nbrSites, 0.0), -1, 0, hamiltonian);
	StrunzSimulatorGrid simulator(singleAlpha, nbrSites, hamiltonian, 0.01, 20, 1, false);
	Eigen::VectorXcd initial(2);
	initial[0] = 1;
	initial[1] = 0;
	std::vector<Eigen::VectorXcd> trajectory;
	StrunzSimulatorGrid::Workspace wksp(simulator.workspace());	
	simulator.simulateTrajectory(wksp, initial, trajectory);
	EXPECT_NEAR(0.97355971494807791, trajectory.back()[0].real(), 1E-15);
	EXPECT_NEAR(0.0, trajectory.back()[0].imag(), 1E-15);
	EXPECT_NEAR(0.0, trajectory.back()[1].real(), 1E-15);
	EXPECT_NEAR(0.19775911747980654, trajectory.back()[1].imag(), 1E-15);
}


