#include <gtest/gtest.h>
#include <complex>
#include <Eigen/Core>
#include <iostream>
#include "evolver_taylor_expansion.h"

using namespace Eigen;

class EvolverTaylorExpansionTest: public testing::Test
{
};

static const unsigned int ORDER = 4;
static const double TIME_STEP = 0.01;
static const double TIME = 10;

TEST_F(EvolverTaylorExpansionTest, Hermitian)
{
	const unsigned int dim = 4;
	MatrixXcd ham(dim,dim);
	ham.fill(0.0);
	ham(0,0) = -0.01;
	ham(1,1) = 0.1;
	ham(2,2) = -1;
	ham(3,3) = 10;
	std::vector<double> times(1);
	times[0] = TIME;
	VectorXcd initState(dim);
	for (unsigned int i = 0; i < dim; ++i) {
		initState[i] = std::complex<double>(1.0, 0.0);
	}
	std::vector<VectorXcd> states(1);
	states[0].resize(dim);
	EvolverTaylorExpansion<ORDER> ev(ham, TIME_STEP);
	ev.evolve(initState, times, states);
	VectorXcd expectedFinalState(dim);
	const double tolerances[] = {1E-14, 1E-14, 1E-9, 1E-4};
	for (unsigned int i = 0; i < dim; ++i) {
		expectedFinalState[i] = initState[i] * exp(- times[0]*std::complex<double>(0.0,1.0)*ham(i,i));
		//std::cout << expectedFinalState[i] << "\t" << states[0][i] << "\t" << std::abs(expectedFinalState[i] - states[0][i]) << "\n";
		EXPECT_NEAR(0.0, std::abs(expectedFinalState[i] - states[0][i]), tolerances[i]);
	}
	std::cout << std::endl;
}

TEST_F(EvolverTaylorExpansionTest, AntiHermitian)
{
	const unsigned int dim = 4;
	MatrixXcd ham(MatrixXcd::Zero(dim,dim));
	ham(0,0) = -0.01;
	ham(1,1) = -0.1;
	ham(2,2) = -1;
	ham(3,3) = -10;
	ham *= std::complex<double>(0.0, 1.0);
	std::vector<double> times(1);
	times[0] = TIME;
	VectorXcd initState(dim);
	for (unsigned int i = 0; i < dim; ++i) {
		initState[i] = std::complex<double>(1.0, 0.0);
	}
	std::vector<VectorXcd> states(1);
	states[0].resize(dim);
	EvolverTaylorExpansion<ORDER> ev(ham, TIME_STEP);
	ev.evolve(initState, times, states);
	VectorXcd expectedFinalState(dim);
	const double tolerances[] = {1E-14, 2E-14, 1E-9, 1E-4};
	for (unsigned int i = 0; i < dim; ++i) {
		expectedFinalState[i] = initState[i] * exp(-times[0]*std::complex<double>(0.0,1.0)*ham(i,i));
		const double rel_error = std::abs(expectedFinalState[i] - states[0][i]) / std::abs(expectedFinalState[i]);
		//std::cout << expectedFinalState[i] << "\t" << states[0][i] << "\t" << rel_error << "\n";
		EXPECT_NEAR(0.0, rel_error, tolerances[i]) << i;
	}
	std::cout << std::endl;
}

TEST_F(EvolverTaylorExpansionTest, AlmostHermitian)
{
	const unsigned int dim = 4;
	MatrixXcd ham(MatrixXcd::Zero(dim,dim));
	ham(0,0) = 0.01;
	ham(1,1) = 0.1;
	ham(2,2) = 1;
	ham(3,3) = 10;	
	for (unsigned int i = 0; i < dim; ++i) {
		ham(i,i).imag(-0.1);
	}
	std::vector<double> times(1);
	times[0] = TIME;
	VectorXcd initState(dim);
	for (unsigned int i = 0; i < dim; ++i) {
		initState[i] = std::complex<double>(1.0, 0.0);
	}
	std::vector<VectorXcd> states(1);
	states[0].resize(dim);
	EvolverTaylorExpansion<ORDER> ev(ham, TIME_STEP);
	ev.evolve(initState, times, states);
	VectorXcd expectedFinalState(dim);
	const double tolerances[] = {1E-13, 2E-13, 1E-9, 1E-4};
	for (unsigned int i = 0; i < dim; ++i) {
		expectedFinalState[i] = initState[i] * exp(-times[0]*std::complex<double>(0.0,1.0)*ham(i,i));
		const double rel_error = std::abs(expectedFinalState[i] - states[0][i]) / std::abs(expectedFinalState[i]);
		//std::cout << expectedFinalState[i] << "\t" << states[0][i] << "\t" << rel_error << "\n";
		EXPECT_NEAR(0.0, rel_error, tolerances[i]) << i;
	}
	std::cout << std::endl;
}

TEST_F(EvolverTaylorExpansionTest, NonDiagonal)
{
	const unsigned int dim = 2;
	MatrixXcd ham(dim,dim);
	ham(0,0) = 0;
	ham(0,1) = 1;
	ham(1,0) = 1;
	ham(1,1) = 1;
	MatrixXcd U(dim,dim);
	const double t = 0.98;
	U(0,0) = std::complex<double>(0.590779121942, 0.135598126886);
	U(0,1) = std::complex<double>(-0.374315718684, -0.701769848176);
	U(1,0) = std::complex<double>(-0.374315718684, -0.701769848176);
	U(1,1) = std::complex<double>(0.216463403258, -0.566171721289);
	std::vector<double> times(1);
	times[0] = t;
	VectorXcd initState(dim);
	initState[0] = 0.1;
	initState[1] = -0.2;
	const VectorXcd expectedFinalState(U*initState);
	EvolverTaylorExpansion<ORDER> ev(ham, TIME_STEP);
	std::vector<VectorXcd> states(1);
	states[0].resize(dim);
	ev.evolve(initState, times, states);
	const double tolerances[] = {1E-9, 1E-9};
	for (unsigned int i = 0; i < dim; ++i) {
		const double rel_error = std::abs(expectedFinalState[i] - states[0][i]) / std::abs(expectedFinalState[i]);
//		std::cout << expectedFinalState[i] << "\t" << states[0][i] << "\t" << rel_error << "\n";
		EXPECT_NEAR(0.0, rel_error, tolerances[i]) << i;
	}
}