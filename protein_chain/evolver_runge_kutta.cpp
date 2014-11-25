#include "evolver_runge_kutta.h"

#include <cassert>
#include <cmath>

EvolverRungeKutta::EvolverRungeKutta(const Eigen::MatrixXcd& hamiltonian, double timeStep)
	: EvolverStepping(timeStep), m_M(hamiltonian*std::complex<double>(0.0,-1.0))
{
	assert(m_M.rows() == m_M.cols());
	for (unsigned int k = 0; k < m_order; ++k) {
		m_k[k].resize(m_M.rows());
	}
}

// Coefficients of the 4th order RK method for linear ODEs from "Runge-Kutta method for Linear Ordinary Differential Equations", D. W. Zingg and T. T. Chisholm
const boost::array<double,EvolverRungeKutta::m_order> EvolverRungeKutta::m_b = { 0.07801567728325, 0.04708870117112, 0.47982272993855, 0.39507289160708 };
const boost::array<double,EvolverRungeKutta::m_order> EvolverRungeKutta::m_c = { 0, 0.69631521002413, 0.29441651741, 0.82502163765 };
const boost::array<double,1> EvolverRungeKutta::m_a1 = { m_c[0] };
const boost::array<double,2> EvolverRungeKutta::m_a2 = { m_b[0], 0.21640084013679 };
const boost::array<double,3> EvolverRungeKutta::m_a3 = { m_b[0], m_b[1], 0.69991725920066 };

void EvolverRungeKutta::step(const double step_size, Eigen::VectorXcd& state)
{
	/*
	// Euler:
	state += step_size*m_M*state;
	*/
	m_k[0].noalias() = m_M*state;
	m_k[1].noalias() = m_M*(state + step_size*m_a1[0]*m_k[0]);
	m_k[2].noalias() = m_M*(state + step_size*(m_a2[0]*m_k[0] + m_a2[1]*m_k[1]));
	m_k[3].noalias() = m_M*(state + step_size*(m_a3[0]*m_k[0] + m_a3[1]*m_k[1] + m_a3[2]*m_k[2]));
	for (unsigned int i = 0; i < m_order; ++i) {
		state.noalias() += (step_size*m_b[i])*m_k[i];
	}
}

void EvolverRungeKutta::step(const double step_size, Eigen::MatrixXcd& states)
{
	Eigen::VectorXcd state(states.rows());
	for (unsigned int i = 0; i < static_cast<unsigned int>(states.cols()); ++i) {
		state = states.col(i);
		step(step_size, state);
		states.col(i) = state;
	}
}


